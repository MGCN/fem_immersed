/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "UserForcing.h"
#include "Function.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "Assembly.h"

registerMooseObject("ImmersedBoundaryApp", UserForcing);

template<>
InputParameters validParams<UserForcing>()
{
    InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("vel_x",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
    params.addCoupledVar("vel_y",
                         "The gradient of this variable will be used as the "
                         "velocity vector.");
    params.addCoupledVar("vel_z",
                         "The gradient of this variable will be used as the "
                         "velocity vector.");
    params.addCoupledVar("pressure","pressure");
    params.addParam<Real>("viscosity", 1.0,
                          "dynamic viscosity");
    params.addRequiredCoupledVar("disp_x",
                                 "The gradient of this variable will be used as the velocity vector.");
    params.addCoupledVar("disp_y",
                         "The gradient of this variable will be used as the velocity vector.");
    params.addCoupledVar("disp_z",
                         "The gradient of this variable will be used as the velocity vector.");
    params.addRequiredParam<unsigned int>("component", "The component for the stress vector");
    params.addParam<FunctionName>("function",
                                  "The function that describes the pressure");


    return params;
}

UserForcing::UserForcing(const InputParameters & parameters) :
Kernel(parameters),
_component(getParam<unsigned int>("component")),
_pressureFluid(parameters.isParamValid("pressure") ? coupledValue("pressure"):_zero),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_y(_mesh.dimension() >= 2 ? coupledGradient("disp_y") : _grad_zero),
_grad_disp_z(_mesh.dimension() >= 3 ? coupledGradient("disp_z") : _grad_zero),
_grad_vel_x(coupledGradient("vel_x")),
_grad_vel_y(_mesh.dimension() >= 2 ? coupledGradient("vel_y") : _grad_zero),
_grad_vel_z(_mesh.dimension() >= 3 ? coupledGradient("vel_z") : _grad_zero),
_function2(getFunction("function")),
_normals(_assembly.normals()),
_mu(getParam<Real>("viscosity"))
{
}



Real
UserForcing::computeQpResidual()
{
    RealTensorValue F(_grad_disp_x[_qp], _grad_disp_y[_qp],  _grad_disp_z[_qp]);
    
    RealVectorValue normal;
    
    for (int i = 0; i < 3; ++i)
    {
        F(i, i) += 1.0;
        normal(i) = _normals[_qp](i);
    }
    
    Real J = F.det();
    
    RealTensorValue invFtr = F.inverse().transpose();
    
    RealTensorValue temp = J * invFtr;
    
    // First row
    sigma(0,0) = 2.0 * _mu*_grad_vel_x[_qp](0)  - _pressureFluid[_qp]; //* _function2.value(_t, _q_point[_qp]);
    sigma(0,1) = 1.0 * _mu*(_grad_vel_x[_qp](1) + _grad_vel_y[_qp](0));
    sigma(0,2) = 1.0 * _mu*(_grad_vel_x[_qp](2) + _grad_vel_z[_qp](0));
    
    // Second row
    sigma(1,0) = 1.0 * _mu*(_grad_vel_y[_qp](0) + _grad_vel_x[_qp](1));
    sigma(1,1) = 2.0 * _mu*_grad_vel_y[_qp](1)  - _pressureFluid[_qp]; // * _function2.value(_t, _q_point[_qp]);
    sigma(1,2) = 1.0 * _mu*(_grad_vel_y[_qp](2) + _grad_vel_z[_qp](1));
    
    // Third row
    sigma(2,0) = 1.0 * _mu*(_grad_vel_z[_qp](0) + _grad_vel_x[_qp](2));
    sigma(2,1) = 1.0 * _mu*(_grad_vel_z[_qp](1) + _grad_vel_y[_qp](2));
    sigma(2,2) = 2.0 * _mu*_grad_vel_z[_qp](2)  - _pressureFluid[_qp]; // * _function2.value(_t, _q_point[_qp]);
    
    
    RealTensorValue  P_fluid = temp * sigma;

    return - 1.0  * ( P_fluid(_component,0) * _grad_test[_i][_qp](0)
                  +   P_fluid(_component,1) * _grad_test[_i][_qp](1)
                  +   P_fluid(_component,2) * _grad_test[_i][_qp](2));
}


Real
UserForcing::computeQpJacobian()
{
    
    RealVectorValue zeros(0.0, 0.0, 0.0);
    
    RealTensorValue H(zeros, zeros, zeros);
    
    for (unsigned k = 0; k < 3; ++k)
        H(_component, k) = _grad_phi[_j][_qp](k);
    
    RealTensorValue F(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);
    RealVectorValue normal;
    
    for (int i = 0; i < 3; ++i)
    {
        F(i, i) += 1.0;
        normal(i) = _normals[_qp](i);
    }
    
    Real J = F.det();
    RealTensorValue invFtr = F.inverse().transpose();
    
    RealTensorValue TEMP =
    invFtr.contract(H) * invFtr - invFtr * H.transpose() * invFtr;
    RealTensorValue temp = J * TEMP;
    
    // First row
    sigma(0,0) = 2.0 * _mu * _grad_vel_x[_qp](0)  - _pressureFluid[_qp]; // * _function2.value(_t, _q_point[_qp]);
    sigma(0,1) = 1.0 * _mu * (_grad_vel_x[_qp](1) + _grad_vel_y[_qp](0));
    sigma(0,2) = 1.0 * _mu * (_grad_vel_x[_qp](2) + _grad_vel_z[_qp](0));
    
    // Second row
    sigma(1,0) = 1.0 * _mu * (_grad_vel_y[_qp](0) + _grad_vel_x[_qp](1));
    sigma(1,1) = 2.0 * _mu * _grad_vel_y[_qp](1)  - _pressureFluid[_qp]; // * _function2.value(_t, _q_point[_qp]);
    sigma(1,2) = 1.0 * _mu * (_grad_vel_y[_qp](2) + _grad_vel_z[_qp](1));
    
    // Third row
    sigma(2,0) = 1.0 * _mu * (_grad_vel_z[_qp](0) + _grad_vel_x[_qp](2));
    sigma(2,1) = 1.0 * _mu * (_grad_vel_z[_qp](1) + _grad_vel_y[_qp](2));
    sigma(2,2) = 2.0 * _mu * _grad_vel_z[_qp](2)  - _pressureFluid[_qp];// * _function2.value(_t, _q_point[_qp]);
    
    RealTensorValue  P_fluid = temp * sigma;
    
    return - 1.0  * (P_fluid(_component,0) * _grad_test[_i][_qp](0)
                  +  P_fluid(_component,1) * _grad_test[_i][_qp](1)
                  +  P_fluid(_component,2) * _grad_test[_i][_qp](2));
}

Real
UserForcing::computeQpOffDiagJacobian(unsigned int jvar)
{
    if (jvar == 3)
        return 0.0;
    
    RealVectorValue zeros(0.0, 0.0, 0.0);
    
    RealTensorValue H(zeros, zeros, zeros);
    
    for (unsigned k = 0; k < 3; ++k)
        H(jvar, k) = _grad_phi[_j][_qp](k);
    
    RealTensorValue F(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);
    RealVectorValue normal;
    
    for (int i = 0; i < 3; ++i)
    {
        F(i, i) += 1.0;
        normal(i) = _normals[_qp](i);
    }
    
    Real J = F.det();
    RealTensorValue invFtr = F.inverse().transpose();
    
    RealTensorValue TEMP =
    (invFtr.contract(H)) * invFtr - (invFtr * (H.transpose())) * invFtr;
    RealTensorValue temp = J * TEMP;
    
    // First row
    sigma(0,0) = 2.0 * _mu * _grad_vel_x[_qp](0)  - _pressureFluid[_qp]; // * _function2.value(_t, _q_point[_qp]);
    sigma(0,1) = 1.0 * _mu * (_grad_vel_x[_qp](1) + _grad_vel_y[_qp](0));
    sigma(0,2) = 1.0 * _mu * (_grad_vel_x[_qp](2) + _grad_vel_z[_qp](0));
    
    // Second row
    sigma(1,0) = 1.0 * _mu * (_grad_vel_y[_qp](0) + _grad_vel_x[_qp](1));
    sigma(1,1) = 2.0 * _mu * _grad_vel_y[_qp](1)  - _pressureFluid[_qp]; // * _function2.value(_t, _q_point[_qp]);
    sigma(1,2) = 1.0 * _mu * (_grad_vel_y[_qp](2) + _grad_vel_z[_qp](1));
    
    // Third row
    sigma(2,0) = 1.0 * _mu * (_grad_vel_z[_qp](0) + _grad_vel_x[_qp](2));
    sigma(2,1) = 1.0 * _mu * (_grad_vel_z[_qp](1) + _grad_vel_y[_qp](2));
    sigma(2,2) = 2.0 * _mu * _grad_vel_z[_qp](2)  - _pressureFluid[_qp];// * _function2.value(_t, _q_point[_qp]);
    
    RealTensorValue  P_fluid = temp * sigma;
    
    return - 1.0  * (P_fluid(_component,0) * _grad_test[_i][_qp](0)
                  +  P_fluid(_component,1) * _grad_test[_i][_qp](1)
                  +  P_fluid(_component,2) * _grad_test[_i][_qp](2));
}


