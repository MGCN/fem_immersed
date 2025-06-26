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
#include "PenaltyNodalAuxKernel.h"
#include "MooseMesh.h"
registerMooseObject("ImmersedBoundaryApp", PenaltyNodalAuxKernel);
template<>
InputParameters validParams<PenaltyNodalAuxKernel>()
{
    InputParameters params = validParams<AuxKernel>();
    params.addRequiredParam<Real>("penalty", "Penalty scalar");
    params.addRequiredParam<unsigned int>("component", "The desired value.");
    params.addRequiredCoupledVar("vel_x", "x-velocity");
    params.addCoupledVar("vel_y", "y-velocity"); //only required in 2D and in 3D;
    params.addCoupledVar("vel_z", "z-velocity"); //only required in 3D;
    params.addRequiredCoupledVar("disp_x", "x-disp");
    params.addCoupledVar("disp_y", "y-disp"); //only required in 2D and in 3D;
    params.addCoupledVar("disp_z", "z-disp"); //only required in 3D;
    
    return params;
}

PenaltyNodalAuxKernel::PenaltyNodalAuxKernel(const InputParameters & parameters) :
 AuxKernel(parameters),
 _component(getParam<unsigned int>("component")),
 _vel_x(coupledValue("vel_x")),
 _vel_y(_mesh.dimension()>=2 ? coupledValue("vel_y"): _zero),
 _vel_z(_mesh.dimension()==3 ? coupledValue("vel_z"): _zero),
 _disp_x_old(coupledValueOld("disp_x")),
 _disp_y_old(_mesh.dimension()>=2 ? coupledValueOld("disp_y"): _zero),
 _disp_z_old(_mesh.dimension()==3 ? coupledValueOld("disp_z"): _zero),
 _p(getParam<Real>("penalty"))

{}

Real
PenaltyNodalAuxKernel::computeValue()
{
    if (_component==0){
        std::cout<< _disp_x_old[_qp] << std::endl;
        return _p * (_disp_x_old[_qp] + _vel_x[_qp] * _dt);}
    if (_component==1){
        return _p * (_disp_y_old[_qp] + _vel_y[_qp] * _dt);}
    if (_component==2){
        return _p * (_disp_z_old[_qp] + _vel_z[_qp] * _dt) ;}
    else
    {
        std::cout<<"qua non dovresti esserci NODLA POSTPROCESSORAUX\n\n";
        return 0.0;
    }
    
}



