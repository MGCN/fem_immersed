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

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#include "DensOnB.h"
#include "MooseMesh.h"

registerMooseObject("ImmersedBoundaryApp", DensOnB);

template<>
InputParameters validParams<DensOnB>()
{
    InputParameters params = validParams<NodalKernel>();
    params.addRequiredParam<Real>("rho_f", "density_fluid");
    params.addRequiredCoupledVar("velocity","velocity variable");
    params.addRequiredCoupledVar("acceleration","acceleration variable");
    params.addRequiredParam<Real>("beta","beta parameter for Newmark Time integration");
    params.addRequiredParam<Real>("gamma","gamma parameter for Newmark Time integration");
    params.addParam<Real>("eta",0,"eta parameter for mass dependent Rayleigh damping");
    params.addParam<Real>("alpha_m",0,"alpha parameter for CH method");
    return params;
}

DensOnB::DensOnB(const InputParameters & parameters) :
NodalKernel(parameters),
_rho_f(getParam<Real>("rho_f")),
_u_old(valueOld()),
_vel_old(coupledValueOld("velocity")),
_accel_old(coupledValueOld("acceleration")),
_beta(getParam<Real>("beta")),
_gamma(getParam<Real>("gamma")),
_eta(getParam<Real>("eta")),
_alpha_m(getParam<Real>("alpha_m"))


{
}

Real
DensOnB::computeQpResidual()
{
    
    if (_dt == 0)
        return 0;
    else
    {
        diff_rho = -_rho_f;
              
        Real accel = 1./_beta * (((_u[0] - _u_old[0])/( _dt * _dt )) - _vel_old[0]/_dt - _accel_old[0] * (0.5 -_beta));
        
        Real vel = _vel_old[0] + ( _dt * ( 1 - _gamma )) * _accel_old[0] + _gamma * _dt * accel;
        
        return (1.0 - _alpha_m) * diff_rho * accel + _alpha_m * _accel_old[0] * diff_rho  +   diff_rho * vel * _eta;
    }
    
}
