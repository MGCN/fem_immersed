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
/*    MECH - ICS Mechanical simulation framework                */
/*                Prepared by Barna Becsek,                     */
/*                                                              */
/*      NodalBC to set the boundary velocity.                   */
/****************************************************************/


#include "ImposeVelocityCN.h"
#include "Function.h"
registerMooseObject("ImmersedBoundaryApp", ImposeVelocityCN);
template<>
InputParameters validParams<ImposeVelocityCN>()
{
    InputParameters params = validParams<DirichletBCBase>();
    params.addCoupledVar("vel_s", "velocity_xs");
    params.addParam<FunctionName>("function", "", "Function describing the velocity.");
    params.addParam<Real>("theta", 1., "theta parameter for CN integration");
    return params;
}


ImposeVelocityCN::ImposeVelocityCN(const InputParameters & parameters) :
DirichletBCBase(parameters),
_u_old(valueOld()),
_u_older(valueOlder()),
_vel_s(coupledValue("vel_s")),
_vel_s_old(coupledValueOld("vel_s")),
_vel_s_older(coupledValueOlder("vel_s")),
//_function(parameters.get<FunctionName>("function") != "" ? &getFunction("function") : NULL),
_theta(getParam<Real>("theta"))
{
}

Real
ImposeVelocityCN::computeQpValue()
{
    if (_t_step==1)
        return  _u_old[_qp] + _dt * _vel_s[_qp];
    else {
        //std::cerr << "timestep: " << _t_step << " v: " << _vel_s[_qp] << " v_old: " << _vel_s_old[_qp] << " v_older[_qp]: " << _vel_s_older[_qp] << std::endl;
        //return _u_older[_qp] +  _dt/3.*(_vel_s[_qp] + 4.* _vel_s_old[_qp] + _vel_s_older[_qp]); // Simpson's rule
        return _u_old[_qp] + _dt * _vel_s[_qp] * _theta + _dt_old * _vel_s_old[_qp] * (1. - _theta);
        }
    
}
