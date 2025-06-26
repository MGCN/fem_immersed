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
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*      NodalBC to set the boundary velocity.                   */
/****************************************************************/


#include "ImposeVelocity.h"
#include "Function.h"
registerMooseObject("ImmersedBoundaryApp", ImposeVelocity);
template<>
InputParameters validParams<ImposeVelocity>()
{
    InputParameters params = validParams<DirichletBCBase>();
    params.addCoupledVar("vel_s", "velocity_xs");
    params.addParam<FunctionName>("function", "", "Function describing the velocity.");
    return params;
}


ImposeVelocity::ImposeVelocity(const InputParameters & parameters) :
DirichletBCBase(parameters),
_u_old(valueOld()),
_u_older(valueOlder()),
_vel_s(coupledValue("vel_s")),
_vel_s_old(coupledValueOld("vel_s")),
_function(parameters.get<FunctionName>("function") != "" ? &getFunction("function") : NULL)
{
}

Real
ImposeVelocity::computeQpValue()
{
    if (_t_step==1)
        return  _u_old[_qp] + _dt * _vel_s[_qp];
    else
        return _u_old[_qp] + _dt * _vel_s[_qp];
//_u_old[_qp] +  3./2 * _dt * _vel_s[_qp] - 1./2 * _dt * _vel_s_old[_qp] ;
    
}
