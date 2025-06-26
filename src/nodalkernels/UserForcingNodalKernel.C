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

#include "UserForcingNodalKernel.h"

registerMooseObject("ImmersedBoundaryApp", UserForcingNodalKernel);

template<>
InputParameters validParams<UserForcingNodalKernel>()
{
  InputParameters params = validParams<NodalKernel>();
  params.addRequiredCoupledVar("u_f", "The coupled variable which provides the velocity");
  return params;
}

UserForcingNodalKernel::UserForcingNodalKernel(const InputParameters & parameters) :
    NodalKernel(parameters),
    _u_f(coupledValue("u_f"))
{}

Real
UserForcingNodalKernel::computeQpResidual()
{
  return - _u_f[0];
}
