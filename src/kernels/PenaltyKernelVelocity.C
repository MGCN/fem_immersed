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
#include "PenaltyKernelVelocity.h"
#include <cmath>
#include <iostream>     
#include <algorithm>

registerMooseObject("ImmersedBoundaryApp", PenaltyKernelVelocity);

template<>
InputParameters validParams<PenaltyKernelVelocity>()
{
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<Real>("penalty", "Penalty scalar");
    params.addRequiredCoupledVar("velocity", "velocity variable");
    return params;
}

PenaltyKernelVelocity::PenaltyKernelVelocity(const InputParameters & parameters) :
Kernel(parameters),
_v(coupledValue("velocity")),
_u_old(valueOld()),
_p(getParam<Real>("penalty"))
{}

Real
PenaltyKernelVelocity::computeQpResidual()
{
  return _p * _test[_i][_qp] * ( _u[_qp] - _u_old[_qp]  - _v[_qp] * _dt);
}


Real
PenaltyKernelVelocity::computeQpJacobian()
{
  return  _p * _phi[_j][_qp] * _test[_i][_qp];
}

