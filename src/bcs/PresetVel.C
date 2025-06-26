//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PresetVel.h"
#include "Function.h"

registerMooseObject("ImmersedBoundaryApp", PresetVel);

template <>
InputParameters
validParams<PresetVel>()
{
  InputParameters params = validParams<DirichletBCBase>();
  params.addClassDescription(
      "Prescribe the displacement on a given boundary in a given direction.");
  params.addParam<Real>("scale_factor", 1, "Scale factor if function is given.");
  params.addRequiredCoupledVar("forcing_vel", "Velocity obtained from the fluid.");
  params.addRequiredCoupledVar("velocity", "The velocity variable.");
  params.addRequiredCoupledVar("acceleration", "The acceleration variable.");
  params.addRequiredParam<Real>("beta", "beta parameter for Newmark time integration.");
  params.addRequiredParam<Real>("gamma", "gamma parameter for Newmark time integration.");
  return params;
}

PresetVel::PresetVel(const InputParameters & parameters)
  : DirichletBCBase(parameters),
    _u_old(valueOld()),
    _scale_factor(parameters.get<Real>("scale_factor")),
    _force_vel(coupledValueOld("forcing_vel")),
    _force_vel_old(coupledValueOlder("forcing_vel")),
    _vel_old(coupledValueOld("velocity")),
    _accel_old(coupledValueOld("acceleration")),
    _beta(getParam<Real>("beta")),
    _gamma(getParam<Real>("gamma"))
{
}

Real
PresetVel::computeQpValue()
{
  Real accel = (_force_vel[_qp] - _vel_old[_qp]) / (_gamma*_dt) + (_gamma - 1.)/_gamma*_accel_old[_qp];

  return _u_old[_qp] + _dt * _vel_old[_qp] +
         ((0.5 - _beta) * _accel_old[_qp] + _beta * accel) * _dt * _dt;
}
