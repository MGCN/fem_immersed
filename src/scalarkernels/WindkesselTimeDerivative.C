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

#include "WindkesselTimeDerivative.h"

registerMooseObject("ImmersedBoundaryApp", WindkesselTimeDerivative);

template<>
InputParameters validParams<WindkesselTimeDerivative>()
{
  InputParameters params = validParams<ODEKernel>();
  return params;
}

WindkesselTimeDerivative::WindkesselTimeDerivative(const InputParameters & parameters) :
    ODEKernel(parameters)
{
}

Real
WindkesselTimeDerivative::computeQpResidual()
{
  return 0.0; //_u_dot[_i];
}

Real
WindkesselTimeDerivative::computeQpJacobian()
{
  if (_i == _j)
    return 0.0; //_du_dot_du[_i];
  else
    return 0;
}

