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

#include "VelocityL2Error.h"
#include "MooseVariable.h"

registerMooseObject("ImmersedBoundaryApp", VelocityL2Error);

template <>
InputParameters
validParams<VelocityL2Error>()
{
  InputParameters params = validParams<ElementIntegralVariablePostprocessor>();
  return params;
}

VelocityL2Error::VelocityL2Error(const InputParameters & parameters)
  : ElementIntegralVariablePostprocessor(parameters)
{
}

Real
VelocityL2Error::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real
VelocityL2Error::computeQpIntegral()
{

//  MooseVariable & _var1 = _fe_problem.getVariable(0,"velocity_sm");
//  MooseVariable & _var2 = _fe_problem.getVariable(0,"velocity_fm");
//  Real diff = _var1.getValue(_current_elem,0) - _var2.getValue(_current_elem,0);
  return 0.0;
}
