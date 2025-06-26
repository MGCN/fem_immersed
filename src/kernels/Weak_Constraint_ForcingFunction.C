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
/*    Immersed_Boundary- ICS Mechanical simulation framework    */	
/*                Prepared by Maria Nestola,                    */	
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */	
/*      Kernel to impose an essential conditions weakly.        */
/*      It implements the forcing term  applied weakly.         */
/*      To be used with Weak_Constraint_multiplier              */
/*     It implements the term _int(_f() * test) * dS            */
/****************************************************************/

#include "Weak_Constraint_ForcingFunction.h"
#include "Function.h"

registerMooseObject("ImmersedBoundaryApp", Weak_Constraint_ForcingFunction);

template<>
InputParameters validParams<Weak_Constraint_ForcingFunction>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<FunctionName>("function", "The forcing function");
  return params;
}

Weak_Constraint_ForcingFunction::Weak_Constraint_ForcingFunction(const InputParameters & parameters) :
    Kernel(parameters),
    _func(getFunction("function"))
{
}

Real
Weak_Constraint_ForcingFunction::f()
{
  return  _func.value(_t, _q_point[_qp]);
}

Real
Weak_Constraint_ForcingFunction::computeQpResidual()
{
  return  -1.0 * _test[_i][_qp] * f();
}

