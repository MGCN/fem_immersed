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
/*                Prepared by Barna Becsek,                     */	
/*                  ARTORG, Uni Bern, Bern                      */	
/*                                                              */		
/*    Auxkernel to weigh a variable with a previous value       */
/****************************************************************/

#include "WeightedVariable.h"
registerMooseObject("ImmersedBoundaryApp", WeightedVariable);
template<>
InputParameters validParams<WeightedVariable>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("var","variable to read");
    params.addRequiredParam<Real>("weight","weight to use");
  return params;
}

WeightedVariable::WeightedVariable(const InputParameters & parameters) :
  AuxKernel(parameters),
   _var_old(coupledValueOld("var")),
   _var(coupledValue("var")),
   _weight(getParam<Real>("weight"))
{
}

Real
WeightedVariable::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");

  //Calculates acceeleration using Newmark time integration method
  return _weight*_var[_qp] + (1 - _weight)*_var_old[_qp];
}
