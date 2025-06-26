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

#ifndef WEIGHTEDVARIABLE_H
#define WEIGHTEDVARIABLE_H

#include "AuxKernel.h"

class WeightedVariable;

template<>
InputParameters validParams<WeightedVariable>();

class WeightedVariable : public AuxKernel
{
public:

  /**
  *Computes Acceleration using Newmark Time integration scheme
  */
  WeightedVariable(const InputParameters & parameters);

  virtual ~WeightedVariable() {}

protected:
  virtual Real computeValue();

  const VariableValue & _var_old;
  const VariableValue & _var;
  Real _weight;

};

#endif //WEIGHTEDVARIABLE_H
