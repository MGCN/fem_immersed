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
/*          Auxkernel to compute acceleration field             */
/****************************************************************/

#ifndef SumPostProc_H
#define SumPostProc_H

#include "AuxKernel.h"

class SumPostProc;

template<>
InputParameters validParams<SumPostProc>();

class SumPostProc : public AuxKernel
{
public:

  /**
  *Computes Acceleration using Newmark Time integration scheme
  */
  SumPostProc(const InputParameters & parameters);

  virtual ~SumPostProc() {}

protected:
  virtual Real computeValue();
  const PostprocessorValue & _postprocessor_value1;
  const PostprocessorValue & _postprocessor_value2;
//
//  const VariableValue & _var1;
//  const VariableValue & _var2;

};

#endif //SumVariables_H
