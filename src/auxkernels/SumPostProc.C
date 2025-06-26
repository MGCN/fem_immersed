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
/*              Prepared by Maria Nestola,                      */
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */		
/*          Auxkernel to compute acceleration field             */
/****************************************************************/

#include "SumPostProc.h"
registerMooseObject("ImmersedBoundaryApp",SumPostProc);
template<>
InputParameters validParams<SumPostProc>()
{
  InputParameters params = validParams<AuxKernel>();
//    params.addRequiredCoupledVar("var1","variable1");
//    params.addRequiredCoupledVar("var2","variable2");
    params.addRequiredParam<PostprocessorName>("postprocessor1", "The postprocessor to set the value to on the boundary.");
    params.addRequiredParam<PostprocessorName>("postprocessor2", "The postprocessor to set the value to on the boundary.");
  return params;
}

SumPostProc::SumPostProc(const InputParameters & parameters) :
  AuxKernel(parameters),
  _postprocessor_value1(getPostprocessorValue("postprocessor1")),
  _postprocessor_value2(getPostprocessorValue("postprocessor2"))
//   _var1(coupledValue("var1")),
//   _var2(coupledValue("var2"))

{
}

Real
SumPostProc::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");

 
  return -1.0 * (_postprocessor_value1 + _postprocessor_value2);
}
