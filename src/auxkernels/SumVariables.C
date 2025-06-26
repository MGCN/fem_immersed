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
/*              Prepared by Maria Nestola,                    */	
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */		
/*          Auxkernel to compute acceleration field             */
/****************************************************************/

#include "SumVariables.h"
registerMooseObject("ImmersedBoundaryApp", SumVariables);
template<>
InputParameters validParams<SumVariables>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("var1","variable1");
    params.addRequiredCoupledVar("var2","variable2");
    params.addCoupledVar("var3","variable3");
  return params;
}

SumVariables::SumVariables(const InputParameters & parameters) :
  AuxKernel(parameters),
   _var1(coupledValue("var1")),
   _var2(coupledValue("var2")),
   _var3(parameters.isParamValid("var3") ? coupledValue("var3"):_zero)


{
}

Real
SumVariables::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");

 
  return _var1[_qp] + _var2[_qp] + _var3[_qp];
}
