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
/*    MECH - ICS Mechanical simulation framework                */		
/*                Prepared by Maria Nestola,                    */	
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */	
/*     NodalBC to set a forcing function as  boundary velocity. */ 
/****************************************************************/

#include "ImposeVelocity_function.h"
#include "Function.h"
registerMooseObject("ImmersedBoundaryApp", ImposeVelocity_function);
template<>
InputParameters validParams<ImposeVelocity_function>()
{
  InputParameters params = validParams<DirichletBCBase>();
  params.addParam<FunctionName>("function", "", "Function describing the velocity.");
  return params;
}


ImposeVelocity_function::ImposeVelocity_function(const InputParameters & parameters) :
  DirichletBCBase(parameters),
  _u_old(valueOld()),
  _function(parameters.get<FunctionName>("function") != "" ? &getFunction("function") : NULL)  

{}

Real
ImposeVelocity_function::computeQpValue()
{
  
 //Real vel( function->value(_t, *_current_node) );
    
    return  0.0;
}
