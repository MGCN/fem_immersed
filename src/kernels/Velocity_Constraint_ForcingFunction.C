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
/*      To be used with Weak_Constraint_multiplier              */
/****************************************************************/

#include "Velocity_Constraint_ForcingFunction.h"
#include "Function.h"

registerMooseObject("ImmersedBoundaryApp", Velocity_Constraint_ForcingFunction);

template<>
InputParameters validParams<Velocity_Constraint_ForcingFunction>()
{
   InputParameters params = validParams<Kernel>();
   params.addRequiredParam<FunctionName>("function", "The forcing function.");
   params.addRequiredCoupledVar("slave", "slave-variable"); 
   return params;
}

Velocity_Constraint_ForcingFunction::Velocity_Constraint_ForcingFunction(const InputParameters & parameters) :
  Kernel(parameters),
  _func(getFunction("function")),
  _slave(coupledValue("slave")),
  _slave_old(coupledValueOld("slave")),
  _slave_var_number(coupled("slave"))  
{
}  

Real
Velocity_Constraint_ForcingFunction::f()
{
  return 0.0;
}

Real
Velocity_Constraint_ForcingFunction::computeQpResidual()
{    
  Real g = 0;
  return ((_slave[_qp] - _slave_old[_qp])/_dt - _func.value(_t, _q_point[_qp]))*_test[_i][_qp];    
}      
    
Real
Velocity_Constraint_ForcingFunction::computeQpJacobian()
{ 
   return  0;
}

Real
Velocity_Constraint_ForcingFunction::computeQpOffDiagJacobian(unsigned int j_var)
{      
   if (j_var == _slave_var_number)
      return _phi[_j][_qp]*_test[_i][_qp]/_dt;
   else 
      return 0.0;
}
