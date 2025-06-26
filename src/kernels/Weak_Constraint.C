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
/*   Kernel to impose an essential conditions weakly.           */
/*   To be used with:                                           */
/*     1)  Weak_Constraint_Forcing_Function for                 */
/*         the non homogeneous case.                            */
/*     2)  Weak_Constraint_multiplier                           */
/*     It implements the term _int(_u * test)                   */
/****************************************************************/


#include "Weak_Constraint.h"
#include "Function.h"

registerMooseObject("ImmersedBoundaryApp", Weak_Constraint);

template<>
InputParameters validParams<Weak_Constraint>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("primal", "primal-variable"); 
  return params;
}

Weak_Constraint::Weak_Constraint(const InputParameters & parameters) :
  Kernel(parameters),
  _primal(coupledValue("primal")),
  _primal_var_number(coupled("primal"))

{}  

Real



Weak_Constraint::computeQpResidual()
{    
  return  _primal[_qp]*_test[_i][_qp]/* + 0.001*_u[_qp]*_test[_i][_qp]*/;  
}      
    
Real
Weak_Constraint::computeQpJacobian()
{ 
  return 0.0; //0.001*_phi[_j][_qp]*_test[_i][_qp];
}

Real
Weak_Constraint::computeQpOffDiagJacobian(unsigned int j_var)
{     
  if (j_var == _primal_var_number)
    
      return _phi[_j][_qp] * _test[_i][_qp];
          
  else 
    
      return 0.0;

}
