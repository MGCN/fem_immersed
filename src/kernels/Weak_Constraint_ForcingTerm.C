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
/*      It implements the forcing function applied weakly.      */
/*      It implements the term _int(_u_f() * test) dS           */
/****************************************************************/


#include "Weak_Constraint_ForcingTerm.h"

registerMooseObject("ImmersedBoundaryApp", Weak_Constraint_ForcingTerm);

template<>
InputParameters validParams<Weak_Constraint_ForcingTerm>()
{
   InputParameters params = validParams<Kernel>();
   params.addRequiredCoupledVar("u_f", "The coupled variable which provides the velocity"); 
   return params;
}

Weak_Constraint_ForcingTerm::Weak_Constraint_ForcingTerm(const InputParameters & parameters) :
  Kernel(parameters),
  _u_f(coupledValue("u_f"))
  
{}  


Real
Weak_Constraint_ForcingTerm::computeQpResidual()
{    

    return  - 1.0  * _u_f[_qp] * _test[_i][_qp] * _dt;
    
}      
    
Real
Weak_Constraint_ForcingTerm::computeQpJacobian()
{ 
    return   0.0;
}

Real
Weak_Constraint_ForcingTerm::computeQpOffDiagJacobian(unsigned int j_var)
{      
    return 0.0;
}
