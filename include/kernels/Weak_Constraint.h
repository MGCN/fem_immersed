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


#ifndef WEAK_CONSTRAINT
#define WEAK_CONSTRAINT

#include "Kernel.h"


// Forward Declarations
class Weak_Constraint;

template<>
InputParameters validParams<Weak_Constraint>();

class Weak_Constraint: public Kernel
{
public:
 Weak_Constraint(const InputParameters & parameters);

protected:

Real f();
virtual Real computeQpResidual();
virtual Real computeQpJacobian();
virtual Real computeQpOffDiagJacobian(unsigned int);

private:
const VariableValue& _primal;
unsigned int _primal_var_number;


};
 
#endif 