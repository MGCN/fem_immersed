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

#ifndef WEAK_CONSTRAINT_FORCINGTERM
#define WEAK_CONSTRAINT_FORCINGTERM


#include "Kernel.h"

class Weak_Constraint_ForcingTerm;

template<>
InputParameters validParams<Weak_Constraint_ForcingTerm>();

class Weak_Constraint_ForcingTerm: public Kernel
{
public:
 Weak_Constraint_ForcingTerm(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);

private:
 
  
 const VariableValue& _u_f;
  
};
 
#endif 