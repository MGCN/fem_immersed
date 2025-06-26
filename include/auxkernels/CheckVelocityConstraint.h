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

#ifndef CheckVelocityConstraint_h
#define CheckVelocityConstraint_h


#include "AuxKernel.h"

class  CheckVelocityConstraint;

template<>
InputParameters validParams<CheckVelocityConstraint>();

class CheckVelocityConstraint: public AuxKernel
{
public:
 CheckVelocityConstraint(const InputParameters & parameters);
  
  virtual ~CheckVelocityConstraint() {}

protected:
  virtual Real computeValue();

private:
 
  
 const VariableValue& _u_vel;
 const VariableValue& _dispOld;
 const VariableValue& _disp;
  
};
 
#endif 