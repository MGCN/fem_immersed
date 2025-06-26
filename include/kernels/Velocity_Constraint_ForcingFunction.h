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

#ifndef VELOCITY_CONSTRAINT_FORCINGFUNCTION
#define VELOCITY_CONSTRAINT_FORCINGFUNCTION

#include "Kernel.h"


// Forward Declarations
class Velocity_Constraint_ForcingFunction;

template<>
InputParameters validParams<Velocity_Constraint_ForcingFunction>();

class Velocity_Constraint_ForcingFunction : public Kernel
{
public:
  Velocity_Constraint_ForcingFunction(const InputParameters & parameters);

protected:
  Real f();
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
  
  
private:
  const Function & _func;
  const  VariableValue& _slave;
  const VariableValue& _slave_old;
  unsigned int _slave_var_number;
};
 
#endif 
