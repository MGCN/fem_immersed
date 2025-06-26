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

#ifndef IMPOSEVELOCITY_FUNCTION_H
#define IMPOSEVELOCITY_FUNCTION_H
#include "DirichletBCBase.h"


class ImposeVelocity_function : public DirichletBCBase
{
public:
  ImposeVelocity_function(const InputParameters & parameters);

protected:
  virtual Real computeQpValue();

  const VariableValue & _u_old;
  const Function * const _function;
 
};

template<>
InputParameters validParams<ImposeVelocity_function>();

#endif /*  IMPOSEVELOCITY_FUNCTION_H */
