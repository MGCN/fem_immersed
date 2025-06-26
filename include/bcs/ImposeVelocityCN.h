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
/*                Prepared by Barna Becsek,                     */	
/*                                                              */	
/*      NodalBC to set the boundary velocity.                   */ 
/****************************************************************/

#ifndef IMPOSEVELOCITYCN_H
#define IMPOSEVELOCITYCN_H

#include "DirichletBCBase.h"


class ImposeVelocityCN : public DirichletBCBase
{
public:
   ImposeVelocityCN(const InputParameters & parameters);

protected:
  virtual Real computeQpValue();

  const VariableValue & _u_old;

  const VariableValue & _u_older;

  const VariableValue& _vel_s;
    
  const VariableValue& _vel_s_old;

  const VariableValue& _vel_s_older;

//  const Function * const _function;
  
  Real _theta;
};

template<>
InputParameters validParams<ImposeVelocityCN>();

#endif /* IMPOSEVELOCITYCN_H */
