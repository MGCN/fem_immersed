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
/*  Auxkernel to compute displacement field from velocity field */
/****************************************************************/

#ifndef FLUIDDISPLACEMENTAUX_H
#define FLUIDDISPLACEMENTAUX_H

#include "AuxKernel.h"

//Forward Declarations
class FluidDisplacementAux;

template<>
InputParameters validParams<FluidDisplacementAux>();
/**
 * Velocity auxiliary value
 */
class FluidDisplacementAux : public AuxKernel
{
public:
 
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
 FluidDisplacementAux(const InputParameters & parameters);
  
 
 virtual ~FluidDisplacementAux() {}

protected:
  virtual Real computeValue();

  
   
   const VariableValue& _u_vel;
   const VariableValue& _v_vel;
   const VariableValue& _w_vel;
   unsigned _component;

};

#endif //VELOCITYAUX_H
