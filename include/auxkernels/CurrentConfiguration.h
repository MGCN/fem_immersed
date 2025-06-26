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
/*         Auxkernel to compute displacement field              */
/*         in the current configuration                         */
/****************************************************************/

#ifndef CurrentConfiguration_H
#define CurrentConfiguration_H

#include "AuxKernel.h"

class CurrentConfiguration;

template<>
InputParameters validParams<CurrentConfiguration>();

class CurrentConfiguration : public AuxKernel
{
public:
    CurrentConfiguration(const InputParameters &parameters);
    
    ~CurrentConfiguration() {};
protected:

   virtual Real computeValue();
   unsigned _component;

   
  const VariableValue & _disp_x;
  const VariableValue & _disp_y;
  const VariableValue & _disp_z;


};

#endif // ELASTICENERGYAUX_H
