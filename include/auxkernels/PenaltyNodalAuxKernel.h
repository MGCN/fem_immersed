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

#ifndef PenaltyNodalAuxKernel_H
#define PenaltyNodalAuxKernel_H

#include "AuxKernel.h"

class PenaltyNodalAuxKernel;

template<>
InputParameters validParams<PenaltyNodalAuxKernel>();

class PenaltyNodalAuxKernel : public AuxKernel
{
public:
   PenaltyNodalAuxKernel(const InputParameters &parameters);
    
    ~PenaltyNodalAuxKernel() {};
protected:

   virtual Real computeValue();
   unsigned _component;
   const VariableValue & _vel_x;
   const VariableValue & _vel_y;
   const VariableValue & _vel_z;
   const VariableValue & _disp_x_old;
   const VariableValue & _disp_y_old;
   const VariableValue & _disp_z_old;
   Real _p;


};

#endif // ELASTICENERGYAUX_H
