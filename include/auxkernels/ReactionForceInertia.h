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
/*          Auxkernel to compute acceleration field             */
/****************************************************************/

#ifndef ReactionForceInertia_H
#define ReactionForceInertia_H

#include "AuxKernel.h"

class ReactionForceInertia;

template<>
InputParameters validParams<ReactionForceInertia>();

class ReactionForceInertia : public AuxKernel
{
public:
    
    /**
     *Computes Acceleration using Newmark Time integration scheme
     */
   ReactionForceInertia(const InputParameters & parameters);
    
    virtual ~ReactionForceInertia() {}
    
protected:
    virtual Real computeValue();
    
    const VariableValue & _var1;
    Real _rho_s;
    Real _rho_f;
    

};

#endif //SumVariables_H
