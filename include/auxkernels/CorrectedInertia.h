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

#ifndef CORRECTEDINERTIA_H
#define CORRECTEDINERTIA_H

#include "AuxKernel.h"

class CorrectedInertia;

template<>
InputParameters validParams<CorrectedInertia>();

class CorrectedInertia : public AuxKernel
{
public:
    
    /**
     *Computes Acceleration using Newmark Time integration scheme
     */
    CorrectedInertia(const InputParameters & parameters);
    
    virtual ~CorrectedInertia() {}
    
protected:
    virtual Real computeValue();
    
    const VariableValue & _inertial_residual;
    Real _fluid_density;
    Real _solid_density;
    std::vector<int> _boundary_tags;

};

#endif //CorrectedInertia_H
