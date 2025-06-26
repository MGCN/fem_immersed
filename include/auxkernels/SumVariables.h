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

#ifndef SumVariables_H
#define SumVariables_H

#include "AuxKernel.h"

class SumVariables;

template<>
InputParameters validParams<SumVariables>();

class SumVariables : public AuxKernel
{
public:
    
    /**
     *Computes Acceleration using Newmark Time integration scheme
     */
    SumVariables(const InputParameters & parameters);
    
    virtual ~SumVariables() {}
    
protected:
    virtual Real computeValue();
    
    const VariableValue & _var1;
    const VariableValue & _var2;
    const VariableValue & _var3;

};

#endif //SumVariables_H
