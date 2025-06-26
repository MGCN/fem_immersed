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
/*      It implements he forcing term applied weakly.           */
/*      To be used with Weak_Constraint_multiplier              */
/*      It implements the term _int(_f() * test) * dS           */
/****************************************************************/



#ifndef  WEAK_CONSTRAINT_FORCINGFUNCTION_H
#define  WEAK_CONSTRAINT_FORCINGFUNCTION_H
#include "Kernel.h"

//Forward Declarations
class Weak_Constraint_ForcingFunction;
class Function;
template<>
InputParameters validParams<Weak_Constraint_ForcingFunction>();

    class Weak_Constraint_ForcingFunction: public Kernel
    {
      public:
   
        Weak_Constraint_ForcingFunction(const InputParameters & parameters);
   
      protected:
        Real f();
        virtual Real computeQpResidual();
        const Function & _func;  
    };

#endif // WEAK_CONSTRAINT_FORCINGFUNCTION_H
