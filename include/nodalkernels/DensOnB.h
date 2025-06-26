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
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/


#ifndef DensOnB_H
#define DensOnB_H

#include "NodalKernel.h"
#include <petscsnes.h>
//#include <passo.h>
//#include <utopia.hpp>
//#include <petsc_libmesh_utilities.h>

//Forward Declarations
class DensOnB;

template<>
InputParameters validParams<DensOnB>();

/**
 * Represents the rate in a simple ODE of du/dt = f
 */
class DensOnB : public NodalKernel
{
public:
    /**
     * Constructor grabs the Function
     */
    DensOnB(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;
    Real _rho_f;
    Real diff_rho; 
private:
    
    const VariableValue & _u_old;
    const VariableValue & _vel_old;
    const VariableValue & _accel_old;
    const Real _beta=0.0;
    const Real _gamma=0.0;
    const Real _eta=0.0;
    const Real _alpha_m=0.0;
};

#endif
