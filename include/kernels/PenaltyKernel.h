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
#ifndef PENALTYKERNEL_H
#define PENALTYKERNEL_H

#include "Kernel.h"

class PenaltyKernel;


template<>
InputParameters validParams<PenaltyKernel>();

/**
 * A different approach to applying Dirichlet BCs
 *
 * uses \f$ \int(p u \cdot \phi)=\int(p f \cdot \phi)\f$ on \f$d\omega\f$
 *
 */

class PenaltyKernel : public Kernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same constructor.
   */
PenaltyKernel(const InputParameters & parameters);

  

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  const VariableValue & _f;

  Real _p;
};

#endif
