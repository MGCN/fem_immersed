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

#ifndef COMPUTEFLUX_H
#define COMPUTEFLUX_H

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
class ComputeFlux;

template<>
InputParameters validParams<ComputeFlux>();

/**
 * This postprocessor computes a side integral of the mass flux.
 */
class ComputeFlux : public SideIntegralVariablePostprocessor
{
public:
  ComputeFlux(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue &_u_old;
//  const VariableValue &_vel_y;
//  const VariableValue &_vel_z;

 };

#endif // SIDEFLUXINTEGRAL_H
