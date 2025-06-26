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

#ifndef USERFORCINGNODALKERNEL_H
#define USERFORCINGNODALKERNEL_H

#include "NodalKernel.h"

//Forward Declarations
class UserForcingNodalKernel;

template<>
InputParameters validParams<UserForcingNodalKernel>();

/**
 * Represents the rate in a simple ODE of du/dt = f
 */
class UserForcingNodalKernel : public NodalKernel
{
public:
  /**
   * Constructor grabs the Function
   */
  UserForcingNodalKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  const VariableValue& _u_f;
};

#endif
