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

#ifndef ELEMENTL2ERROR_H
#define ELEMENTL2ERROR_H

#include "ElementIntegralVariablePostprocessor.h"

class Function;

// Forward Declarations
class VelocityL2Error;

template <>
InputParameters validParams<VelocityL2Error>();

class VelocityL2Error : public ElementIntegralVariablePostprocessor
{
public:
  VelocityL2Error(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  virtual Real computeQpIntegral() override;

};

#endif // ELEMENTL2ERROR_H
