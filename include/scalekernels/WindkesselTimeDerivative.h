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

#ifndef WindkesselTimeDerivative_H
#define WindkesselTimeDerivative_H

#include "ODEKernel.h"

// Forward Declaration
class WindkesselTimeDerivative;
class Function;

template<>
InputParameters validParams<WindkesselTimeDerivative>();

class WindkesselTimeDerivative : public ODEKernel
{
public:
  WindkesselTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

};

#endif //WindkesselTimeDerivative_H
