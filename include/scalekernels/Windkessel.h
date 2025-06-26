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

#ifndef Windkessel_H
#define Windkessel_H


#include "ODEKernel.h"

//Forward Declarations
class Windkessel;

template<>
InputParameters validParams<Windkessel>();

/**
 *
 */
class Windkessel : public ODEKernel
{
public:
  Windkessel(const InputParameters & parameters);
  virtual ~Windkessel();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  Real _tau;
  const Function & _func;

};


#endif /* Windkessel_H */
