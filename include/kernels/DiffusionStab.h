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

#ifndef DIFFUSIONSTAB_H
#define DIFFUSIONSTAB_H

#include "Kernel.h"

class DiffusionStab;

template<>
InputParameters validParams<DiffusionStab>();


class DiffusionStab : public Kernel
{
public:
  DiffusionStab(const InputParameters & parameters);


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
  Real _Dstab;
};


#endif /* DIFFUSION_H */
