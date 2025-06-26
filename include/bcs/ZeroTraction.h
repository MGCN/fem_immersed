/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ZeroTraction_H
#define ZeroTraction_H

#include "IntegratedBC.h"

// Forward Declarations
class ZeroTraction;

template<>
InputParameters validParams<ZeroTraction>();

/**
 * This class implements the "No BC" boundary condition
 * discussed by Griffiths, Papanastiou, and others.
 */
class ZeroTraction : public IntegratedBC
{
public:
  ZeroTraction(const InputParameters & parameters);

  virtual ~ZeroTraction(){}

protected:
  virtual Real computeQpResidual();


  unsigned _component;
  const Function & _func;

};


#endif // ZeroTraction_H
