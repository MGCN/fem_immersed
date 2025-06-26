/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef Outflow_H
#define Outflow_H


#include "NodalBC.h"

// Forward Declarations
class Outflow;
class GeneralPostprocessor;

template<>
InputParameters validParams<Outflow>();

/**
 * This class implements the "No BC" boundary condition
 * discussed by Griffiths, Papanastiou, and others.
 */
class Outflow : public NodalBC
{
public:
  Outflow(const InputParameters & parameters);

  virtual ~Outflow(){}

protected:
  virtual Real computeQpResidual();
  const PostprocessorValue & _postprocessor_value;
  Real _resistence;



};


#endif // Outflow_H
