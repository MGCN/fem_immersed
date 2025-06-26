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

#ifndef NODALSUMR_H
#define NODALSUMR_H

#include "NodalVariablePostprocessor.h"

// Forward Declarations
class NodalSumR;

template<>
InputParameters validParams<NodalSumR>();

/**
 * Computes a sum of the nodal values of the coupled variable.
 */
class NodalSumR : public NodalVariablePostprocessor
{
public:
  NodalSumR(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();

  /**
   * This will return the degrees of freedom in the system.
   */
  virtual Real getValue();

  void threadJoin(const UserObject & y);
    

protected:
  MeshBase & _mesh;
  Real _volume;
  const Real & _current_side_volume;  
  Real _sum;

    

};

#endif //NODALSUM_H
