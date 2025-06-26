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

#ifndef USERFORCING_H
#define USERFORCING_H

#include "Kernel.h"

//Forward Declarations
class UserForcing;
class Function;

template<>
InputParameters validParams<UserForcing>();

/**
 * Define the Kernel for a user defined forcing function that looks like:
 *
 * test function * forcing function
 */
class UserForcing : public Kernel
{
public:

  UserForcing(const InputParameters & parameters);

protected:
  /**
   * Evaluate f at the current quadrature point.


  /**
   * Computes test function * forcing function.
   */
protected:
    virtual Real computeQpResidual();
    virtual Real computeQpJacobian();
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);
    
    const int _component;
    
    const VariableValue & _pressureFluid;
    const VariableGradient &_grad_disp_x;
    const VariableGradient &_grad_disp_y;
    const VariableGradient &_grad_disp_z;
    const VariableGradient &_grad_vel_x;
    const VariableGradient &_grad_vel_y;
    const VariableGradient &_grad_vel_z;
    const Function &_function2;
    const MooseArray<Point> & _normals;
    
    Real _mu;
    RealTensorValue sigma;

    
};

#endif //USERFORCINGFUNCTION_H
