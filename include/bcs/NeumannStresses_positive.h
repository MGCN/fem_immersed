/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NeumannStresses_positive_H
#define NeumannStresses_positive_H

#include "IntegratedBC.h"


class NeumannStresses_positive;

template<>
InputParameters validParams<NeumannStresses_positive>();

/**
 * Pressure applies a pressure on a given boundary in the direction defined by component
 */
class NeumannStresses_positive : public IntegratedBC
{
public:
    NeumannStresses_positive(const InputParameters & parameters);
    
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
    Real _mu;
    RealTensorValue sigma;

};

#endif //PRESSURE_H
