/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NeumannStresses_H
#define NeumannStresses_H

#include "IntegratedBC.h"


class NeumannStresses;

template<>
InputParameters validParams<NeumannStresses>();

/**
 * Pressure applies a pressure on a given boundary in the direction defined by component
 */
class NeumannStresses : public IntegratedBC
{
public:
    NeumannStresses(const InputParameters & parameters);
    
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
    //    RealTensorValue D;
    //    RealTensorValue p_I;
    //    RealVectorValue zeros;
};

#endif //PRESSURE_H
