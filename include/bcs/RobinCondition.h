/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RobinCondition_H
#define RobinCondition_H

#include "IntegratedBC.h"
#include "ElasticMaterial.h"
#include "LinearElastic.h"

class RobinCondition;

template<>
InputParameters validParams<RobinCondition>();

/**
 * Pressure applies a pressure on a given boundary in the direction defined by component
 */
class RobinCondition : public IntegratedBC
{
public:
    RobinCondition(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual();
    virtual Real computeQpJacobian();
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);
    
    const int _component;
    Real _alpha_disp;
    Real _alpha_stress;
    MaterialProperty<RealTensorValue> const & _P;
    const MaterialProperty<ElasticMaterial *> &_materialPointer;
    //    RealTensorValue D;
    //    RealTensorValue p_I;
    //    RealVectorValue zeros;
};

#endif //PRESSURE_H
