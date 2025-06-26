/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ZeroTraction.h"
#include "MooseMesh.h"
#include "Function.h"
registerMooseObject("ImmersedBoundaryApp", ZeroTraction);
template<>
InputParameters validParams<ZeroTraction>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<FunctionName>("function", "The function.");
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");

  return params;
}



ZeroTraction::ZeroTraction(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _component(getParam<unsigned>("component")),
  _func(getFunction("function"))

{
}



Real ZeroTraction::computeQpResidual()
{
    return -_test[_i][_qp] * _func.value(_t, _q_point[_qp])*_normals[_qp](_component);
}




