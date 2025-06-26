/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "Outflow.h"
registerMooseObject("ImmersedBoundaryApp", Outflow);
template<>
InputParameters validParams<Outflow>()
{
  InputParameters params = validParams<NodalBC>();
  params .addRequiredParam<PostprocessorName>("postprocessor", "The postprocessor to set the value to on the boundary.");
  params.addRequiredParam<Real>("Resistence", "Resistence");
  return params;
}



Outflow::Outflow(const InputParameters & parameters) :
NodalBC(parameters),

  // Coupled variables
  _postprocessor_value(getPostprocessorValue("postprocessor")),
  _resistence(getParam<Real>("Resistence"))

{
}



Real Outflow::computeQpResidual()
{
  
  return _u[_qp] - _resistence *_postprocessor_value ;
}




