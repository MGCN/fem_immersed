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

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       HART - ICS Electromechanical simulation framework      */
/*                Prepared by Mathias Winkel,                   */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* AuxKernel to compute the cardiac thickness i.e. point's      */
/* location inside the left ventricular mesh.                   */
/****************************************************************/

#include "CardiacThicknessParameterFiber.h"
registerMooseObject("ImmersedBoundaryApp", CardiacThicknessParameterFiber);
template <>
InputParameters
validParams<CardiacThicknessParameterFiber>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar(
      "distance_LV_inner",
      "Aux variable of the distance to the inner wall of the left ventricle.");
  params.addRequiredCoupledVar(
      "distance_outer",
      "Aux variable of the distance to the outer wall of the heart.");

  return params;
}

CardiacThicknessParameterFiber::CardiacThicknessParameterFiber(
    const InputParameters &parameters)
    : AuxKernel(parameters), _d_lv(coupledValue("distance_LV_inner")),
      _d_o(coupledValue("distance_outer"))
{
}

Real
CardiacThicknessParameterFiber::computeValue()
{
  // @todo: maybe a cleaner regularization instead of adding 1.e-16 to the
  // denominator might make sense...
  // left ventricle --> positive thickness parameter
  return _d_o[_qp] /
         (_d_o[_qp] + _d_lv[_qp]); // here there was 1.-16 at denominator
}
