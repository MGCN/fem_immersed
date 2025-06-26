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

#ifndef CardiacThicknessParameterFiber_H
#define CardiacThicknessParameterFiber_H

#include "AuxKernel.h"

class CardiacThicknessParameterFiber;

template <>
InputParameters validParams<CardiacThicknessParameterFiber>();

/**
 * AuxKernel for computing a thickness parameter for a cardiac model.
 *
 * See the decription of CardiacThicknessParameterAux.
 */
class CardiacThicknessParameterFiber : public AuxKernel
{
public:
    CardiacThicknessParameterFiber(const InputParameters &parameters);

  // virtual ~CardiacThicknessParameterAuxLV() {}

protected:
  virtual Real computeValue();

  const VariableValue &_d_lv, &_d_o;
};
#endif
