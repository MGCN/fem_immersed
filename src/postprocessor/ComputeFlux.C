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
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#include "ComputeFlux.h"

registerMooseObject("ImmersedBoundaryApp", ComputeFlux);

template<>
InputParameters validParams<ComputeFlux>()
{
  InputParameters params = validParams<SideIntegralVariablePostprocessor>();


  return params;
}

ComputeFlux::ComputeFlux(const InputParameters & parameters) :
    SideIntegralVariablePostprocessor(parameters),
    _u_old(valueOld())

{}

Real
ComputeFlux::computeQpIntegral()
{
    return _u_old[_qp];
}

