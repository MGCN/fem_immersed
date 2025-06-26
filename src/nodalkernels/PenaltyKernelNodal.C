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
#include "PenaltyKernelNodal.h"
#include <cmath>
#include <iostream>     
#include <algorithm>

registerMooseObject("ImmersedBoundaryApp", PenaltyKernelNodal);

template<>

InputParameters validParams<PenaltyKernelNodal>()
{
    InputParameters params = validParams<NodalKernel>();
    params.addRequiredParam<Real>("penalty", "Penalty scalar");
    params.addRequiredCoupledVar("control", "force variable");
    return params;
}

PenaltyKernelNodal::PenaltyKernelNodal(const InputParameters & parameters) :
NodalKernel(parameters),
_f(coupledValue("control")),
_p(getParam<Real>("penalty"))
{}

Real
PenaltyKernelNodal::computeQpResidual()
{
    if ((std::abs(_f[0]))>=0.7)
    {
      return _p * _u[0];
    }
    else{
    return 0.0;
    }
}


Real
PenaltyKernelNodal::computeQpJacobian()
{
  return  0.0;
}
