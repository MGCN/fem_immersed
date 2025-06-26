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
#include "PenaltyKernel.h"
#include <cmath>
#include <iostream>     
#include <algorithm>
registerMooseObject("ImmersedBoundaryApp", PenaltyKernel);
template<>
InputParameters validParams<PenaltyKernel>()
{
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<Real>("penalty", "Penalty scalar");
    params.addRequiredCoupledVar("control", "force variable");
    return params;
}

PenaltyKernel::PenaltyKernel(const InputParameters & parameters) :
Kernel(parameters),
_f(coupledValue("control")),
_p(getParam<Real>("penalty"))
{}

Real
PenaltyKernel::computeQpResidual()
{
    if ((std::abs(_f[_qp]))>=0.7)
    {
      return _p * _test[_i][_qp] * _u[_qp];
    }
    else{
    return 0.0;
    }
}

Real
PenaltyKernel::computeQpJacobian()
{
    
    if ((std::abs(_f[_qp]))>=0.7)
    { 
       return  _p * _phi[_j][_qp] * _test[_i][_qp];
    }
    else{
     return 0.0;
    }
    
}

