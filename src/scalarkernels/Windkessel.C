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

#include "Windkessel.h"
#include "Function.h"

registerMooseObject("ImmersedBoundaryApp", Windkessel);

template<>
InputParameters validParams<Windkessel>()
{
  InputParameters params = validParams<ODEKernel>();
  params.addRequiredParam<Real>("time_constant", "time constant");
  params.addRequiredParam<FunctionName>("function", "The forcing function");
  return params;
}

Windkessel::Windkessel(const InputParameters & parameters) :
    ODEKernel(parameters),
    _tau(getParam<Real>("time_constant")),
    _func(getFunction("function"))
  
{
}

Windkessel::~Windkessel()
{
}

Real
Windkessel::computeQpResidual()
{
  return  _tau * _u[_i] - _func.value(_t, 0);
}

Real
Windkessel::computeQpJacobian()
{
  return 0.0;
}

Real
Windkessel::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
