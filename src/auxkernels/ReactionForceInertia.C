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
/*    Immersed_Boundary- ICS Mechanical simulation framework    */
/*              Prepared by Maria Nestola,                      */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*          Auxkernel to compute acceleration field             */
/****************************************************************/

#include "ReactionForceInertia.h"
registerMooseObject("ImmersedBoundaryApp", ReactionForceInertia);
template<>
InputParameters validParams<ReactionForceInertia>()
{
    InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("var1","variable1");
    params.addRequiredParam<Real>("rho_f", "density_fluid");
    params.addRequiredParam<Real>("rho_s", "density_solid");
    
    return params;
}

ReactionForceInertia::ReactionForceInertia(const InputParameters & parameters) :
AuxKernel(parameters),
_var1(coupledValue("var1")),
_rho_s(getParam<Real>("rho_s")),
_rho_f(getParam<Real>("rho_f"))

{
}

Real
ReactionForceInertia::computeValue()
{
    
    return ( _rho_s - _rho_f ) * _var1[_qp];
}
