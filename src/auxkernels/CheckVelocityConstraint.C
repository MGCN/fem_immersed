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
/*                Prepared by Maria Nestola,                    */	
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */		
/*  Auxkernel to compute displacement field from velocity field */
/****************************************************************/

#include "CheckVelocityConstraint.h"
#include "MooseMesh.h"
registerMooseObject("ImmersedBoundaryApp", CheckVelocityConstraint);
template<>
InputParameters validParams<CheckVelocityConstraint>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("u", "velocity");
  params.addRequiredCoupledVar("disp", "disp");
  return params;
}

CheckVelocityConstraint::CheckVelocityConstraint(const InputParameters & parameters) :
    AuxKernel(parameters),
    _u_vel(coupledValue("u")),
    _dispOld(coupledValueOld("disp")),
    _disp(coupledValue("disp"))
{}

Real CheckVelocityConstraint::computeValue()
{    
   
    return _disp[_qp] - (_dispOld[_qp] + _u_vel[_qp] * _dt);
   
}

