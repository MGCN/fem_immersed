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
/*     Immersed_Boundary-ICS Mechanical simulation framework    */
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*      Kernel to add e diffusion  stabilization term           */
/*                                                              */
/****************************************************************/
#include "DiffusionStab.h"
registerMooseObject("ImmersedBoundaryApp", DiffusionStab);
template<>
InputParameters validParams<DiffusionStab>()
{
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<Real>("Dstab", "Dstab");
    return params;
}

DiffusionStab::DiffusionStab(const InputParameters & parameters):
    Kernel(parameters),
    _Dstab(getParam<Real>("Dstab"))
{}

Real DiffusionStab::computeQpResidual()
{
   return 1.0 * _Dstab * _grad_test[_i][_qp]*_grad_u[_qp];
}

Real DiffusionStab::computeQpJacobian()
{
   return 1.0 * _Dstab * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}
