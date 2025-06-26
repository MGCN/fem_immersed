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
/*              Prepared by Barna Becsek                        */	
/*                                                              */		
/*          Auxkernel to compute acceleration field             */
/****************************************************************/

#include "CorrectedInertia.h"
registerMooseObject("ImmersedBoundaryApp", CorrectedInertia);
template<>
InputParameters validParams<CorrectedInertia>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("inertial_res","inertial_residual");
  params.addRequiredParam<Real>("fluid_density","density of the IBM fluid");
  params.addRequiredParam<Real>("solid_density","density of the IBM solid");
  params.addParam<std::vector<int>>("boundary_ids", 
    "The list of boundary IDs from the mesh where the density of the fluid will be removed");
  return params;
}

CorrectedInertia::CorrectedInertia(const InputParameters & parameters)
  : AuxKernel(parameters),
    _inertial_residual(coupledValue("inertial_res")),
    _fluid_density(getParam<Real>("fluid_density")),
    _solid_density(getParam<Real>("fluid_density"))
{
  if (isParamValid("boundary_ids")){
      _boundary_tags=getParam<std::vector<int>>("boundary_ids");
  }
}

Real
CorrectedInertia::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");

  Real value = _inertial_residual[_qp];
 
  if (isParamValid("boundary_ids")){
    for(auto t : _boundary_tags)
    {
      if (_mesh.getMesh().get_boundary_info().has_boundary_id(_current_node,t))
        value = _inertial_residual[_qp]*(1. - _fluid_density/_solid_density);
    }
  } else {
        value = _inertial_residual[_qp]*(1. - _fluid_density/_solid_density);
  }
  
  return value;
}
