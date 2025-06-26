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
#include "NodalSumR.h"
#include "MooseMesh.h"
#include "SubProblem.h"
#include "Assembly.h"

registerMooseObject("ImmersedBoundaryApp",NodalSumR);

template<>
InputParameters validParams<NodalSumR>()
{
  InputParameters params = validParams<NodalVariablePostprocessor>();
  params.set<bool>("unique_node_execute") = true;
  return params;
}

NodalSumR::NodalSumR(const InputParameters & parameters) :
    NodalVariablePostprocessor(parameters),
    _mesh(_subproblem.mesh().getMesh()),
    _volume(0),
    _current_side_volume(_assembly.sideElemVolume()),
    _sum(0)
{
}

void
NodalSumR::initialize()
{
  _sum = 0;
  _volume = 0;
}

void
NodalSumR::execute()
{
  _sum += _u[_qp];
  _volume += 1;
}

Real
NodalSumR::getValue()
{
  gatherSum(_sum);
  
  gatherSum(_volume);
    
    BoundaryInfo & bi = _mesh. get_boundary_info();
    
    std::vector<dof_id_type> node_list;
    
    std::vector<boundary_id_type> bc_id_list;
    
    bi.build_node_list (node_list, bc_id_list);
    
    std::cout<<"bc_id_list"<<_volume<<std::endl;
    

  return _sum/_volume;
}

void
NodalSumR::threadJoin(const UserObject & y)
{

    
  const NodalSumR & pps = static_cast<const NodalSumR &>(y);
    
  _sum += pps._sum;
}
