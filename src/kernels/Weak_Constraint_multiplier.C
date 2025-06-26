/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

/****************************************************************/	
/*               DO NOT MODIFY THIS HEADER                      */	
/*       MECH - ICS Mechanical simulation framework             */	
/*                Prepared by Maria Nestola,                    */	
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */	
/*     Kernel to impose an essential conditions weakly.         */
/*     To be used with Weak_Constraint.                         */ 
/*     It implements the term _int(_lambda * _test_function )dS */
/****************************************************************/



#include "Weak_Constraint_multiplier.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"
#include "DisplacedProblem.h"
#include "GeneratedMesh.h"
#include "MooseMesh.h"
#include "libmesh/mesh_tools.h"

registerMooseObject("ImmersedBoundaryApp", Weak_Constraint_multiplier);

template<>
InputParameters validParams<Weak_Constraint_multiplier>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("lambda", "lagrange multiplier");
  return params;      
  }

Weak_Constraint_multiplier::Weak_Constraint_multiplier(const InputParameters & parameters) :
    Kernel(parameters),
    _lambda(coupledValue("lambda")),
    _lambda_var_number(coupled("lambda"))
{}  


Real
Weak_Constraint_multiplier::computeQpResidual()
{ 
  return   _lambda[_qp]*_test[_i][_qp];

}      
    
Real
Weak_Constraint_multiplier::computeQpJacobian()
{ 

  return 0.0;
    
}

Real
Weak_Constraint_multiplier::computeQpOffDiagJacobian(unsigned int j_var)
{        
    if (j_var == _lambda_var_number)
    
          return  _phi[_j][_qp]*_test[_i][_qp]; 
    else
          return 0.0;
          
}


