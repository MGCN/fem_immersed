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


#ifndef WEAK_CONSTRAINT_MULTIPLIER
#define WEAK_CONSTRAINT_MULTIPLIER


#include "Kernel.h"


class Weak_Constraint_multiplier;
class MooseMesh;
class FEProblem;
template<>
InputParameters validParams<Weak_Constraint_multiplier>();

class Weak_Constraint_multiplier : public Kernel
{

public:
  Weak_Constraint_multiplier(const InputParameters & parameters);


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
  private:
//   FEProblem & _fe_problem;
//   FEProblem * _problem_ptr;
//   MooseMesh * _mesh_ptr;
  const VariableValue& _lambda;
  unsigned int _lambda_var_number;
//   MooseMesh * _mesh;
 
//  MeshBase & getMesh();
};
 
#endif 