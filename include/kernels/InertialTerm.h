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
/*      Kernel to consider the intertial term.                  */
/*      It implements the term (_rho_s - _rho_f) * u_tt         */
/*      by using the Newmark method                             */
/****************************************************************/

#ifndef INERTIALTERM_H
#define INERTIALTERM_H

#include "Kernel.h"

//Forward Declarations
class InertialTerm;

template<>
InputParameters validParams<InertialTerm>();

class InertialTerm: public Kernel
{
public:

  InertialTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();
  Real _rho_f;
  Real _rho_s;

private:
  
  const VariableValue & _u_old;
  const VariableValue & _vel_old;
  const VariableValue & _accel_old;
  const Real _beta;
  const Real _gamma;
  const Real _eta;
  const Real _alpha_m;
  std::vector<int> _boundary_tags;
};

#endif //INERTIALTERMH
