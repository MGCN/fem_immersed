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
/*        MECH - ICS Imm Bound  simulation framework            */
/*                Prepared by Maria Nestola                     */
/*                                                              */
/*      Kernel to consider the intertial term.                  */
/*      It implements the term                                  */
/*     _int ((_rho_s - _rho_f) * u_tt * test * dV )             */
/*      by using the Newmark method                             */
/****************************************************************/


#include "InertialTerm.h"
#include "SubProblem.h"
#include "MooseMesh.h"
registerMooseObject("ImmersedBoundaryApp", InertialTerm);
template<>
InputParameters validParams<InertialTerm>()
{
    InputParameters params = validParams<Kernel>();
    params.addParam<std::vector<int>>("boundary_ids", "The list of boundary IDs from the mesh where the density of the fluid will be removed");
    params.addRequiredParam<Real>("rho_f", "density_fluid");
    params.addRequiredParam<Real>("rho_s", "density_solid");
    params.addRequiredCoupledVar("velocity","velocity variable");
    params.addRequiredCoupledVar("acceleration","acceleration variable");
    params.addRequiredParam<Real>("beta","beta parameter for Newmark Time integration");
    params.addRequiredParam<Real>("gamma","gamma parameter for Newmark Time integration");
    params.addParam<Real>("eta",0,"eta parameter for mass dependent Rayleigh damping");
    params.addParam<Real>("alpha_m",0,"alpha parameter for CH method");
    return params;
}

InertialTerm::InertialTerm(const InputParameters & parameters) :
Kernel(parameters),
_rho_f(getParam<Real>("rho_f")),
_rho_s(getParam<Real>("rho_s")),
_u_old(valueOld()),
_vel_old(coupledValueOld("velocity")),
_accel_old(coupledValueOld("acceleration")),
_beta(getParam<Real>("beta")),
_gamma(getParam<Real>("gamma")),
_eta(getParam<Real>("eta")),
_alpha_m(getParam<Real>("alpha_m"))

{
    if (isParamValid("boundary_ids")){
        _boundary_tags=getParam<std::vector<int>>("boundary_ids");
    }
}

Real
InertialTerm::computeQpResidual()
{
        if (_dt == 0)
            return 0;
        else
        {
            Real diff_rho = _rho_s;
            
            
           for(uint nodeno = 0; nodeno < _current_elem->n_nodes(); ++nodeno)
                for(auto t : _boundary_tags)
                    if (_mesh.getMesh().get_boundary_info().has_boundary_id(_current_elem->node_ptr(nodeno),t))
                        diff_rho =_rho_s - _rho_f;
         
            Real accel = ((_u[_qp] - _u_old[_qp])/( _dt * _dt )); // - _vel_old[_qp]/_dt - _accel_old[_qp] * (0.5 -_beta));
            
            Real vel = _vel_old[_qp] + ( _dt * ( 1 - _gamma )) * _accel_old[_qp] + _gamma * _dt * accel;
            
            return _test[_i][_qp] * diff_rho * accel; // + _alpha_m * _test[_i][_qp] * _accel_old[_qp] * diff_rho  + _test[_i][_qp] * diff_rho * vel * _eta;
        }
}
        
Real
InertialTerm::computeQpJacobian()
{

        Real diff_rho = _rho_s;

        for(uint nodeno = 0; nodeno < _current_elem->n_nodes(); ++nodeno)
            for(auto t : _boundary_tags)
                if (_mesh.getMesh().get_boundary_info().has_boundary_id(_current_elem->node_ptr(nodeno),t))
                    diff_rho =_rho_s - _rho_f;

        if (_dt == 0)
            return 0;
        else
            return _test[_i][_qp] *  diff_rho/ (_dt * _dt) * _phi[_j][_qp]; // + _eta * _test[_i][_qp] * diff_rho * _gamma / _beta / _dt * _phi[_j][_qp];
}
