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

#include "FluidDisplacementAux.h"
#include "MooseMesh.h"
registerMooseObject("ImmersedBoundaryApp", FluidDisplacementAux);
template<>
InputParameters validParams<FluidDisplacementAux>()
{ 

  InputParameters params = validParams<AuxKernel>();

  MooseEnum component("disp_xf disp_yf disp_zf");
  params.addRequiredParam<MooseEnum>("component", component, "The desired value.");
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); //only required in 2D and in 3D;
  params.addCoupledVar("w", "z-velocity"); //only required in 3D;
 
  
  return params;
}

FluidDisplacementAux::FluidDisplacementAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    
    
    _u_vel(coupledValue("u")),
    _v_vel(_mesh.dimension()>=2 ? coupledValue("v"): _zero),
    _w_vel(_mesh.dimension()==3 ? coupledValue("w"): _zero),
    _component(getParam<MooseEnum>("component"))
{}

Real FluidDisplacementAux::computeValue()
{    
      Real disp_xf=_u_vel[_qp]*_dt;
      Real disp_yf=_v_vel[_qp]*_dt;  
      Real disp_zf=_w_vel[_qp]*_dt;    

   if (_component==0)  return disp_xf;  
   if (_component==1)  return disp_yf;  
   if (_component==2)  return disp_zf;  
      
 std::cout<<"qua non dovresti esserci NODLA POSTPROCESSORAUX\n\n";
    return 0.0;
   
}

