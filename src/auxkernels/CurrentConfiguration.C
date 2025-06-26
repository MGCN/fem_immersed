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
/*         Auxkernel to compute displacement field              */
/*         in the current configuration                         */
/****************************************************************/

#include "CurrentConfiguration.h"
#include "MooseMesh.h"
registerMooseObject("ImmersedBoundaryApp",CurrentConfiguration);
template<>
InputParameters validParams< CurrentConfiguration>()
{
  InputParameters params = validParams<AuxKernel>();
    
  MooseEnum component("disp_xc disp_yc disp_zc");
    
  params.addRequiredParam<MooseEnum>("component", component, "The desired value.");
  params.addRequiredCoupledVar("disp_x", "x-disp");
  params.addCoupledVar("disp_y", "y-disp"); //only required in 2D and in 3D;
  params.addCoupledVar("disp_z", "z-disp"); //only required in 3D;
  return params;
}

 CurrentConfiguration:: CurrentConfiguration(const InputParameters & parameters):
AuxKernel( parameters ),
 _component(getParam<MooseEnum>("component")),
 _disp_x(coupledValue("disp_x")),
 _disp_y(_mesh.dimension()>=2 ? coupledValue("disp_y"): _zero),
 _disp_z(_mesh.dimension()==3 ? coupledValue("disp_z"): _zero)
{}

Real  CurrentConfiguration::computeValue()
{
      Point p = (*_current_node);
    
      Real disp_xc = _disp_x[_qp] + p(0);
      Real disp_yc = _disp_y[_qp] + p(1);  
      Real disp_zc = _disp_z[_qp] + p(2);   

     
   if (_component==0) return disp_xc;  
   if (_component==1) return disp_yc;  
   if (_component==2) return disp_zc;  
      
 std::cout<<"qua non dovresti esserci NODLA POSTPROCESSORAUX\n\n";
    return 0.0;   

    

}


