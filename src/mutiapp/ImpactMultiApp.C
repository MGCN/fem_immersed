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
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*         Kernel for the Nearly Incompressibility              */
/****************************************************************/


// MOOSE includes
#include "ImpactMultiApp.h"

#include "AllLocalDofIndicesThread.h"
#include "AuxiliarySystem.h"
#include "Console.h"
#include "MooseMesh.h"
#include "Output.h"
#include "TimeStepper.h"
#include "Transient.h"

// libMesh includes
#include "libmesh/mesh_tools.h"

registerMooseObject("ImmersedBoundaryApp", ImpactMultiApp);

template <>
InputParameters
validParams<ImpactMultiApp>()
{
    InputParameters params = validParams<TransientMultiApp>();
    return params;
}

ImpactMultiApp::ImpactMultiApp(const InputParameters & parameters)
  : TransientMultiApp(parameters)

{
  
}

ImpactMultiApp::~ImpactMultiApp()
{
    if (!_has_an_app)
        return;
    
    MPI_Comm swapped = Moose::swapLibMeshComm(_my_comm);
    

    
    for (unsigned int i = 0; i < _my_num_apps; i++)
    {
        auto & app = _apps[i];
        
        Transient * ex = dynamic_cast<Transient *>(app->getExecutioner());
        
        ex->postExecute();
    }
    
    // Swap back
    Moose::swapLibMeshComm(swapped);
}


