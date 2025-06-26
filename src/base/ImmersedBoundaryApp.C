//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ImmersedBoundaryApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
//#include "MooseUtopia_Solver.h"

// So we can register objects from the fluid_properties module.


template <>
InputParameters
validParams<ImmersedBoundaryApp>()
{
    InputParameters params = validParams<MooseApp>();
    
    return params;
}

registerKnownLabel("ImmersedBoundaryApp");

ImmersedBoundaryApp::ImmersedBoundaryApp(InputParameters parameters) : MooseApp(parameters)
{
    Moose::registerObjects(_factory);
    ImmersedBoundaryApp::registerObjects(_factory);
    
    Moose::associateSyntax(_syntax, _action_factory);
    ImmersedBoundaryApp::associateSyntax(_syntax, _action_factory);
    
}

ImmersedBoundaryApp::~ImmersedBoundaryApp() {}

extern "C" void
ImmersedBoundaryApp_registerApps()
{
    ImmersedBoundaryApp::registerApps();
}
void
ImmersedBoundaryApp::registerApps()
{
    registerApp(ImmersedBoundaryApp);
}

// External entry point for dynamic object registration
extern "C" void
ImmersedBoundaryApp__registerObjects(Factory & factory)
{
    ImmersedBoundaryApp::registerObjects(factory);
}
void
ImmersedBoundaryApp::registerObjects(Factory & factory)
{
    Registry::registerObjectsTo(factory, {"ImmersedBoundaryApp"});
    // PetscErrorCode ierr;  
    // ierr = SNESRegister("stabilized_newmark", SNESCreate_ContactStabilizedNewmark);
    // CHKERRV(ierr);
}

// External entry point for dynamic syntax association
extern "C" void
ImmersedBoundaryApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
    ImmersedBoundaryApp::associateSyntax(syntax, action_factory);
}
void
ImmersedBoundaryApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
    Registry::registerActionsTo(action_factory, {"ImmersedBoundaryApp"});
}


