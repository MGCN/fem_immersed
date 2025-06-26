/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef IMMERSEDBOUNDARYAPP_H
#define IMMERSEDBOUNDARYAPP_H

#include "MooseApp.h"

class ImmersedBoundaryApp;

template<>
InputParameters validParams<ImmersedBoundaryApp>();

class ImmersedBoundaryApp : public MooseApp
{
public:
  ImmersedBoundaryApp(InputParameters parameters);
  virtual ~ImmersedBoundaryApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* IMMERSEDBOUNDARYAPP_H */
