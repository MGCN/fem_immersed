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
#ifndef ImpactMultiApp_H
#define ImpactMultiApp_H

#include "MultiApp.h"

#include "TransientMultiApp.h"

// Forward declarations
class ImpactMultiApp;
class Transient;

template <>
InputParameters validParams<ImpactMultiApp>();


class ImpactMultiApp : public TransientMultiApp
{
public:
    ImpactMultiApp(const InputParameters & parameters);
    virtual ~ ImpactMultiApp();
    
private:
    
    std::vector<Transient *> _transient_executioners;
    
};

#endif // ImpactMultiApp_H
