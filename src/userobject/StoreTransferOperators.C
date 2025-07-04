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
/*    Immersed_Boundary - ICS Mechanical simulation framework   */		
/*                Prepared by Maria Nestola,                    */	
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */	
/****************************************************************/

#include "StoreTransferOperators.h"
#include "MooseVariable.h"

registerMooseObject("ImmersedBoundaryApp", StoreTransferOperators);

template<>
InputParameters validParams<StoreTransferOperators>()
{
  InputParameters params = validParams<GeneralUserObject>();
  return params;
}

StoreTransferOperators::StoreTransferOperators(const InputParameters & parameters) :
    GeneralUserObject(parameters)
{
   _B           = std::make_shared<SparseMatT>();
   _B_reverse   = std::make_shared<SparseMatT>();
   _P           = std::make_shared<SparseMatT>();
   _P_disp      = std::make_shared<SparseMatT>();
  // _M_fluid     = std::make_shared<SparseMatrix<Number>>();
}
