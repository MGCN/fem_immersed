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
/*                Prepared by Barna Becseck                     */  
/*                                                              */  
/****************************************************************/

#ifndef STORETRANSFEROPERATORS_H
#define STORETRANSFEROPERATORS_H

// MOOSE includes
#include "GeneralUserObject.h"

// utopia includes
#include "utopia.hpp"
#include "utopia_fe_base.hpp"
#include <memory>
// Forward Declarations
class StoreTransferOperators;

template<>
InputParameters validParams<StoreTransferOperators>();

class StoreTransferOperators :
  public GeneralUserObject
{
public:
  typedef utopia::USparseMatrix SparseMatT;  

  // constructor
  StoreTransferOperators(const InputParameters & parameters);

  // returns pointer to stored transfer operators
  std::shared_ptr<SparseMatT> &
  getTransferOperator()
  {
    return _B;
  };
  std::shared_ptr<SparseMatT> &
  getTransferOperatorReverse()
  {
    return _B_reverse;
  };

  std::shared_ptr<SparseMatT> &
  setTransferOperator()
  {
    return _B;
  };
  std::shared_ptr<SparseMatT> &
  setTransferOperatorReverse()
  {
    return _B_reverse;
  };
  std::shared_ptr<SparseMatT> &
  getPermuteOperator()
  {
    return _P;
  };

  std::shared_ptr<SparseMatT> &
  setPermuteOperator()
  {
    return _P;
  };

  std::shared_ptr<SparseMatT> &
  getPermuteOperatorDisp()
  {
    return _P_disp;
  };

  std::shared_ptr<SparseMatT> &
  setPermuteOperatorDisp()
  {
    return _P_disp;
  };

	
  std::shared_ptr<void> & getVoidPointer()
  {
        return void_ptr;
  }
 

  /*std::shared_ptr<SparseMatrix<Number>> &
  setFluidMass()
  {
    return _M_fluid;
  };


  std::shared_ptr<SparseMatrix<Number>> &
  getFluidMass()
  {
    return _M_fluid;
  };*/


    std::shared_ptr<void> void_ptr;

  void execute(){};
  void initialize(){};
  void finalize(){};

protected:
  std::shared_ptr<SparseMatT> _B;
  std::shared_ptr<SparseMatT> _B_reverse;
  std::shared_ptr<SparseMatT> _P;
  std::shared_ptr<SparseMatT> _P_disp;
  //std::shared_ptr<SparseMatrix<Number>> _M_fluid;
};

#endif
