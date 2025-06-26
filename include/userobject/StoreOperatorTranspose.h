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
#ifndef StoreOperatorTranspose_H
#define StoreOperatorTranspose_H

#include "GeneralUserObject.h"
#include "libmesh/fparser.hh"
#include "MooseApp.h"
#include "MooseMesh.h"

#include "MultiAppTransfer.h"
#include <utopia_fe.hpp>
#include <utopia_LibMeshBackend.hpp>
#include <utopia.hpp>
#include <memory>
//Forward Declarations
class StoreOperatorTranspose;
class MultiApp;



template<>
InputParameters validParams<StoreOperatorTranspose>();


class StoreOperatorTranspose : public GeneralUserObject
{
public:
    StoreOperatorTranspose(const InputParameters & parameters);
    /// The Terminator DEFINITELY needs a destructor!
    virtual ~StoreOperatorTranspose();
    typedef utopia::USparseMatrix SparseMatT;
    typedef utopia::UVector VecT;
    virtual void initialize() override;
    virtual void execute() override;
    virtual void finalize() override {}
    const  std::shared_ptr<SparseMatT> & computeOpT1(){ return _B;}
    const  std::shared_ptr<SparseMatT> & computeOpT2(){ return _B;}
    
protected:
    
    
    FEProblem * _fe_problem;
    
    /// The name of the variable to transfer to
    
    std::string _to_var_name;
    
    /// Name of variable transfering from
    std::string  _from_var_name;
    
    bool _impact, _biorthogonal;
    
    bool _fixed_meshes;
    
    bool _moved_source_mesh;
    
    bool _moved_target_mesh;
    
    std::string _multiapp_name;
    
    bool _cache_matrices = false;
    
    bool _pseudoL2 = false;
    
    bool _displaced_source_mesh;
    
    bool _displaced_target_mesh;
    
    //enum _direction;

    enum
    {
        TO_MULTIAPP,
        FROM_MULTIAPP
    };
    
   // std::shared_ptr<MultiApp>  _multi_app;
    /// The MultiApp this Transfer is transferring data to or from
    //std::shared_ptr<MultiApp> _multi_app_name
    // std::shared_ptr<MultiApp> * _multi_app;
    // FEProblem * _fe_problem;
    std::shared_ptr<SparseMatT> _B;
    
    void assemble_mass_matrix(EquationSystems & es,
                              const std::string & system_name, SparseMatT &D_inv);
    
};

#endif //TERMINATOR_H
