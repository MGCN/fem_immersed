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
/*    Immersed_Boundary-ICS Mechanical simulation framework    */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/


#ifndef IBTransfer_H
#define IBTransfer_H

#include "MultiAppTransfer.h"
#include <utopia_fe.hpp>
#include <utopia_LibMeshBackend.hpp>
#include <utopia.hpp>

#include "UserObjectInterface.h"
#include "StoreTransferOperators.h"
#include <memory>
#include <unordered_set>
class IBTransfer;


 // *                                  _.._
 // *                                .' .-'`
 // *                               /  /
 // *                               |  |
 // *                               \  '.___.;
 // *                                '._  _.'
 // *
 // *                 ) \     /_(
 // *                  )_`-)-_(_
 // *                   `,' `__,)
 // *                  _/   ((
 // *         ______,-'    )
 // *        (            ,
 // *         \  )    (   |
 // *        /| /`-----` /|
 // *        \ \        / |
 // *        |\|\      /| |\


template <>
InputParameters validParams<IBTransfer>();

/**
 * Project values from one domain to another
 */
class IBTransfer : public MultiAppTransfer,
   public UserObjectInterface
{
public:
  IBTransfer(const InputParameters & parameters);

    virtual void initialSetup();
  
    virtual void execute();

    virtual void assembleSys(EquationSystems & es,const std::string & system_name);  
    
    typedef utopia::USparseMatrix SparseMatT;
  
    typedef utopia::UVector VecT;

protected:
    AuxVariableName _to_var_name_1;
  
    VariableName _from_var_name_1;

    AuxVariableName _to_var_name_2;
  
    VariableName _from_var_name_2;

    AuxVariableName _to_var_name_3;
  
    VariableName _from_var_name_3;

    std::string _to_var_name_r;
  
    std::string _from_var_name_r;

    //AuxVariableName _to_var_name_main;
    
    //AuxVariableName _from_var_name_main_1;

    MooseEnum _proj_type;

    bool _pseudoL2 = false;

    bool _impactEFSI, _reverse;

    bool _compute_operators;
    
    unsigned int _n_var_r;

    unsigned int _n_var;
    
    //AuxVariableName _from_var_name_main_2, _from_var_name_main_3;
    
    // userobject to store our operators
    const StoreTransferOperators & _operator_storage;
 
    unsigned int _recompute_interval;
  
    Real _sys_tol;

    bool _displaced_source_mesh, _displaced_target_mesh;

    Real _scale_factor;    

    LinearImplicitSystem * _proj_sys;
  
    unsigned int _proj_var_num1,_proj_var_num2,_proj_var_num3;
    
    std::shared_ptr<SparseMatT> _B = NULL;
    std::shared_ptr<SparseMatT> _B_reverse = NULL;
    std::shared_ptr<SparseMatT> _P = NULL;
    std::shared_ptr<SparseMatT> _P_efsi = NULL;
 
    EquationSystems * _equationSystemsT;
    LinearImplicitSystem * _linearImplicitSystemT; 


    unsigned int _test_var_1;
    unsigned int _test_var_2;
    unsigned int _test_var_3;
    void buildTransferOperators();

    virtual void get_solution(FEProblemBase * problem, std::string var_name, Vec & vector);
 
    virtual void set_pseudo_solution(VecT &sol);
  
    //virtual void assembleIBSys(EquationSystems & es,const std::string & system_name);
  
    friend void  assemblyIBVector(EquationSystems & es, const std::string & system_name);
  
    virtual void projectSolution();
  
    virtual void set_solution(FEProblemBase &problem, EquationSystems &proj_es);

    virtual void Pseudol2projection();
    
    // dofs participating in the transfer
    std::unordered_set<utopia::SizeType> _participating_slave_dofs;

    void Permutation(FEProblemBase *problem, std::shared_ptr<SparseMatT> P_matrix, unsigned int number);

    void PermutationEFSI(FEProblemBase *problem, std::shared_ptr<SparseMatT> P_matrix, unsigned int number);

    bool _first_recover = true;

    bool _compute_FSI=true;
};

#endif /* IBTransfer_H */
