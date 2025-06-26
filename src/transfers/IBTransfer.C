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
/*    Immersed_Boundary- ICS/ARTORG Mechanical                  */
/*                simulation framework                          */
/*           Prepared by Maria Nestola, Barna Becsek            */
/*                 ARTORG, UniBe, 3008 Bern                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/


#include "IBTransfer.h"
#include "FEProblem.h"
#include "AddVariableAction.h"
#include "MooseMesh.h"
#include "AuxiliarySystem.h"
#include "NonlinearSystem.h"
#include "DisplacedProblem.h"
#include <opencl_adapter.hpp>
#include "libmesh/linear_partitioner.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/system_subset_by_subdomain.h"
//#include <LibmeshTransferForMooseReverse.hpp>
#include <petscksp.h>
#include <petscmat.h>
#include "utopia_TransferAssembler.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include <utopia.hpp>

#include "MooseVariable.h"
#include "Transient.h"
#include "MultiApp.h"
#include <cmath>
#include <queue>
#include <unordered_set>

#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_solver_type.h"


typedef utopia::USparseMatrix SparseMatT;
typedef utopia::UVector VecT;


//registerMooseObject("ImmersedBoundaryApp", IBTransfer);
registerMooseObject("ImmersedBoundaryApp", IBTransfer);

void
assemblyVector(EquationSystems & es, const std::string & system_name)
{
    IBTransfer * transferObj =
    es.parameters.get< IBTransfer *>("transfer");
    transferObj->assembleSys(es, system_name);
}

template <>
InputParameters
validParams<IBTransfer>()
{
    MooseEnum proj_type("L2 pseudo-L2");
    
    InputParameters params = validParams<MultiAppTransfer>();
    params.addRequiredParam<AuxVariableName>("variable_1",
                                             "The auxiliary variable to store the transferred values in.");
    
    params.addRequiredParam<VariableName>("source_variable_1",
                                          "The variable to transfer from.");
    
    params.addParam<AuxVariableName>("variable_2",
                                     "The auxiliary variable to store the transferred values in.");
    
    params.addParam<VariableName>("source_variable_2",
                                  "The variable to transfer from.");
    
    params.addParam<AuxVariableName>("variable_3",
                                     "The auxiliary variable to store the transferred values in.");
    
    params.addParam<VariableName>("source_variable_3",
                                  "The variable to transfer from.");
    
    
    params.addParam<AuxVariableName>("variable_r",
                                     "The auxiliary variable to store the transferred values in.");
    
    params.addParam<VariableName>("source_variable_r",
                                  "The variable to transfer from.");
    
    params.addParam<MooseEnum>("proj_type", proj_type,
                               "The type of the projection.");
    
    params.addParam<bool>("pseudoL2",false,
                          "Set to false when you want to pure L2-projection.");
    
    params.addParam<bool>("impactEFSI",false,"impactEFSI");
    
    params.addParam<bool>("reverse",false,
                          "reverse matrix.");
    
    params.addParam<bool>("displaced_source_mesh", false,
                          "Use the source mesh after displacement");
    
    params.addParam<bool>("displaced_target_mesh", false,
                          "Use the target mesh after displacement");
    
    params.addParam<bool>("compute_operators",true,
                          "if true, this class will assemble the operators (including the reverse), "
                          "if false the operator is read from a userobject");
    
    params.addRequiredParam<unsigned int>("n_var_r","number of reversible variable");
    
    params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    
    params.addRequiredParam<unsigned int>("recompute_interval","the interval of timesteps after which to recompute the transfer operators");

    params.addRequiredParam<Real>("sys_tol","tolerance to which to solve the system to.");    

    params.addParam<Real>("system_scale_factor", 1.,
                          "factor to multiply the system with (for better convergence behavior)");

    return params;
    
}

IBTransfer::IBTransfer(const InputParameters & parameters)
:  MultiAppTransfer(parameters), UserObjectInterface(this),
_to_var_name_1(getParam<AuxVariableName>("variable_1")),
_from_var_name_1(getParam<VariableName>("source_variable_1")),
_to_var_name_r(getParam<AuxVariableName>("variable_r")),
_from_var_name_r(getParam<VariableName>("source_variable_r")),
_proj_type(getParam<MooseEnum>("proj_type")),
_pseudoL2(getParam<bool>("pseudoL2")),
_impactEFSI(getParam<bool>("impactEFSI")),
_reverse(getParam<bool>("reverse")),
_compute_operators(getParam<bool>("compute_operators")),
_n_var_r(getParam<unsigned int>("n_var_r")),
_operator_storage(getUserObject<StoreTransferOperators>("operator_userobject")),
_recompute_interval(getParam<unsigned int>("recompute_interval")),
_sys_tol(getParam<Real>("sys_tol")),
_displaced_source_mesh(getParam<bool>("displaced_source_mesh")),
_displaced_target_mesh(getParam<bool>("displaced_target_mesh")),
_scale_factor(getParam<Real>("system_scale_factor"))
{
    _n_var=1;
    
    if (isParamValid("variable_2")){
        _to_var_name_2 = getParam<AuxVariableName>("variable_2");
        _from_var_name_2 = getParam<VariableName>("source_variable_2");
        _n_var++;}
    
    if (isParamValid("variable_3")){
        _to_var_name_3 = getParam<AuxVariableName>("variable_3");
        _from_var_name_3 = getParam<VariableName>("source_variable_3");
        _n_var++;}
    
//    std::cout<<"n_var"<<_n_var<<std::endl;
   
    _B = NULL;
    _B_reverse = NULL;
    
}




void
IBTransfer::initialSetup()
{
    _console << "Initial Setup of Transfer System!" << std::endl;
   
    getAppInfo();
    
    if (_to_problems.size() > 1 || _from_problems.size() > 1) {
        mooseError("NOT IMPLEMENTED: Right now we can only handle one master and one slave");
    }
    
   
     if (!_app.isRecovering()){ 

	     _first_recover = false;

     }
   
   
    unsigned int zero = 0;
    
    if(_pseudoL2==false){
        
        //_console << "L2-projection" << std::flush << std::endl;
        
        FEProblemBase & to_problem = *_to_problems[zero];
        
        MeshBase *to_mesh = &to_problem.mesh().getMesh();
        
        int  to_dim=to_mesh->mesh_dimension();
        
        EquationSystems & to_es = to_problem.es();
        
        // Add the projection system.
        
        FEType fe_type = to_problem.getStandardVariable(0, _to_var_name_1).feType();
        _proj_sys = &to_es.add_system<LinearImplicitSystem>("proj-sys-" + name());
        _proj_var_num1 = _proj_sys->add_variable("var1", fe_type);

	if (_n_var>1) {
            fe_type = to_problem.getStandardVariable(0, _to_var_name_2).feType();
            _proj_var_num2 = _proj_sys->add_variable("var2", fe_type); 
        }
        if (_n_var>2) { 
            fe_type = to_problem.getStandardVariable(0, _to_var_name_3).feType();
            _proj_var_num3 = _proj_sys->add_variable("var3", fe_type); 
        }

	_proj_sys->attach_assemble_function(assemblyVector);
        _proj_sys->hide_output() = true;
   
   	// when we transfer to IMPACT we do not want to zero out the operator matrix
	// since it is unchanged over the simulation and we want to cache it.
        // we'll take care of zeroing out the rhs manually during the assembly
	_proj_sys->zero_out_matrix_and_rhs = false;
	std::cout<<"n_var"<<_n_var<<std::endl;
        to_es.reinit();
    }
    
}


void
IBTransfer::execute()
{
    
    //_console << "L2ProjectionLibMeshTransfer transfer " << std::endl;
    
    
    IBTransfer::buildTransferOperators();
    
    
    if(_pseudoL2==false){
        
      //_console << "Projecting solution." << std::endl;
        
      IBTransfer::projectSolution();
        
    }
    
    else{
          
      IBTransfer::Pseudol2projection();
    }
    
    
    //_console << "L2ProjectionLibMeshTransfer " << name() << std::endl;
    
}

void
IBTransfer::buildTransferOperators()
{
    using namespace utopia;
    //_console << " BuildTransferOperators for all the variables " << std::endl;
    
    
    FEProblemBase & from_problem = *_from_problems[0];
    
    MeshBase *from_mesh = &from_problem.mesh().getMesh();
    
    int  from_dim=from_mesh->mesh_dimension();
    
    // displaced meshes are special
    //
    //std::cout<<"from_problem.getDisplacedProblem()"<<from_problem.getDisplacedProblem()<<std::endl;


    if (_displaced_source_mesh && from_problem.getDisplacedProblem())
    {
        from_mesh = &from_problem.getDisplacedProblem()->mesh().getMesh();
        
        _console << "_displaced_source_mesh"<<std::endl;
    }
    
    FEProblemBase & to_problem = *_to_problems[0];
    
    MeshBase * to_mesh = &to_problem.mesh().getMesh();
    
    int  to_dim=to_mesh->mesh_dimension();
    
    if (_displaced_target_mesh && to_problem.getDisplacedProblem())
    {
        
        to_mesh = &to_problem.getDisplacedProblem()->mesh().getMesh();
        
        _console << "_displaced_target_mesh"<<std::endl;
        
    }
    
    
    moonolith::Communicator expressComm(_communicator.get());
    
    MooseVariable & _from_var_1 = _from_problems[0]->getStandardVariable(0, _from_var_name_1);
    
    MooseVariable & _to_var_1 = _to_problems[0]->getStandardVariable(0, _to_var_name_1);
    
    // check for correct variable numbers, so that we can use them for the reverse operators.
    // Essentially, the variables we transfer from and to ought to have the same number in the
    // system respectively, such that rows and columns of the operator match also in the reversed
    // state. For simplicity, we require they be the first variables. It doesn't check for the
    // other components but for now let's assume that those come consecutively.
    // if (_to_var_1.number() != 0 || _from_var_1.number() != 0)
    // mooseError("Reorder your variables and auxvariables for the transfer such that they start at 0");
    
    
    // monitor Picard iterations and only recompute operators once we advance a timestep
    Real picIts = static_cast<Transient*>(getMooseApp().getExecutioner())->picardSolve().numPicardIts();
    
    if( ( (_compute_operators && !( (from_problem.timeStep()-1) % _recompute_interval) && !(picIts-1) ) || !_B || !_B_reverse ) && !_reverse)
    {
        // this is necessary after a restart where _recompute_interval is not equal 1
        if (!_B || !_B_reverse) {
          _console << "Assembling operators due to them being non-existent \n";
        }
        // dof sizes
        // this is a main system
         //DofMap &master_dof = _from_var_1.sys().system().get_dof_map(); //solid_problem
        // this is an auxiliary system
      
       	 DofMap &slave_dof  = _to_var_1.sys().system().get_dof_map();  //fluid_problem
        
        //MooseVariable & _master_aux = _from_problems[0]->getStandardVariable(0, _to_var_name_r);
        //MooseVariable & _slave_main = _to_problems[0]->getStandardVariable(0, _from_var_name_r);
        
        //DofMap &slave_dof_r = _slave_main.sys().system().get_dof_map(); //solid_problem
        //DofMap &master_dof_r = _master_aux.sys().system().get_dof_map();  //fluid_problem
        
       _console << "Assembling operators with MOONoLith\n";
       _B            = std::make_shared<SparseMatT>();
       _B_reverse    = std::make_shared<SparseMatT>();
        
       utopia::InputParameters params;
      
       TransferOptions opts;

       opts.n_var  = _n_var;

       //MooseVariable & _fr_v = _from_problems[0]->getStandardVariable(0, "res_elastic_scaled_sub_x");

       //MooseVariable & _to_v = _to_problems[0]->getStandardVariable(0, "residual_comp_x");

       MooseVariable & _from_var_d = _from_problems[0]->getStandardVariable(0,_from_var_name_r);    /*disp*/

       DofMap &master_dof = _from_var_d.sys().system().get_dof_map(); //solid_problem

       opts.from_var_num     = _from_var_d.number(); //master_system.variable_number("var_1_s"); //_master_aux.number();

       opts.to_var_num       = _to_var_1.number();

       //from_mesh->print_info();

       //to_mesh->print_info();

       //std::cout<<"number of local dofs"<<slave_dof.n_dofs()<<std::endl;
	
       auto cs = std::make_shared<utopia::MeshTransferOperator>(make_ref(*from_mesh),
                                                                 make_ref(master_dof),
                                                                 make_ref(*to_mesh),
                                                                 make_ref(slave_dof),
                                                                 opts);


        //std::unordered_set<utopia::SizeType> &particip_slave_dofs,

        params.set("coupling-matrix-only", true);
        cs->read(params);
        cs->initialize(BIDIRECTIONAL_L2_PROJECTION);
        auto &matrices = cs->matrices();
        assert(matrices.size() == 2);
        _B = matrices[0];
        _B_reverse = matrices[1];
       
       	_console << "Saving operators in UserObject\n";
        // cast away constantness to be able to access members
        
	_participating_slave_dofs.clear();
        //utopia::Write<USparseMatrix> w_B(*_B,utopia::GLOBAL_INSERT);
        utopia::each_read(*_B, [&](const utopia::SizeType i, const utopia::SizeType j, const double value) {
            //for(utopia::SizeType d = 0; d < _n_var; ++d) {
                //(*_B).set(i+d, j+d, value);
                _participating_slave_dofs.emplace(i); // here we collect all dofs of the slave grid to be able to
                                                  // restrict our mass matrix
           // }
        });
	
	const_cast<StoreTransferOperators&>(_operator_storage).setTransferOperator()        = _B;
        const_cast<StoreTransferOperators&>(_operator_storage).setTransferOperatorReverse() = _B_reverse;
        // _console << "operator sizes\n";
        // disp((*_B).size());
        // disp((*_B_reverse).size());
        // SparseMatT B_t = transpose(*_B);
        // write("B_reverse.m",*_B_reverse);
        // write("B_rev.m",*_B_reverse);
      
        
    } else
    {
        _console << "Reading operator from UserObject\n";
        _B            = std::make_shared<SparseMatT>();
        _B_reverse    = std::make_shared<SparseMatT>();
       
       	if (_reverse)
            _B = const_cast<StoreTransferOperators&>(_operator_storage).getTransferOperatorReverse();
        
        else     
            _B = const_cast<StoreTransferOperators&>(_operator_storage).getTransferOperator();
    }



    if( (_reverse && _direction == FROM_MULTIAPP) && (_to_problems[0]->timeStep()== 1 || _compute_FSI)){


          _P = std::make_shared<SparseMatT>();

          _console << "Compute FSI Permutation" << _first_recover << std::flush << std::endl;

          Permutation(_to_problems[0], _P, 0);

          const_cast<StoreTransferOperators&>(_operator_storage).setPermuteOperator() = _P;

	  _compute_FSI=false;
    }

    if((_direction == TO_MULTIAPP) && (_to_problems[0]->timeStep()== 1 || _first_recover)){
	 
          
            _console << "Adding Mass Slave System" << std::flush << std::endl;
	 
	    _equationSystemsT=&_to_problems[0]->es();
	
	    _linearImplicitSystemT = &_equationSystemsT->add_system<LinearImplicitSystem> ("MassFluidSys");
	  
	    _test_var_1 = _linearImplicitSystemT[0].add_variable ("test_var_1", FIRST);
	  
	    if(_n_var > 1) {
	    
		    _test_var_2 = _linearImplicitSystemT[0].add_variable ("test_var_2", FIRST);
	    }

	  
	    if(_n_var > 2){
               
		    _test_var_3 = _linearImplicitSystemT[0].add_variable ("test_var_3", FIRST);
	    }
      }

     if((_direction == TO_MULTIAPP) && (_to_problems[0]->timeStep()== 1 || _first_recover)){

             _console << "Adding Fluid Mass Matrix" << std::flush << std::endl;

	     _linearImplicitSystemT[0].add_matrix("FluidMassMatrix");
	    
	    _equationSystemsT[0].reinit();
     }

   

    if(_impactEFSI && _direction == TO_MULTIAPP){
	    
	    if(_from_problems[0]->timeStep()== 1 || _first_recover){

            //_P_1 = std::make_shared<SparseMatT>();

             _P_efsi = std::make_shared<SparseMatT>();

             _console << "Compute Permutation EFSI" << std::flush << std::endl;

             PermutationEFSI(_from_problems[0], _P_efsi, 0);

              const_cast<StoreTransferOperators&>(_operator_storage).setPermuteOperatorDisp() = _P_efsi;
           }
      }
    
// already called by execute()   
//    if(_pseudoL2)
//    {
//        Pseudol2projection();
//    }
    
}

void
IBTransfer::Pseudol2projection()
{
    using namespace utopia;
    
    
    _console << "Pseudol2projection" << std::flush << std::endl;
    
    Vec tmp;
    
    UVector  diag_elem;
    
    VecT from_solution, solution;
    
    Vec w;
    
    get_solution(_from_problems[0], _from_var_name_1, tmp);
    
    utopia::USparseMatrix T;
    
    utopia::convert(tmp, from_solution);
        
    utopia::USparseMatrix Dinv;
    
    utopia::USparseMatrix D;
    
     if(_reverse)
     {
         _console << "Pseudol2projection Reverse" << std::flush << std::endl;
        
         D = diag(sum(*_B_reverse, 1));

         diag_elem = 1./sum((*_B_reverse),1);
        
         Dinv = diag(diag_elem);
        
         T = Dinv*(*_B_reverse);
     }
    else{
        
        D = diag(sum(*_B, 1));

        diag_elem = 1./sum((*_B),1);
        
        Dinv = diag(diag_elem);
        
        T = Dinv*(*_B);
    }
    
    
    solution = T * from_solution;
    
    set_pseudo_solution(solution);
    
    
}



void
IBTransfer::get_solution(FEProblemBase * problem, std::string var_name,  Vec &vector)
{
    
    //_console << "Get Solution Vector" << std::flush << std::endl;
    
    MooseVariable & var = problem->getStandardVariable(0, var_name);
    
    System & sys = var.sys().system();
    
    vector = static_cast<PetscVector<Number> *>(sys.current_local_solution.get())->vec();
    
}



void
IBTransfer::set_pseudo_solution(VecT &sol)
{
    
    _console << "Set PSEUDO-L2-projection Solution Vector" << std::flush << std::endl;
    
    FEProblemBase & to_problem = *_to_problems[0];
    MeshBase & to_mesh = _to_problems[0]->mesh();
    
    MooseVariable & slave_var_1 = to_problem.getStandardVariable(0, _to_var_name_1);
    
    System & slave_sys = slave_var_1.sys().system();
    NumericVector<Number> * to_solution = slave_sys.solution.get();
    
    PetscInt       rstart,rend;
    PetscScalar    tmp_sol;
    VecGetOwnershipRange(utopia::raw_type(sol),&rstart,&rend);
    
    
    
    {
        MeshBase::const_node_iterator it = to_mesh.local_nodes_begin();
        const MeshBase::const_node_iterator end_it = to_mesh.local_nodes_end();
        for ( ; it != end_it; ++it)
        {
            const Node * node = *it;
            
            for (unsigned int comp = 0;comp < node->n_comp(slave_sys.number(), slave_var_1.number()); comp++)
            {
                const dof_id_type to_index = node->dof_number(slave_sys.number(), slave_var_1.number(), comp);
                
                to_solution->set(to_index,  sol.get(to_index));
                
            }
        }
    }
    
    if (_n_var == 2)
    {
        
        MooseVariable & slave_var_2 = to_problem.getStandardVariable(0, _to_var_name_2);
        
        MeshBase::const_node_iterator it = to_mesh.local_nodes_begin();
        
        const MeshBase::const_node_iterator end_it = to_mesh.local_nodes_end();
        for ( ; it != end_it; ++it)
        {
            const Node * node = *it;
            
            for (unsigned int comp = 0;comp < node->n_comp(slave_sys.number(), slave_var_2.number()); comp++)
            {
                const dof_id_type to_index = node->dof_number(slave_sys.number(), slave_var_2.number(), comp);
                
                to_solution->set(to_index,  sol.get(to_index));
                
            }
        }
    }
    
    if (_n_var == 3)
    {
        
        MooseVariable & slave_var_3 = to_problem.getStandardVariable(0, _to_var_name_3);
        
        MeshBase::const_node_iterator it = to_mesh.local_nodes_begin();
        
        const MeshBase::const_node_iterator end_it = to_mesh.local_nodes_end();
        for ( ; it != end_it; ++it)
        {
            const Node * node = *it;
            
            for (unsigned int comp = 0;comp < node->n_comp(slave_sys.number(), slave_var_3.number()); comp++)
            {
                const dof_id_type to_index = node->dof_number(slave_sys.number(), slave_var_3.number(), comp);
                
                to_solution->set(to_index,  sol.get(to_index));
                
            }
        }
    }
    
    
    to_solution->close();
    slave_sys.update();
    
    {
        MeshBase::const_element_iterator it = to_mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator end_it = to_mesh.active_local_elements_end();
        for ( ; it != end_it; ++it)
        {
            const Elem * elem = *it;
            for (unsigned int comp = 0; comp < elem->n_comp(slave_sys.number(), slave_var_1.number()); comp++)
            {
                const dof_id_type to_index = elem->dof_number(slave_sys.number(), slave_var_1.number(), comp);
                to_solution->set(to_index, sol.get(to_index));
            }
        }
    }
    
    
    if (_n_var == 2)
    {
        
        MooseVariable & slave_var_2 = to_problem.getStandardVariable(0, _to_var_name_2);
        
        MeshBase::const_element_iterator it = to_mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator end_it = to_mesh.active_local_elements_end();
        for ( ; it != end_it; ++it)
        {
            const Elem * elem = *it;
            for (unsigned int comp = 0; comp < elem->n_comp(slave_sys.number(), slave_var_2.number()); comp++)
            {
                const dof_id_type to_index = elem->dof_number(slave_sys.number(), slave_var_2.number(), comp);
                to_solution->set(to_index, sol.get(to_index));
            }
        }
    }
    
    if (_n_var == 3)
    {
        
        MooseVariable & slave_var_3 = to_problem.getStandardVariable(0, _to_var_name_3);
        
        MeshBase::const_element_iterator it = to_mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator end_it = to_mesh.active_local_elements_end();
        for ( ; it != end_it; ++it)
        {
            const Elem * elem = *it;
            for (unsigned int comp = 0; comp < elem->n_comp(slave_sys.number(), slave_var_3.number()); comp++)
            {
                const dof_id_type to_index = elem->dof_number(slave_sys.number(), slave_var_3.number(), comp);
                to_solution->set(to_index, sol.get(to_index));

            }
        }
    }
    
    
    to_solution->close();
    slave_sys.update();
    
}



void
IBTransfer::assembleSys(EquationSystems & es,
                                               const std::string & system_name)
{
    _console << "Assembling the EquationSystem for L2-pojection." << std::flush << std::endl;
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to our system.
    LinearImplicitSystem & slave_system = es.get_system<LinearImplicitSystem>(system_name);
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = slave_system.get_dof_map().variable_type(0);
    
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    
    //SparseMatrix<Number> & matrix_A = *slave_system.matrix;
    auto & matrix_A = *slave_system.matrix;


    //_console << "is matrix closed: " << matrix_A.closed() << std::endl;

    // The element mass matrix.
    DenseMatrix<Number> Me;
    
    // A  Gauss quadrature rule for numerical integration.
    QGauss qrule (dim, fe_type.default_quadrature_order());
    
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);
    
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();
    
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    
    const DofMap & dof_map = slave_system.get_dof_map();
    
    std::vector<dof_id_type> dof_indices;
 
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    // monitor Picard iterations and only recompute operators once we advance a timestep
    Real picIts = static_cast<Transient*>(getMooseApp().getExecutioner())->picardSolve().numPicardIts();

    // --- Set Up System Matrix --- 
    // loop over elements in the slave system 
    // this is relatively expensive, especially for impact where there are many elements!
    // we could cache this matrix, since impact does not deform its grid.
    if (_direction == FROM_MULTIAPP && 
       ((_direction == FROM_MULTIAPP && !matrix_A.closed()) || 
       (_direction == FROM_MULTIAPP && !(((*_to_problems[0]).timeStep()-1) % _recompute_interval) && !(picIts-1)) ) ) {
            _console << "Assembling Solid System Matrix Begin" << std::endl;
            // first we need to manually zero the matrix
            matrix_A.zero();   
	    
	    for ( ; el != end_el; ++el)
	    {
		const Elem * elem = *el;
		Elem * ele = *el;
		
		dof_map.dof_indices (elem, dof_indices);
		
		fe->reinit (elem);
		
		Me.resize (dof_indices.size(), dof_indices.size());
		
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		    for (unsigned int i=0; i<phi.size(); i++)
			for (unsigned int j=0; j<phi.size(); j++) 
			    for (unsigned int k=0; k<_n_var; k++)
				Me(i + k*phi.size(), j + k*phi.size()) += JxW[qp]*phi[i][qp]*phi[j][qp];
			
		//scale if wished
		Me.scale(_scale_factor);
		 
		dof_map.constrain_element_matrix(Me, dof_indices, false);
	        //Me.print();

		matrix_A.add_matrix (Me, dof_indices);
		
	    }    
	      _console << "Assembling Solid System Matrix end" << std::endl;
    }

    if((_direction == TO_MULTIAPP)  && (_to_problems[0]->timeStep()== 1 || _first_recover)){
	   

            std::cout<<"timeToMultiapp=> "<< _to_problems[0]->timeStep()<<std::endl;

	    //if( _to_problems[0]->timeStep()==0 || _app.isRecovering()){

             _console << "Assembling Fluid System Matrix Begin" << std::endl;	


	    // Get a constant reference to the mesh object.
             const MeshBase & mesh = es.get_mesh();

            // The dimension that we are running.
            const unsigned int dim = mesh.mesh_dimension();

	    libMesh::SparseMatrix<Number> & _M_fluid = _linearImplicitSystemT->get_matrix("FluidMassMatrix");

	    // The element mass matrix.
            DenseMatrix<Number> Me_f;

            // A  Gauss quadrature rule for numerical integration.
            QGauss qrule (dim, fe_type.default_quadrature_order());
 
            // Tell the finite element object to use our quadrature rule.
            fe->attach_quadrature_rule (&qrule);

            // The element Jacobian * quadrature weight at each integration point.
            const std::vector<Real> & JxW = fe->get_JxW();

            // The element shape functions evaluated at the quadrature points.
            const std::vector<std::vector<Real> > & phi = fe->get_phi();

            const DofMap & dof_map = slave_system.get_dof_map();

            std::vector<dof_id_type> dof_indices;

            MeshBase::const_element_iterator       el_1     = mesh.active_local_elements_begin();
            const MeshBase::const_element_iterator end_el_1 = mesh.active_local_elements_end();

	    _M_fluid.zero(); /// _M_fluid

            for ( ; el_1 != end_el_1; ++el_1)
            {
                const Elem * elem = *el_1;
                Elem * ele = *el_1;

                dof_map.dof_indices (elem, dof_indices);

                fe->reinit (elem);

                Me_f.resize (dof_indices.size(), dof_indices.size());

                for (unsigned int qp=0; qp<qrule.n_points(); qp++)
                    for (unsigned int i=0; i<phi.size(); i++)
                        for (unsigned int j=0; j<phi.size(); j++)
                            for (unsigned int k=0; k<_n_var; k++)
                                Me_f(i + k*phi.size(), j + k*phi.size()) += JxW[qp]*phi[i][qp]*phi[j][qp];

                //scale if wished
                Me_f.scale(_scale_factor);

                dof_map.constrain_element_matrix(Me, dof_indices, false);
                //Me.print();

                
		 _M_fluid.add_matrix (Me_f, dof_indices);

		//_console << "Assembling Fluid System Matrix end" << std::endl;
		

            }

	       _console << "Assembling Fluid System Matrix end" << std::endl;

	        _M_fluid.close();

	       //_M_fluid.print_matlab("Ciao_ci_potrei_essere.m");
	// }  
    }

  if(_direction == TO_MULTIAPP){

	    _console << "Copy Fluid System Matrix into Slave System" << std::endl;

	    matrix_A.zero();
	    
	    //_M_fluid.print_matlab("Ciao_ci_sono.m");

	    //matrix_A = _M_fluid;

            libMesh::SparseMatrix<Number> & _M_fluid = _linearImplicitSystemT->get_matrix("FluidMassMatrix");

            //_M_fluid.print_matlab("Ciao_ci_sono.m");
	    
	    matrix_A.attach_dof_map(dof_map);
	    
	    matrix_A.close();
	    
	    matrix_A.add(1,_M_fluid);

	    //matrix_A.close();

	    //_M_fluid.print_matlab("Ciao_ci_sono.m");
	    //
	    if (_app.isRecovering()){ 
		    _first_recover = false;
	    }

    }

   // const_cast<StoreTransferOperators&>(_operator_storage).setFluidMass() = std::make_shared<SparseMatrix<Number>>(*slave_system.matrix);

    Vec tmp;
    
    VecT from_solution, tmp2;
   
    // This solution vector should technically contain the entire solution of the
    // system with all variables. 
    get_solution(_from_problems[0], _from_var_name_1, tmp);

    utopia::convert(tmp, tmp2);

    //utopia::disp(from_solution);
    // scale if wished 

   // if(_reverse)
   //{
         from_solution = tmp2;
   // }
   
      if (_impactEFSI && _direction == TO_MULTIAPP){

         _P_efsi = const_cast<StoreTransferOperators&>(_operator_storage).getPermuteOperatorDisp();

         from_solution = (*_P_efsi) * tmp2;
       }

   /* else{
        
      from_solution= tmp2;
    
   }*/
    
    VecT Bv =  (*_B) * from_solution * _scale_factor;

    VecT Bv_rev = Bv;
 

    
    if (_reverse==true)   {

       std::cout<<"Breverse:We need to permute"<<std::endl;	    
       
       _P = const_cast<StoreTransferOperators&>(_operator_storage).getPermuteOperator();
       
       Bv_rev  = (*_P) * Bv;

       //utopia::disp(Bv_rev);
    }




    //std::cout << "Bv size: " << Bv.size().get(0) << std::endl;   

    PetscInt local_begin, local_end, ncols;
    
    const PetscInt * cols, one = 1;
    
    const PetscScalar * vals;
    
    MooseVariable & to_var_1 = _to_problems[0]->getStandardVariable(0, _to_var_name_1);
   
    // this should be the same for all variables 
    System & to_sys = to_var_1.sys().system();
   
    // get iterators 
    MeshBase::const_node_iterator it = mesh.local_nodes_begin();
    const MeshBase::const_node_iterator end_it = mesh.local_nodes_end();
    
    // --- Set Up RHS ---
    // We said we'd manually take care of zeroing the RHS when transferring to IMPACT
    // (not the matrix however!)
    
    _console << "Assembling RHS Begin" << std::endl;
    
    slave_system.rhs->zero();

    for ( ; it != end_it; ++it){
        unsigned int slave_var_num_1 = slave_system.variable_number("var1");
        
        const Node * node = *it;
        
        for (unsigned int comp = 0;comp < node->n_comp(slave_system.number(), slave_var_num_1); comp++)
        {
            
            // Set RHS for first variable
            const dof_id_type to_index = node->dof_number(slave_system.number(), slave_var_num_1, comp);
            const dof_id_type from_index = node->dof_number(to_sys.number(), to_var_1.number(), comp);
            PetscScalar rhs;
            PetscInt index = from_index;
	    VecGetValues(utopia::raw_type(Bv_rev), one, &index, &rhs);
            slave_system.rhs->set(to_index, rhs);
        }
        if (_n_var > 1) {
           unsigned int slave_var_num_2 = slave_system.variable_number("var2");
           MooseVariable & to_var_2 = _to_problems[0]->getStandardVariable(0, _to_var_name_2);
           for (unsigned int comp = 0;comp < node->n_comp(slave_system.number(), slave_var_num_2); comp++)
           {
               // Set RHS for second variable
               const dof_id_type to_index = node->dof_number(slave_system.number(), slave_var_num_2, comp);
               const dof_id_type from_index = node->dof_number(to_sys.number(), to_var_2.number(), comp);
               PetscScalar rhs;
               PetscInt index = from_index;
               VecGetValues(utopia::raw_type(Bv_rev), one, &index, &rhs);
               slave_system.rhs->set(to_index, rhs);
           }
        }
        if (_n_var > 2) {
           unsigned int slave_var_num_3 = slave_system.variable_number("var3");
           MooseVariable & to_var_3 = _to_problems[0]->getStandardVariable(0, _to_var_name_3);
           for (unsigned int comp = 0;comp < node->n_comp(slave_system.number(), slave_var_num_3); comp++)
           {
               // Set RHS for third variable
               const dof_id_type to_index = node->dof_number(slave_system.number(), slave_var_num_3, comp);
               const dof_id_type from_index = node->dof_number(to_sys.number(), to_var_3.number(), comp);
               PetscScalar rhs;
               PetscInt index = from_index;
               VecGetValues(utopia::raw_type(Bv_rev), one, &index, &rhs);
               slave_system.rhs->set(to_index, rhs);
           }
        }
    }

    // this calls internal assembly routines
    matrix_A.close();  
    slave_system.rhs->close();
    
    _console << "Assembling System Matrix end" << std::endl;
    
    /*if (_direction == TO_MULTIAPP) {
       matrix_A.print_matlab("NewmatrixA_force.m");
      slave_system.rhs->print_matlab("Newrhs_force.m");
    } 
    else {
       matrix_A.print_matlab("NewmatrixA_vel.m");
       slave_system.rhs->print_matlab("Newrhs_vel.m");
    }*/

//_console << "Finished assembling the EquationSystem." << std::flush<< std::endl;
_console << "Matrix infinity norm: " << matrix_A.linfty_norm() << std::endl;
_console << "Matrix L1 norm: " << matrix_A.l1_norm() << std::endl;
_console << "Matrix size: " << matrix_A.m() << "x" << matrix_A.n() << std::endl;
_console << "RHS size: " << slave_system.rhs->size() << std::endl;
_console << "Finished assembling the EquationSystem." << std::flush<< std::endl;

}


void
IBTransfer::projectSolution()
{
    _console << "Updating L2-projection solution." << std::endl;
    
    FEProblemBase & to_problem = *_to_problems[0];
    EquationSystems & proj_es = to_problem.es();
    LinearImplicitSystem & linear_system = *_proj_sys;
    
    //linear_system.rhs->print();

    proj_es.parameters.set<IBTransfer *>("transfer") = this;
   

    //retrieve dofs from the elements and check whether their dofs participate in the
    //transfer. If they do at that element to the subdomain
   if (_direction == TO_MULTIAPP)
    {
      System & to_sys = (*_to_problems[0]).getStandardVariable(0,_to_var_name_1).sys().system();
      // Get a constant reference to the mesh object.
      const MeshBase & mesh = proj_es.get_mesh();
      
      // dofs are checked relative to this system
      const DofMap & dof_map = to_sys.get_dof_map();

      MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.local_elements_end();

      for ( ; el != end_el; ++el)
      {
        Elem * elem = *el;

        std::vector<dof_id_type> dof_indices;
        dof_map.dof_indices (elem, dof_indices);
      
        //elem->subdomain_id() = elem->invalid_subdomain_id; 
        elem->subdomain_id() = 0; //invalid_subdomain_id;; 
        for (std::vector<dof_id_type>::iterator it = dof_indices.begin() ; it != dof_indices.end() ; ++it)
        {
          if (_participating_slave_dofs.count(*it)) // do if they are participating
          {
            elem->subdomain_id() = 1; 
          } 
        }
      }
      //create the subdomain and restrict the solve to it
      std::set<subdomain_id_type> id_list;
      id_list.insert(1);
      SystemSubsetBySubdomain::SubdomainSelectionByList selection(id_list);
      SystemSubsetBySubdomain subset(linear_system, selection);
      linear_system.restrict_solve_to(&subset, SUBSET_ZERO);
    
      //Real tol = proj_es.parameters.get<Real>("linear solver tolerance");
      //proj_es.parameters.set<Real>("linear solver tolerance") = 1e-20;
      
      Real tol = linear_system.get_equation_systems().parameters.get<Real>("linear solver tolerance");
      //default here is 1e-5
      linear_system.get_equation_systems().parameters.set<Real>("linear solver tolerance") = _sys_tol;
      LinearSolver<Number> * linear_solver = linear_system.get_linear_solver();
      //_console << "old solver type: " << linear_solver->solver_type() << std::endl;
      linear_solver->set_solver_type(CG);
      //_console << "preconditioner type: " << linear_solver->preconditioner_type() << std::endl;

      //_console << "Solving our EquationSystem." << std::flush << std::endl;
      //_console << "solving with Solver type: " << linear_solver->solver_type() << std::endl;
      linear_system.solve();
      //_console << "Finished solving our EquationSystem." << std::flush << std::endl;
      // this just gives me 1, why?
      //_console << "Linear iterations used: "<< linear_system.n_linear_iterations() << std::endl;
      // this segfaults, why?
      //_console << "Solver convergence reason: " << linear_solver->get_converged_reason() << std::endl;
      // this just gives me 0, why?
      //_console << "Final linear residual: " << linear_system.final_linear_residual() << std::endl;
      //check solution for NaNs
      //_console << "sum of solution vector: " << linear_system.solution->sum() << std::endl;
 
      //proj_es.parameters.set<Real>("linear solver tolerance") = tol;
      // revert back to default
      linear_system.get_equation_systems().parameters.set<Real>("linear solver tolerance") = tol;
      linear_solver->set_solver_type(CG);
      
      set_solution(to_problem, proj_es);
      //_console << "Finished updating solution." << std::flush << std::endl;
    }
    else
    {
      Real tol = linear_system.get_equation_systems().parameters.get<Real>("linear solver tolerance");
      //default here is 1e-5
      linear_system.get_equation_systems().parameters.set<Real>("linear solver tolerance") = _sys_tol;
    
      LinearSolver<Number> * linear_solver = linear_system.get_linear_solver();
    
      linear_solver->set_solver_type(CG);

      _console << "Solving our EquationSystem." << std::flush << std::endl;
    
      linear_system.solve();
    
      _console << "Finished solving our EquationSystem." << std::flush << std::endl;
    
      // this just gives me 1, why?
      //_console << "Linear iterations used (vel): "<< linear_system.n_linear_iterations() << std::endl;
      // this segfaults, why?
      //_console << "Solver convergence reason: " << linear_solver->get_converged_reason() << std::endl;
      // this just gives me 0, why?
      //_console << "Final linear residual (vel): " << linear_system.final_linear_residual() << std::endl;
     
      //revert back to default 
      linear_system.get_equation_systems().parameters.set<Real>("linear solver tolerance") = tol;

      //linear_solver->set_solver_type(CG);
      _console << "Starting updating solution." << std::endl;      

      set_solution(to_problem, proj_es);
    
      _console << "Finished updating solution." << std::endl;

    }

}



void
IBTransfer::set_solution(FEProblemBase &problem,
                                                EquationSystems &proj_es)
{

    _console << "Set solution." << std::endl;
    // copy projected solution into target es
    MeshBase & mesh = proj_es.get_mesh();
    MooseVariable & var1 = problem.getStandardVariable(0, _to_var_name_1);
   
    // solution of the original system 
    System & sys = var1.sys().system();
    NumericVector<Number> * solution = sys.solution.get();
    
    LinearImplicitSystem & linear_system = *_proj_sys;

    { // loop through our local elements and set the solution from the projection
        
        MeshBase::const_node_iterator it = mesh.local_nodes_begin();
        const MeshBase::const_node_iterator end_it = mesh.local_nodes_end();
        for (; it != end_it; ++it)
        {
            // set solution for variable 1
            const Node * node = *it;
            for (unsigned int comp = 0;
                 comp < node->n_comp(sys.number(), var1.number()); comp++)
            {
                const dof_id_type proj_index =
                node->dof_number(linear_system.number(), _proj_var_num1, comp);
                const dof_id_type index =
                node->dof_number(sys.number(), var1.number(), comp);
		//std::cout<<"value"<<(*linear_system.solution)(proj_index)<<std::endl;
                solution->set(index, (*linear_system.solution)(proj_index));
            }
            // set solution for variable 2
            if (_n_var > 1) {
               MooseVariable & var2 = problem.getStandardVariable(0, _to_var_name_2);
               for (unsigned int comp = 0;
                    comp < node->n_comp(sys.number(), var2.number()); comp++)
               {
                   const dof_id_type proj_index =
                   node->dof_number(linear_system.number(), _proj_var_num2, comp);
                   const dof_id_type index =
                   node->dof_number(sys.number(), var2.number(), comp);
                   solution->set(index, (*linear_system.solution)(proj_index));
               }
            }
            // set solution for variable 3
            if (_n_var > 2) {
               MooseVariable & var3 = problem.getStandardVariable(0, _to_var_name_3);
               for (unsigned int comp = 0;
                    comp < node->n_comp(sys.number(), var3.number()); comp++)
               {
                   const dof_id_type proj_index =
                   node->dof_number(linear_system.number(), _proj_var_num3, comp);
                   const dof_id_type index =
                   node->dof_number(sys.number(), var3.number(), comp);
                   solution->set(index, (*linear_system.solution)(proj_index));
               }
            }
        }
    }
    
    
    solution->close();
//    solution->print_matlab("ciao.m");
    sys.update();
}


void
IBTransfer::Permutation(FEProblemBase *problem, std::shared_ptr<SparseMatT> P_matrix, unsigned int number){


    _console<<"Permutation FSI::begin"<<std::endl;

    MeshBase & mesh = problem->mesh();

    MeshBase::const_element_iterator it = mesh.local_elements_begin();

    const MeshBase::const_element_iterator end_it = mesh.local_elements_end();

    std::vector<dof_id_type> temp_main, temp_aux;

    MooseVariable & aux_var = problem->getStandardVariable(0, _to_var_name_1);

    MooseVariable & main_var = problem->getStandardVariable(0, _to_var_name_r);

    
    DofMap & main_dof = main_var.sys().system().get_dof_map(); //solid_problem


    DofMap & auxi_dof = aux_var.sys().system().get_dof_map(); //solid_problem
    
    //DofMap & auxi_dof = problem->es().get_system(problem_name).get_dof_map(); //auxi_var.sys().system().get_dof_map();  //fluid_problem

    (*P_matrix)=utopia::local_sparse(auxi_dof.n_local_dofs(),main_dof.n_local_dofs(),1.0);

   {
     utopia::Write<SparseMatT> d(*(P_matrix));
     
     auto RowRange=utopia::row_range(*P_matrix);
     
     auto ColRange=utopia::col_range(*P_matrix);

     for ( ; it != end_it; ++it)
     {
         const Elem * elem = *it;

         main_dof.dof_indices(elem, temp_main, main_var.number());

         auxi_dof.dof_indices(elem, temp_aux, number);

         for(int i=0; i<temp_aux.size(); i++){

            const long dof_I = temp_aux[i];

            const long dof_J = temp_main[i];

            if (RowRange.inside(dof_I) && ColRange.inside(dof_J)){

                for(utopia::SizeType d = 0; d < _n_var; ++d){
     
	  		(*P_matrix).set(dof_I+d,dof_J+d,1.0);
                }
            }
         }
      }
   }

   //utopia:: disp(*P_matrix);
   _console<<"Permutation FSI::out"<<std::endl;
}


void
IBTransfer::PermutationEFSI(FEProblemBase *problem, std::shared_ptr<SparseMatT> P_matrix, unsigned int number){


    _console<<"Permutation EFSI::begin"<<std::endl;

    MeshBase & mesh = problem->mesh();

    MeshBase::const_element_iterator it = mesh.local_elements_begin();

    const MeshBase::const_element_iterator end_it = mesh.local_elements_end();

    std::vector<dof_id_type> temp_d, temp_f;

    MooseVariable & f_var = problem->getStandardVariable(0, "res_elastic_scaled_sub_x");

    MooseVariable & d_var = problem->getStandardVariable(0, "disp_x");

    DofMap & d_dof = f_var.sys().system().get_dof_map(); //solid_problem


    DofMap & f_dof = f_var.sys().system().get_dof_map(); //solid_problem
   
    //DofMap & auxi_dof = problem->es().get_system(problem_name).get_dof_map(); //auxi_var.sys().system().get_dof_map();  //fluid_problem

    (*P_matrix)=utopia::local_sparse(d_dof.n_local_dofs(),f_dof.n_local_dofs(),1.0);

   {
     utopia::Write<SparseMatT> d(*(P_matrix));
     auto RowRange=utopia::row_range(*P_matrix);
     auto ColRange=utopia::col_range(*P_matrix);

     for ( ; it != end_it; ++it)
     {
         const Elem * elem = *it;

         d_dof.dof_indices(elem, temp_d, d_var.number());

         f_dof.dof_indices(elem, temp_f, f_var.number());

         for(int i=0; i<temp_f.size(); i++){

            const long dof_I = temp_d[i];

            const long dof_J = temp_f[i];

            if (RowRange.inside(dof_I) && ColRange.inside(dof_J)){

                for(utopia::SizeType d = 0; d < _n_var; ++d){
                   (*P_matrix).set(dof_I+d,dof_J+d,1.0);
                }
            }
        }
     }
   }

   //utopia:: disp(*P_matrix);
   _console<<"Permutation EFSI::out"<<std::endl;
}
