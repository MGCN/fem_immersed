// /****************************************************************/
// /*               DO NOT MODIFY THIS HEADER                      */
// /* MOOSE - Multiphysics Object Oriented Simulation Environment  */
// /*                                                              */
// /*           (c) 2010 Battelle Energy Alliance, LLC             */
// /*                   ALL RIGHTS RESERVED                        */
// /*                                                              */
// /*          Prepared by Battelle Energy Alliance, LLC           */
// /*            Under Contract No. DE-AC07-05ID14517              */
// /*            With the U. S. Department of Energy               */
// /*                                                              */
// /*            See COPYRIGHT for full restrictions               */
// /****************************************************************/

// /****************************************************************/
// /*               DO NOT MODIFY THIS HEADER                      */
// /*    Immersed_Boundary- ICS Mechanical simulation framework    */
// /*                Prepared by Maria Nestola,                    */
// /*                  ICS, USI, 6900 Lugano                       */
// /*                                                              */
// /****************************************************************/


// #include "TransferVector.h"
// #include "FEProblem.h"
// #include "AddVariableAction.h"
// #include "NonlinearSystemBase.h"
// #include "MooseMesh.h"
// #include "DisplacedProblem.h"
// #include <opencl_adapter.hpp>
// #include "libmesh/linear_partitioner.h"
// #include "libmesh/petsc_matrix.h"
// #include "utopia_assemble_volume_transfer.hpp"
// #include <petscksp.h>
// #include <petscmat.h>
// #include "MooseVariable.h"
// #include <utopia.hpp>
// #include "utopia_TransferAssembler.hpp"
// #include "utopia_MeshTransferOperator.hpp"

// #include <cmath>
// #include <queue>

// typedef utopia::DSMatrixd SparseMatT;
// typedef utopia::DVectord VecT;




// registerMooseObject("ImmersedBoundaryApp", TransferVector);

// template <>
// InputParameters
// validParams<TransferVector>()
// {
//   InputParameters params = validParams<MultiAppTransfer>();
//   params.addRequiredParam<AuxVariableName>("variable",
//                                            "The auxiliary variable to store the transferred values in.");
//   params.addRequiredParam<VariableName>("source_variable",
//                                         "The variable to transfer from.");
//   params.addParam<bool>("displaced_source_mesh", false,
//                         "Use the source mesh after displacement");
//   params.addParam<bool>("displaced_target_mesh", false,
//                         "Use the target mesh after displacement");
//   params.addParam<bool>("impact",false,
//                           "impact.");
//   params.addRequiredParam<int>("num_of_variables",false,
//                           "num_of_variables");
//   params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
//   params.addParam<std::string>("operator_type","opeartor_type" /*INTERPOLATION| L2_PROJECTION| PSEUDO_L2_PROJECTION | APPROX_L2_PROJECTION*/);
//   return params;

// }

// TransferVector::TransferVector(const InputParameters & parameters)
//   :  MultiAppTransfer(parameters),
//      UserObjectInterface(this),
//     _to_var_name(getParam<AuxVariableName>("variable")),
//     _from_var_name(getParam<VariableName>("source_variable")),
//     _impact(getParam<bool>("impact")),
//     _operator_storage(getUserObject<StoreTransferOperators>("operator_userobject")),
//     _operator_type(getParam<std::string>("operator_type")),
//     _n_var(getParam<int>("num_of_variables"))

// {
// }

// void
// TransferVector::initialSetup()
// {
//     _console << "Initial Setup of Transfer System!" << std::endl;
//     getAppInfo();

//     if (_to_problems.size() > 1 || _from_problems.size() > 1) {
//         mooseError("NOT IMPLEMENTED: Right now we can only handle one master and one slave");
//     }
    
//     auto _proj_sys_to = &_to_problems[0]->es().add_system<LinearImplicitSystem>("proj-sys-slave");
    
//     auto _proj_sys_from = &_from_problems[0]->es().add_system<LinearImplicitSystem>("proj-sys-master");
    
//     _var_1_from  = _proj_sys_from->add_variable("var_1",FIRST);
    
//     if (_n_var > 1){
//        _var_2_from =  _proj_sys_from->add_variable("var_2",FIRST);
//     }
    
//     if (_n_var > 2){
//        _var_3_from = _proj_sys_from->add_variable("var_3",FIRST);
//     }
    
//     _from_problems[0]->es().reinit();
    
//     _var_1_to  =  _proj_sys_to->add_variable("var_1",FIRST);
    
//     if (_n_var > 1){
//       _var_2_to =  _proj_sys_to->add_variable("var_2",FIRST);
//     }
    
//     if (_n_var > 2){
//       _var_3_to =  _proj_sys_to->add_variable("var_3",FIRST);
//     }
//     _to_problems[0]->es().reinit();
    

// }

// void
// TransferVector::execute()
// {
    
//     _console << "L2ProjectionLibMeshTransfer transfer " << name() << std::endl;
    
//     buildTransferOperators();

// }


// void
// TransferVector::buildTransferOperators()
// {
//     using namespace utopia;

//     _console << " BuildTransferOperators for variable " << name() << std::endl;

    
//     FEProblemBase & from_problem = *_from_problems[0];
   
//     MeshBase *from_mesh = &from_problem.mesh().getMesh();

//     // displaced meshes are special
//     if (_displaced_source_mesh && from_problem.getDisplacedProblem()){from_mesh = &from_problem.getDisplacedProblem()->mesh().getMesh();}

//     FEProblemBase & to_problem = *_to_problems[0];
    
//     MeshBase * to_mesh = &to_problem.mesh().getMesh();

//     if (_displaced_target_mesh && to_problem.getDisplacedProblem()){to_mesh = &to_problem.getDisplacedProblem()->mesh().getMesh();}

//     moonolith::Communicator expressComm(_communicator.get());

//     MooseVariable & _from_var = _from_problems[0]->getStandardVariable(0, _from_var_name);

//     MooseVariable & _to_var = _to_problems[0]->getStandardVariable(0, _to_var_name);

//     DofMap &master_dof = _from_var.sys().system().get_dof_map(); //_from_problems[0]->es().get_system("nl0").get_dof_map();//

//     DofMap &slave_dof  = _to_var.sys().system().get_dof_map(); // _to_problems[0]->es().get_system("aux0").get_dof_map();//

//     unsigned int n_var = 1;

//     TransferOptions opts;
//     opts.from_var_num = _from_var.number();
//     opts.to_var_num   = _to_var.number();
//     opts.n_var        = _n_var;
//     opts.tags         = {};
    
//     _P = std::make_shared<SparseMatT>();
    
//     _P_reverse    = std::make_shared<SparseMatT>();
    
//     if(! const_cast<StoreTransferOperators&>(_operator_storage).getTransferOperator()){
        
//         Permutation(_from_problems[0], _P, _from_var_name, "proj-sys-master", _var_1_from); //from solid to fluid
        
//     }
    
//     if(! const_cast<StoreTransferOperators&>(_operator_storage).getTransferOperatorReverse()){
        
//         Permutation(_from_problems[0], _P_reverse, _from_var_name, "proj-sys-slave", _var_1_from); //from fluid to solid
        
//     }

 
    
//     {
//         auto cs = std::make_shared<utopia::MeshTransferOperator>(make_ref(*from_mesh),make_ref(master_dof),make_ref(*to_mesh),make_ref(slave_dof),opts);
        
//         cs->initialize(_operator_type);
        
//        const_cast<StoreTransferOperators&>(_operator_storage).getVoidPointer() = cs;
//     }

//     if (_direction == TO_MULTIAPP){
//         apply_operator();

//     }

//     if (_direction == FROM_MULTIAPP){
//         if(_impact) {
//             apply_operator();
//         }
//         else{
//             apply_transpose_operator();
//         }
//     }
// }

// void
// TransferVector::get_solution_vec(FEProblemBase * problem, std::string var_name,
//                                                  Vec &vector)
// {

//     _console << "Get Solution Vector" << std::flush << std::endl;

//     MooseVariable & var = problem->getStandardVariable(0, var_name);
//     System & sys = var.sys().system();

//     vector = static_cast<PetscVector<Number> *>(sys.current_local_solution.get())->vec();

// }

// void
// TransferVector::apply_operator()
// {

//     Vec tmp;
    
//     VecT from_sol, to_sol;

//     get_solution_vec(_from_problems[0], _from_var_name, tmp);
    
//     utopia::convert(tmp, from_sol);

//     std::shared_ptr<utopia::MeshTransferOperator> cs;

//     cs = std::static_pointer_cast<utopia::MeshTransferOperator>(const_cast<StoreTransferOperators&>(_operator_storage).getVoidPointer());

//     //cs->apply(from_sol,to_sol);
    
//     _P = const_cast<StoreTransferOperators&>(_operator_storage).getTransferOperator();
    
//     auto from_sol_per=(*_P)*from_sol;
    
//     cs->apply(from_sol_per,to_sol);
   
//     copy_sol(to_sol);

// }

// void
// TransferVector::apply_transpose_operator()
// {

//     Vec tmp;
//     get_solution_vec(_from_problems[0], _from_var_name, tmp);

//     std::cout<<" _from_var_name " << _from_var_name <<std::endl;

//     VecT from_sol, to_sol;

//     utopia::convert(tmp, from_sol);

//     std::shared_ptr<utopia::MeshTransferOperator> cs;

//     cs = std::static_pointer_cast<utopia::MeshTransferOperator>(const_cast<StoreTransferOperators&>(_operator_storage).getVoidPointer());

//     cs->apply_transpose(from_sol,to_sol);


//     _console << "Set PSEUDO-L2-projection Solution Vector" << std::flush << std::endl;

//     FEProblemBase & to_problem = *_to_problems[0];

//     MeshBase & to_mesh = _to_problems[0]->mesh();

//     MooseVariable & to_var = to_problem.getStandardVariable(0, _to_var_name);

//     //std::cout<<"to_var_name"<<_to_var_name<<std::endl;
//     //disp(sol);

//     System & to_sys = to_var.sys().system();
//     NumericVector<Number> * to_solution = to_sys.solution.get();
//     PetscInt       rstart,rend;
//     PetscScalar    tmp_sol;
//     VecGetOwnershipRange(utopia::raw_type(to_sol),&rstart,&rend);
//     utopia::Read<VecT> w_d(to_sol);


//     {
//         MeshBase::const_node_iterator it = to_mesh.local_nodes_begin();
//         const MeshBase::const_node_iterator end_it = to_mesh.local_nodes_end();
//         for ( ; it != end_it; ++it)
//         {
//             const Node * node = *it;

//             for (unsigned int comp = 0;comp < node->n_comp(to_sys.number(), to_var.number()); comp++)
//             {
//                 const dof_id_type to_index = node->dof_number(to_sys.number(), to_var.number(), comp);

//                 to_solution->set(to_index,  to_sol.get(to_index));

//             }
//         }
//     }

//     {
//         MeshBase::const_element_iterator it = to_mesh.active_local_elements_begin();
//         const MeshBase::const_element_iterator end_it = to_mesh.active_local_elements_end();
//         for ( ; it != end_it; ++it)
//         {
//             const Elem * elem = *it;
//             for (unsigned int comp = 0; comp < elem->n_comp(to_sys.number(), to_var.number()); comp++)
//             {
//                 const dof_id_type to_index = elem->dof_number(to_sys.number(), to_var.number(), comp);
//                 to_solution->set(to_index, to_sol.get(to_index));
//             }
//         }
//     }


//     to_solution->close();
//     to_sys.update();

// }

// void
// TransferVector::Permutation(FEProblemBase *problem, std::shared_ptr<SparseMatT> P_matrix, std::string main_var_name, std::string problem_name, unsigned int number){
    
    
//     std::cout<<"Permutation::begin"<<std::endl;
    
//     MeshBase & mesh = problem->mesh();
    
//     MeshBase::const_element_iterator it = mesh.local_elements_begin();
    
//     const MeshBase::const_element_iterator end_it = mesh.local_elements_end();
    
//     std::vector<dof_id_type> temp_main, temp_aux;
    
//     MooseVariable & main_var = problem->getStandardVariable(0, main_var_name);
    
//     DofMap & main_dof = main_var.sys().system().get_dof_map(); //solid_problem
    
//     DofMap & auxi_dof = problem->es().get_system(problem_name).get_dof_map(); //auxi_var.sys().system().get_dof_map();  //fluid_problem
    
//     //std::cout<<"Permutation Loop"<<std::endl;
    
//     (*P_matrix)=utopia::local_sparse(auxi_dof.n_local_dofs(),main_dof.n_local_dofs(),1.0);
    
//     {
//         utopia::Write<SparseMatT> d(*(P_matrix));
//         auto RowRange=utopia::row_range(*P_matrix);
//         auto ColRange=utopia::col_range(*P_matrix);
        
//         for ( ; it != end_it; ++it)
//         {
//             const Elem * elem = *it;
            
//             main_dof.dof_indices(elem, temp_main, main_var.number());
            
//             auxi_dof.dof_indices(elem, temp_aux, number);
            
//             for(int i=0; i<temp_aux.size(); i++){
                
//                 const long dof_I = temp_aux[i];
                
//                 const long dof_J = temp_main[i];
                
//                 if (RowRange.inside(dof_I) && ColRange.inside(dof_J)){
                    
//                     for(utopia::SizeType d = 0; d < _n_var; ++d){
//                         (*P_matrix).set(dof_I+d,dof_J+d,1.0);
//                     }
//                 }
//             }
//         }
//     }
//     std::cout<<"Permutation::end"<<std::endl;
// }


// void
// L2ProjectionLibMeshTransferBidirectional::copy_sol(VecT to_sol)
// {
//     _console << "Set Solution Vector" << std::flush << std::endl;
    
//     FEProblemBase & to_problem = *_to_problems[0];
    
//     MeshBase & to_mesh = _to_problems[0]->mesh();
    
//     std::cout<<"_to_var_name " << _to_var_name <<std::endl;
    
//     MooseVariable & to_var = to_problem.getStandardVariable(0, _to_var_name);
    
//     utopia::Read<VecT> w_d(to_sol);
    
//     System & to_sys = to_var.sys().system();
    
//     NumericVector<Number> * to_solution = to_sys.solution.get();
    
//     auto & proj_sys = to_problem.es().get_system("proj-sys-master");
    
//     {
//         MeshBase::const_node_iterator it = to_mesh.local_nodes_begin();
//         const MeshBase::const_node_iterator end_it = to_mesh.local_nodes_end();
//         for ( ; it != end_it; ++it)
//         {
//             const Node * node = *it;
            
//             for (unsigned int comp = 0;comp < node->n_comp(to_sys.number(), to_var.number()); comp++)
//             {
//                 const dof_id_type to_index = node->dof_number(to_sys.number(), to_var.number(), comp);
//                 const dof_id_type pr_index = node->dof_number(proj_sys.number(), _var_1_from , comp);
//                 to_solution->set(to_index,  to_sol.get(pr_index));
                
//             }
//         }
//     }
    
//     {
//         MeshBase::const_element_iterator it = to_mesh.active_local_elements_begin();
//         const MeshBase::const_element_iterator end_it = to_mesh.active_local_elements_end();
//         for ( ; it != end_it; ++it)
//         {
//             const Elem * elem = *it;
//             for (unsigned int comp = 0; comp < elem->n_comp(to_sys.number(), to_var.number()); comp++)
//             {
//                 const dof_id_type to_index = elem->dof_number(to_sys.number(), to_var.number(), comp);
//                 const dof_id_type pr_index = elem->dof_number(proj_sys.number(), _var_1_from , comp);
//                 to_solution->set(to_index, to_sol.get(pr_index));
//             }
//         }
//     }
    
//     to_solution->close();
//     to_sys.update();
    
    
//     if (_n_var > 1){
        
//         MooseVariable & slave_var_2 = to_problem.getStandardVariable(0, _to_var_name_2);
        
//         {
//             MeshBase::const_node_iterator it = to_mesh.local_nodes_begin();
//             const MeshBase::const_node_iterator end_it = to_mesh.local_nodes_end();
//             for ( ; it != end_it; ++it){
//                 const Node * node = *it;
//                 for (unsigned int comp = 0;comp < node->n_comp(to_sys.number(), slave_var_2.number()); comp++){
//                     const dof_id_type to_index = node->dof_number(to_sys.number(), slave_var_2.number(), comp);
//                     const dof_id_type pr_index = node->dof_number(proj_sys.number(), _var_2_from , comp);
//                     to_solution->set(to_index,  to_sol.get(pr_index));
//                 }
//             }
//         }
        
//         {
//             MeshBase::const_element_iterator it = to_mesh.active_local_elements_begin();
//             const MeshBase::const_element_iterator end_it = to_mesh.active_local_elements_end();
//             for ( ; it != end_it; ++it)
//             {
//                 const Elem * elem = *it;
//                 for (unsigned int comp = 0; comp < elem->n_comp(to_sys.number(), slave_var_2.number()); comp++)
//                 {
//                     const dof_id_type to_index = elem->dof_number(to_sys.number(), slave_var_2.number(), comp);
//                     const dof_id_type pr_index = elem->dof_number(proj_sys.number(), _var_2_from , comp);
//                     to_solution->set(to_index, to_sol.get(pr_index));
//                 }
//             }
//         }
//     }
    
    
//     to_solution->close();
//     to_sys.update();
    
    
//     if (_n_var > 2){
        
//         MooseVariable & slave_var_3 = to_problem.getStandardVariable(0, _to_var_name_3);
        
//         {
//             MeshBase::const_node_iterator it = to_mesh.local_nodes_begin();
//             const MeshBase::const_node_iterator end_it = to_mesh.local_nodes_end();
//             for ( ; it != end_it; ++it){
//                 const Node * node = *it;
//                 for (unsigned int comp = 0;comp < node->n_comp(to_sys.number(), slave_var_3.number()); comp++){
//                     const dof_id_type to_index = node->dof_number(to_sys.number(), slave_var_3.number(), comp);
//                     const dof_id_type pr_index = node->dof_number(proj_sys.number(), _var_3_from , comp);
//                     to_solution->set(to_index,  to_sol.get(pr_index));
//                 }
//             }
//         }
        
//         {
//             MeshBase::const_element_iterator it = to_mesh.active_local_elements_begin();
//             const MeshBase::const_element_iterator end_it = to_mesh.active_local_elements_end();
//             for ( ; it != end_it; ++it)
//             {
//                 const Elem * elem = *it;
//                 for (unsigned int comp = 0; comp < elem->n_comp(to_sys.number(), slave_var_3.number()); comp++)
//                 {
//                     const dof_id_type to_index = elem->dof_number(to_sys.number(), slave_var_3.number(), comp);
//                     const dof_id_type pr_index = elem->dof_number(proj_sys.number(), _var_1_from , comp);
//                     to_solution->set(to_index, to_sol.get(pr_index));
//                 }
//             }
//         }
//     }
    
//     to_solution->close();
//     to_sys.update();
    
    
// }
