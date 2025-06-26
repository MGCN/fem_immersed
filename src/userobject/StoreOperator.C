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
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#include "StoreOperator.h"
//#include "L2ProjectionLibMeshTransfer.h"
#include "Executioner.h"
#include "MultiApp.h"
//#include "L2ProjectionLibMeshTransfer.h"
#include "FEProblem.h"
#include "AddVariableAction.h"
#include "MooseMesh.h"
#include "DisplacedProblem.h"
#include <opencl_adapter.hpp>
#include "libmesh/linear_partitioner.h"
#include "libmesh/petsc_matrix.h"
#include <petscksp.h>
#include <petscmat.h>
#include "MooseVariable.h"
#include <utopia.hpp>


#include <cmath>
#include <queue>
#include <memory>

//utopia::USparseMatrix SparseMatT;

//using namespace utopia;

registerMooseObject("ImmersedBoundaryApp", StoreOperator);

template<>
InputParameters validParams<StoreOperator>()
{
    InputParameters params = validParams<GeneralUserObject>();
    params.addRequiredParam<MultiAppName>("multi_app", "The name of the MultiApp to use.");
    params.addRequiredParam<AuxVariableName>("solid_variable",
                                             "The auxiliary variable to store the transferred values in.");
    params.addRequiredParam<VariableName>("fluid_variable",
                                          "The variable to transfer from.");
    
    params.addParam<bool>("fixed_meshes", false,
                          "Set to true when the meshes are not changing (ie, no "
                          "movement or adaptivity).  This will cache some "
                          "information to speed up the transfer.");
    params.addParam<bool>("moved_source_meshes", false,
                          "Set to true when the meshes are not changing (ie, no "
                          "movement or adaptivity).  This will cache some "
                          "information to speed up the transfer.");
    
    params.addParam<bool>("moved_target_meshes", false,
                          "Set to true when the meshes are not changing (ie, no "
                          "movement or adaptivity).  This will cache some "
                          "information to speed up the transfer.");
    params.addParam<bool>("impact",false,
                          "impact");
    params.addParam<bool>("pseudoL2",false,
                          "Set to false when you want to pure L2-projection.");
    params.addParam<bool>("displaced_source_mesh", false,
                          "Use the source mesh after displacement");
    params.addParam<bool>("displaced_target_mesh", false,
                          "Use the target mesh after displacement");
    params.addRequiredParam<MultiAppName>("multi_app", "The MultiApp's name in your input file!");
    
    
    
    
    return params;
}

StoreOperator::StoreOperator(const InputParameters & parameters) :
GeneralUserObject(parameters),
_fe_problem(parameters.get< FEProblem *>("_fe_problem")),
_to_var_name(getParam<AuxVariableName>("solid_variable")),
_from_var_name(getParam<VariableName>("fluid_variable")),
_impact(getParam<bool>("impact")),
_biorthogonal(getParam<bool>("pseudoL2")),
_fixed_meshes(getParam<bool>("fixed_meshes")),
_moved_source_mesh(getParam<bool>("moved_source_meshes")),
_moved_target_mesh(getParam<bool>("moved_target_meshes")),
_multiapp_name(getParam<MultiAppName>("multi_app"))

{
    
    
    if (_fixed_meshes) {
        _cache_matrices = true;
    }
    else{
        _displaced_source_mesh = getParam<bool>("displaced_source_mesh");
        _displaced_target_mesh = getParam<bool>("displaced_target_mesh");
    }
    
    if (_biorthogonal) {
        _pseudoL2= true;
    }
    
    
    
    
}



void
StoreOperator::initialize()
{
    //    MultiApp &  _multi_app = * _fe_problem->getMultiApp(_multiapp_name);
    //
    //    FEProblemBase & to_problem = _multi_app.problemBase();
    //
    //    FEProblemBase & from_problem = _multi_app.appProblemBase(0);
    
    //    FEType fe_type = to_problem.getVariable(0, _to_var_name).feType();
    //
    //    LinearImplicitSystem & mass_matrix_sys =to_problem.es().add_system<LinearImplicitSystem>("Mass_matrix_system");
    //
    //    unsigned int _proj_var_num =  mass_matrix_sys.add_variable("var", fe_type);
    
    // Reinitialize EquationSystems since we added a system.
    //    to_problem.es().reinit();
}

StoreOperator::~StoreOperator()
{}

void
StoreOperator::execute()
{
   
    using namespace utopia;

    _B = std::make_shared<SparseMatT>();
    //_direction = MultiAppTransfer::directions();    
    // std::cout<<"direction = "<<MultiAppTransfer::direction()<<std::endl;
    
    MultiApp &  _multi_app = * _fe_problem->getMultiApp(_multiapp_name);
    
    //std::shared_ptr<MultiAppTransfer> multi_app_transfer = std::dynamic_pointer_cast<MultiAppTransfer>(transfer);
    
    // int a=MultiAppTransfer::direction();
    
    
    MeshBase  *_mesh_moved = & _subproblem.mesh().getMesh();

    FEProblemBase & to_problem = _multi_app.problemBase(); //SOLID
    
    FEProblemBase & from_problem = _multi_app.appProblemBase(0);   //FLUID
    
    MeshBase *from_mesh = &from_problem.mesh().getMesh();
    
    MeshBase *to_mesh = &to_problem.mesh().getMesh();
    
    if (_moved_source_mesh){
        from_mesh = &from_problem.getDisplacedProblem()->mesh().getMesh();
        std::cout<<"I am using the source displaced Mesh"<<std::endl;
    }
    
    
    if (_moved_target_mesh && to_problem.getDisplacedProblem()){
        to_mesh = &to_problem.getDisplacedProblem()->mesh().getMesh();
        std::cout<<"I am using the target displaced Mesh"<<std::endl;
    }
    
    
    PetscErrorCode ierr;
    
    
    //moonolith::Communicator expressComm(_communicator.get());
    
    MooseVariable & _from_var = from_problem.getStandardVariable(0,  _from_var_name);
    
    DofMap &master_dof = _from_var.sys().system().get_dof_map();
    
    MooseVariable & _to_var = to_problem.getStandardVariable(0, _to_var_name);
    
    DofMap &slave_dof  = _to_var.sys().system().get_dof_map();
    
    
}

void
StoreOperator::assemble_mass_matrix(EquationSystems & es,const std::string & system_name, SparseMatT &D)
{
    //   // using namespace utopia;
    //
    //    // Get a constant reference to the mesh object.
    //    const MeshBase & mesh = es.get_mesh();
    //
    //    // The dimension that we are running.
    //    const unsigned int dim = mesh.mesh_dimension();
    //
    //    // Get a reference to our system.
    //    LinearImplicitSystem & slave_system = es.get_system<LinearImplicitSystem>("Mass_matrix_system");
    //
    //    // Get a constant reference to the Finite Element type
    //    // for the first (and only) variable in the system.
    //    FEType fe_type = slave_system.get_dof_map().variable_type(0);
    //
    //    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    //
    //    // A  Gauss quadrature rule for numerical integration.
    //    QGauss qrule (dim, fe_type.default_quadrature_order());
    //
    //    // Tell the finite element object to use our quadrature rule.
    //    fe->attach_quadrature_rule (&qrule);
    //
    //    // The element Jacobian * quadrature weight at each integration point.
    //    const std::vector<Real> & JxW = fe->get_JxW();
    //
    //    // The element shape functions evaluated at the quadrature points.
    //    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    //
    //
    //    const DofMap & dof_map = slave_system.get_dof_map();
    //
    //
    //    std::vector<dof_id_type> dof_indices;
    //
    //    std::vector<dof_id_type> dof_indices_size;
    //
    //
    //    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    //    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    //
    //    const Elem * elem = *el;
    //
    //    dof_map.dof_indices( elem, dof_indices_size);
    //
    //    express::MapSparseMatrix<double> mat_buffer(dof_indices_size.size(), dof_indices_size.size());
    //
    //    for ( ; el != end_el; ++el)
    //    {
    //
    //        const Elem * elem = *el;
    //
    //
    //        dof_map.dof_indices(elem, dof_indices);
    //
    //
    //        fe->reinit (elem);
    //
    //        libMesh::DenseMatrix<libMesh::Real> elemmat;
    //
    //        elemmat.resize(phi.size(), phi.size());
    //
    //        elemmat.zero();
    //
    //        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    //            for (unsigned int i=0; i<phi.size(); i++)
    //                for (unsigned int j=0; j<phi.size(); j++)
    //                     elemmat(i,j)+=JxW[qp]*phi[i][qp]*phi[j][qp];
    //
    //
    //        for (int i=0; i< elemmat.m(); ++i)
    //            for (int j=0; j< elemmat.n(); ++j)
    //                if (i!=j)
    //                {
    //                     elemmat(i,i)+=elemmat(i,j);
    //                     elemmat(i,j)=0.0;
    //
    //                }
    //
    //
    //
    //        for(int i = 0; i < dof_indices.size(); ++i) {
    //
    //            const long dof_I =dof_indices[i];
    //
    //            for(int j = 0; j < dof_indices.size(); ++j) {
    //
    //                const long dof_J =dof_indices[j];
    //
    //                mat_buffer.add(dof_I, dof_J, elemmat(i, j));
    //            }
    //        }
    //
    //
    //   //     matrix_A.add_matrix (Me, dof_indices);
    //
    //    }
    //
    //
    //    express::Communicator comm(_communicator.get());
    //
    //    express::Array<express::SizeType>  ownershipRanges(comm.size()+1);
    //
    //    ownershipRanges.allSet(0);
    //
    //    ownershipRanges[comm.rank()+1] += static_cast<unsigned int>(dof_map.n_local_dofs());
    //
    //    comm.allReduce(&ownershipRanges[0], ownershipRanges.size(), express::MPISum());
    //
    //    std::partial_sum(ownershipRanges.begin(), ownershipRanges.end(), ownershipRanges.begin());
    //
    //
    //
    //    if(comm.isRoot()) {
    //        std::cout <<ownershipRanges << std::endl;
    //
    //    }
    //
    //    express::Redistribute< express::MapSparseMatrix<double> > redist(comm.getMPIComm());
    //
    //    redist.apply(ownershipRanges, mat_buffer, express::AddAssign<double>());
    //
    //    express::RootDescribe("petsc assembly begin", comm, std::cout);
    //
    //    const utopia::SizeType local_range_raw_range  = ownershipRanges [comm.rank()+1] - ownershipRanges[comm.rank()];
    //
    //    std::cout<<"local_range_raw_range"<<local_range_raw_range<<std::endl;
    //
    //    const utopia::SizeType local_range_col_range = ownershipRanges[comm.rank()+1] - ownershipRanges[comm.rank()];
    //
    //    std::cout<<"local_range_col_range"<<local_range_col_range<<std::endl;
    //
    //    utopia::SizeType  mMaxRowEntries = mat_buffer.maxEntriesXCol();
    //
    //    comm.allReduce(&mMaxRowEntries, 1, express::MPIMax());
    //
    //    D = utopia::local_sparse(local_range_raw_range, local_range_col_range,mMaxRowEntries);
    //
    //    {
    //        utopia::Write<utopia::SparseMatT> write(D);
    //        for (auto it = mat_buffer.iter(); it; ++it) {
    //            //std::cout<<"*it"<<  *it <<std::endl;
    //            //if (*it!=0)
    //            //D_inv.set(it.row(), it.col(), 1./(*it));
    //            ///else
    //            D.set(it.row(), it.col(),*it);
    //            //B.set(it.row()+1, it.col()+1, *it);
    //
    //        }
    //    }
    //
    //
    //    utopia::write("D.m", D);
    
    
    //
    //    utopia::SparseMatT D_inv = utopia::local_sparse( matrix_A.m(), matrix_A.n(), 1);
    //
    //    {
    //        utopia::Write<utopia::SparseMatT> write(D);
    //        for (auto it = mat_buffer.iter(); it; ++it) {
    //            //std::cout<<"*it"<<  *it <<std::endl;
    //            if (*it!=0)
    //                D_inv.set(it.row(), it.col(), 1./(*it));
    //            else
    //                D_inv.set(it.row(), it.col(), 0.0);
    //            //B.set(it.row()+1, it.col()+1, *it);
    //
    //        }
    //    }
    //
    
}





