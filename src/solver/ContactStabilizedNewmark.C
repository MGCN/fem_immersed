#include <ContactStabilizedNewmark.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "utopia.hpp"
#include "FEProblemBase.h"
#include "FEProblem.h"
#include "MooseMesh.h"
//#include "utopia_assemble_contact.hpp"
#include "DisplacedProblem.h"
#include "Transient.h"
#include "TransientMultiApp.h"
#include "InputParameters.h"
#include <memory>
#include <utopia.hpp>
//#include "nonlinear_functions.h"

//#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_Newmark.hpp"
//#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"

PETSC_EXTERN PetscErrorCode SNESComputeBndStresses(SNES_MOOSEUTOPIA snes_MooseUtopia);

/* Our own solver, to be registered at runtime-solver */
#undef __FUNCT__
#define __FUNCT__ "SNESReset_ContactStabilizedNewmark"
PetscErrorCode SNESReset_ContactStabilizedNewmark(SNES /*snes*/)
{
    PetscFunctionBegin;
    PetscFunctionReturn(0);
}

/*
 SNESDestroyContactStabilizedNewmark - Destroys the private SNESContactStabilizedNewmark" context that was created with SNESCreateContactStabilizedNewmark"().
 
 Input Parameter:
 . snes - the SNES context
 
 Application Interface Routine: SNESDestroy()
 */
#undef __FUNCT__
#define __FUNCT__ "SNESDestroy_ContactStabilizedNewmark"
PetscErrorCode SNESDestroy_ContactStabilizedNewmark(SNES snes)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = SNESReset_ContactStabilizedNewmark(snes);CHKERRQ(ierr);
    
    
    ierr = PetscFree(snes->data);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/*
 SNESSetUpContactStabilizedNewmark - Sets up the internal data structures for the later use
 of the SNESSetUpContactStabilizedNewmark nonlinear solver.
 
 Input Parameters:
 +  snes - the SNES context
 -  x - the solution vector
 
 Application Interface Routine: SNESSetUp()
 */
#undef __FUNCT__
#define __FUNCT__ "SNESSetUp_ContactStabilizedNewmark"
PetscErrorCode SNESSetUp_ContactStabilizedNewmark(SNES /*snes*/)
{
    PetscFunctionBegin;
    PetscFunctionReturn(0);
}

/*
 SNESSetFromOptions_ContactStabilizedNewmark - Sets various parameters for the SNESlinearMGLS method.
 
 Input Parameter:
 . snes - the SNES context
 
 Application Interface Routine: SNESSetFromOptions()
 */
#undef __FUNCT__
#define __FUNCT__ "SNESSetFromOptions_ContactStabilizedNewmark"
static PetscErrorCode SNESSetFromOptions_ContactStabilizedNewmark(PetscOptionItems *PetscOptionsObject,SNES snes)
{
    PetscFunctionBegin;
    
    SNESSetFromOptions_MooseUtopia_options(PetscOptionsObject, snes);
    
    PetscFunctionReturn(0);
}

/*
 SNESViewContactStabilizedNewmark - Prints info from the SNESContactStabilizedNewmark data structure.
 
 Input Parameters:
 + SNES - the SNES context
 - viewer - visualization context
 
 Application Interface Routine: SNESView()
 */
#undef __FUNCT__
#define __FUNCT__ "SNESView_ContactStabilizedNewmark"
static PetscErrorCode SNESView_ContactStabilizedNewmark(SNES /*snes*/, PetscViewer viewer)
{
    PetscBool      iascii;
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERASCII, &iascii);CHKERRQ(ierr);
    if (iascii) {
    }
    PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "SNESSolve_ContactStabilizedNewmark"
PetscErrorCode SNESSolve_ContactStabilizedNewmark(SNES snes)
{
    using namespace utopia;
    typedef utopia::ContactStabilizedNewmark<utopia::USparseMatrix, utopia::UVector> ContactSolver;
    typedef utopia::LibMeshFunctionSpace V;
    typedef utopia::ProductFunctionSpace<V> VV;
    typedef utopia::FunctionMaterial<utopia::USparseMatrix, utopia::UVector> FunctionMaterialT;
    typedef PETSCUtopiaNonlinearFunctionMoose<utopia::USparseMatrix, utopia::UVector> FunctionT;
    //typedef PETSCUtopiaNonlinearFunctionNormalEq<utopia::USparseMatrix, utopia::UVector> FunctionT;
    utopia::UVector is_normal_component, zero_normals;
    USparseMatrix zero_normal_trafo;
    bool trafo;
    
    
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    
    
    
    {
        SNES_MOOSEUTOPIA *snes_MooseUtopia;
        PetscPrintf(PETSC_COMM_WORLD, "SNESSolve_ContactStabilizedNewmark: \n");
        
        snes_MooseUtopia = (SNES_MOOSEUTOPIA*)snes->data;

        Transient * ex = dynamic_cast<Transient *>(snes_MooseUtopia->_fe_problem->getMooseApp().getExecutioner());

        double dt = ex->getDT();
        
        auto equations_systems = make_ref(snes_MooseUtopia->_fe_problem->es());
        
        int system_num = snes_MooseUtopia->libmesh_system->number();
        
        unsigned int var_num_x = snes_MooseUtopia->libmesh_system->variable_number("disp_x");
        
        auto contact_problem = snes_MooseUtopia->_fe_problem;

        std::string userobject_name_1 = "Newmark";

        std::string userobject_name_2 = "CurrentTimeStep";
        
        if(!contact_problem->hasUserObject(userobject_name_1)){
            std::string class_name = "StoreTransferOperators";
            auto params = snes_MooseUtopia->_fe_problem->getMooseApp().getFactory().getValidParams(class_name);
            params.set<bool>("use_displaced_mesh") = false;
            params.set<ExecFlagEnum>("execute_on") = "initial";

            contact_problem->addUserObject("StoreTransferOperators", userobject_name_1, params);
            contact_problem->addUserObject("StoreTransferOperators", userobject_name_2, params);

        }
        
        

        auto &store_cs           = const_cast<StoreTransferOperators&>(contact_problem->getUserObject<StoreTransferOperators>(userobject_name_1));
        auto &store_current_time = const_cast<StoreTransferOperators&>(contact_problem->getUserObject<StoreTransferOperators>(userobject_name_2));

        if (!store_current_time.getVoidPointer()) {
            store_current_time.getVoidPointer() = std::make_shared<int>(contact_problem->timeStep());
        }
        
        int current_time_step = *std::static_pointer_cast<int>(store_current_time.getVoidPointer());
        

        std::cout << "newmark_time: " << contact_problem->timeStep() << ", solver_time: " << current_time_step << std::endl;

        bool time_step_changed = false;
        if(current_time_step != contact_problem->timeStep()) {
            time_step_changed = true;
            store_current_time.getVoidPointer() = std::make_shared<int>(contact_problem->timeStep());
        }

        std::shared_ptr<ContactSolver> cs;
        if (!store_cs.getVoidPointer()){
            
            const auto &b_slave = snes_MooseUtopia->_contact_slave_id;
            const auto &b_master = snes_MooseUtopia->_contact_master_id;

            assert(b_slave.size() == b_master.size());
            
            if(b_slave.size() != b_master.size()) {
                std::cerr << "[Error] contact pairs do not have same size: " << b_slave.size() << " != " << b_master.size() << std::endl;
            }

            ContactParams contact_params;
            contact_params.search_radius = snes_MooseUtopia->_search_radius;
            contact_params.contact_pair_tags;

            for(std::size_t i = 0; i < b_slave.size(); ++i) {
                std::cout <<"Contacts IDs= " << b_master[i] << ", " << b_slave[i] << std::endl;
                contact_params.contact_pair_tags.push_back(std::make_pair(int(b_master[i]), int(b_slave[i])));
            }
  
            auto vv = std::make_shared<VV>();
            
            int dims = snes_MooseUtopia->libmesh_system->get_mesh().mesh_dimension();
            
            for(int i = 0; i < dims; ++i) {
                vv->add_subspace(std::make_shared<V>(equations_systems, system_num, var_num_x + i));
            }
            
           // auto function = std::make_shared<FunctionT>(snes);
            MeshBase &mesh=snes_MooseUtopia->libmesh_system->get_mesh();
            
            DofMap & dof_map = snes_MooseUtopia->libmesh_system->get_dof_map();
            
            std::cout<<"density_value=> "<<snes_MooseUtopia->_density_value<<std::endl;

            auto function = std::make_shared<FunctionT>(snes);


            auto material = std::make_shared<FunctionMaterialT>(function);
            
            cs =  std::make_shared<ContactSolver>(vv, material, dt, contact_params);
 
            MooseVariable & to_var_x = snes_MooseUtopia->_fe_problem->getStandardVariable(0,"disp_x");

            MooseVariable & to_var_y = snes_MooseUtopia->_fe_problem->getStandardVariable(0,"disp_y");

            utopia::UVector initial_data;

            initial_data=local_zeros(dof_map.n_local_dofs());

            System & to_sys = to_var_x.sys().system();

            if(snes_MooseUtopia->_fe_problem->timeStep()==1 && snes_MooseUtopia->init_value!=0.0)
            {

                std::cout<<"snes_MooseUtopia->_init_value_x==>"<<snes_MooseUtopia->init_value<<std::endl;
                //std::cout<<"snes_MooseUtopia->_init_value_y==>"<<snes_MooseUtopia->init_value_y<<std::endl;

                    MooseVariable & var = snes_MooseUtopia->_fe_problem->getStandardVariable(0, "disp_x");

                    System & sys = var.sys().system();

                    Vec vector = static_cast<libMesh::PetscVector<libMesh::Number> *>(sys.current_local_solution.get())->vec();

                    utopia::convert(vector, initial_data);

                    //utopia::disp(initial_data);

                    cs->initial_condition(snes_MooseUtopia->_density_value,initial_data);
            }

            else
            {

             	cs->initial_condition(snes_MooseUtopia->_density_value);

            }
            
            //cs->set_use_ssn(true);
            cs->initialize();

            cs->set_tol(snes_MooseUtopia->tol_solver);
            
            //cs->set_tol(1e-3);

            const bool use_mg = false;
           
            
            if(use_mg) {
                //begin: multigrid
                cs->set_tol(5e-6);
                auto linear_solver = std::make_shared<BiCGStab<utopia::USparseMatrix, utopia::UVector>>();
                auto smoother = std::make_shared<ProjectedGaussSeidel<utopia::USparseMatrix, utopia::UVector> >();
                auto mg = std::make_shared<SemiGeometricMultigrid>(smoother, linear_solver);
                
                mg->verbose(true);
                mg->init(vv->subspace(0), 4);
                
                mg->algebraic().atol(1e-11);
                mg->algebraic().rtol(1e-11);
                mg->algebraic().stol(1e-11);

                cs->set_linear_solver(mg);
                //end: multigrid
            }
            
            store_cs.getVoidPointer() = cs;
            
        } else {
            cs = std::static_pointer_cast<ContactSolver>(store_cs.getVoidPointer());
            MeshBase &mesh=snes_MooseUtopia->libmesh_system->get_mesh();
            DofMap & dof_map = snes_MooseUtopia->libmesh_system->get_dof_map();
            auto function = std::make_shared<FunctionT>(snes);
            auto material = std::make_shared<FunctionMaterialT>(function);
            cs->set_material(material);
        }

        if(time_step_changed) {
            cs->next_step();
            cs->reset();
        }

        cs->set_dt_and_update(dt);
        //solve only for increment
        cs->solve_contact();
        utopia::UVector sol_c;
       
        if(trafo){
            std::cout<<"Apply Transpose Orthogonal Trnsformation"<<'\n';
            sol_c = zero_normal_trafo * cs->displacement();
        }
        else{
            sol_c = cs->displacement();
            
        }
        convert(sol_c, snes->vec_sol);
    }
    
    snes->reason = SNES_CONVERGED_FNORM_ABS;
    
    
    PetscPrintf(PETSC_COMM_WORLD, "End of SNESSolve_ContactStabilizedNewmark: \n");
    PetscFunctionReturn(0);
    
}
//
#undef __FUNCT__
#define __FUNCT__ "SNESCreate_ContactStabilizedNewmark"
PETSC_EXTERN PetscErrorCode SNESCreate_ContactStabilizedNewmark(SNES snes)
{
    PetscErrorCode   ierr;
    SNES_MOOSEUTOPIA *neP;
    
    PetscFunctionBegin;
    snes->ops->destroy        = SNESDestroy_ContactStabilizedNewmark;
    snes->ops->setup          = SNESSetUp_ContactStabilizedNewmark;
    snes->ops->setfromoptions = SNESSetFromOptions_ContactStabilizedNewmark;
    snes->ops->view           = SNESView_ContactStabilizedNewmark;
    snes->ops->solve          = SNESSolve_ContactStabilizedNewmark;
    snes->ops->reset          = SNESReset_ContactStabilizedNewmark;
    
    snes->usesksp = PETSC_FALSE;
    snes->usesnpc  = PETSC_TRUE;
    
    //snes->pcside = PC_LEFT;
    
    ierr       = PetscNewLog(snes,&neP);CHKERRQ(ierr);
    snes->data = (void*) neP;
    
    snes->reason = SNES_CONVERGED_ITS;
    SNESMonitorSet(snes,MooseUtopiaMonitor,neP,0);
    
    PetscPrintf(PETSC_COMM_WORLD, "In Function SNESCreate_ContactStabilizedNewmark \n");
    PetscFunctionReturn(0);
}

