#ifndef MooseUtopia_Solver_h
#define MooseUtopia_Solver_h

#include <numeric>
#include <string>
#include <iostream>
#include <vector>
#include <petscsnes.h>
#include <petsc/private/snesimpl.h>
#include <Eigen/Dense>
#include <utopia.hpp>

#include <libmesh/mesh_base.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/system.h>
#include <libmesh/dof_object.h>
#include <libmesh/node.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include "FEProblem.h"



PETSC_EXTERN PetscErrorCode SNESCreate_ContactStabilizedNewmark(SNES snes);




#define MOOSEUTOPIA_SOLVER_INFTY 1e30
#define MOOSEUTOPIA_SOLVER_MINFTY -1e30
#define MOOSEUTOPIA_SOLVER_ASSOTOL 1e-12
#define MOOSEUTOPIA_SOLVER_MIN 1e-30

// extern PetscErrorCode v_nan_value(Vec v, PetscScalar replace);
// extern PetscErrorCode v_value(Vec v, PetscScalar replace);
// extern PetscErrorCode xtAx(Mat A, Vec x, PetscScalar * _enorm);
// extern PetscErrorCode energy_norm(Mat A, Vec f, Vec c, PetscScalar * _enorm);
// extern PetscErrorCode lin_residual(Mat A, Vec b, Vec x, Vec res);
// extern PetscErrorCode passo_linesearch(Vec g, Mat A, Vec c, PetscScalar * alpha);
// extern PetscErrorCode inactive_set_residual_norm(Mat A, Vec f, Vec c, Vec active_set, PetscScalar * norm);

extern PetscErrorCode MooseUtopiaMonitor(SNES, PetscInt, PetscReal, void *);
// extern PetscErrorCode GetSolveInfo(SNES &snes, const utopia::SolutionStatus &params);
// extern PetscErrorCode SNESConvergencePasso(SNES snes, PetscInt it,
//                                            PetscReal xnorm, PetscReal snorm,
//                                            PetscReal fnorm,
//                                            SNESConvergedReason * reason,
//                                            void * dummy);

PETSC_EXTERN  PetscErrorCode SNESSetFromOptions_MooseUtopia_options(PetscOptionItems *PetscOptionsObject,SNES snes);

#define MOOSEUTOPIA_SOLVERS  \
"stabilized_newmark"


typedef struct{
  
  SNES parent_snes;
  Vec parent_form_fct;
  Mat parent_form_jac;
  
  IS child_range;
  
  /* hack for linear elasticity */
  Vec tmp_f;
  Mat tmp_J;
  
} CONTACT_CTX;


/*
 passo solver parameters
 */
typedef struct
{
    /* general */
    PetscBool verbose;
    PetscBool ls_verbose;
    PetscBool time_statistics;
    
    PetscReal tol;  // to be deleted
    
    PetscReal abs_tol;
    PetscReal rtol;
    PetscReal stol;
    
    PetscInt max_it = 100;
    
    bool normal_eq = true; // determine which kind of minimization is happening
    
    /* multigrid */
    PetscInt mg_presmoothing_steps = 1;
    PetscInt mg_postsmoothing_steps = 1;
    PetscInt mg_type; /*!< V, W  cycle  */
    PetscBool mg_native;
    char smoother_type[256];
    char coarse_solve_type[256];
    
    /* linesearch parameters */
    char ls_alg[256]; /*!< choice of LS alg  */
    PetscReal ls_rho;    /*!<   */
    PetscReal c1;    /*!<   */
    PetscReal c2;    /*!<   */
    
    /* trust region parameters */
    char tr_alg[256];   /*!< choice of TR alg  */
    char path_to_log[256]; 
    char consistency_level[256]; 

    
    PetscReal delta_max; /*!< max TR radius delta  */
    PetscReal delta0;    /*!< initial TR radius  */
    
    PetscReal gamma1; /*!< amount of shrinking TR */
    PetscReal gamma2; /*!< amount of enlarging TR */
    
    PetscReal eta1; /*!< amount of shrinking TR */
    PetscReal eta2; /*!< amount of enlarging TR  */
    
    PetscReal rho_tol;
    PetscReal ST_tol;
    
    PetscInt local_max_it; /*!< amount of enlarging TR  */
    PetscInt overlap; /*!< amount of enlarging TR  */
    
    // here we are mostly following PETSc standard
    char ksp_type[256];
    char pc_type[256];
    char pc_factor_mat_solver_package[256];
    
    PetscReal ksp_atol;
    PetscReal ksp_rtol;
    PetscReal ksp_dtol;
    PetscReal ksp_max_it;
    PetscBool linear_solver_verbose;
    
    PetscReal nl_atol;  /*!< absolute tolerance for residual in nonlinear solve  */
    
    /* datatypes for logging to matlab using PETSc */
    PetscViewer viewer_iterates;
    Vec * iterates;
    PetscBool MooseUtopia_log_iterates = PETSC_FALSE;
    
    PetscViewer viewer_system;
    PetscBool MosseUtopia_log_system = PETSC_FALSE;
    
    PetscViewer viewer_norms;
    Vec fnorm_out, snorm_out, enorm_out, rate_out;
    PetscBool MooseUtopia_log_norms = PETSC_FALSE;
    
    bool initialized = false;
    PetscBool static_time_step;

    PetscInt hessian_update = 1 ; 

    PetscInt blocksize;  /*!< Blocksize for block Gauss-Seidel Smoother  */
    PetscReal omega; /*!< not used: For doing partial updates in NBGS  */
} MOOSEUTOPIA_SOLVERS_PARAMS;

/* ctx for logging */
typedef struct
{
    PetscScalar fnorm; // |r_{k+1}|
    PetscScalar snorm; // |x_{k+1} - x_k|
    PetscScalar enorm; // Linear case: x'Ax - f'x
    PetscScalar rate;  // |x_{k+1} - x_k|/|x_{k} - x_{k-1}|
    PetscScalar alpha; // for the alpha that is set when using nbgs with LS
    PetscScalar active_coords;
    PetscScalar changed_coords;
    
    PetscInt dof=-1;
    PetscInt procs=-1;

    PetscBool log_to_file;
    std::string * logfile_base;
    
    FILE *logfile;
    
}MOOSEUTOPIA_LOG;





/* the snes struct for all solvers */
typedef struct
{
    // Mat A = nullptr; // for convenience in case of field splitting
    // Mat B = nullptr; // for convenience in case of field splitting
    // Mat C = nullptr; // for convenience in case of field splitting
    // Mat D = nullptr; // for convenience in case of field splitting
    
    // Vec b = nullptr; // right-hand-side
    Vec solution = nullptr; // the a priori known solution. This is for testing.
    // Vec nodal_coordinates = nullptr;
    
    // // variables for constraints
    // Vec nodal_normals = nullptr;
    
    libMesh::System * libmesh_system;
    libMesh::System * libmesh_aux_system;
    FEProblem * _fe_problem;

    std::vector<BoundaryID> _contact_master_id;
    std::vector<BoundaryID> _contact_slave_id;
    double _search_radius;
    double _density_value;
    double init_value;
    double tol_solver;


    //bool init_value_bool;
  
//    
    std::vector<std::string> variables; // names of solution variables. This can sometimes be a subset of the variables libmesh sets up.
    
    MOOSEUTOPIA_SOLVERS_PARAMS MooseUtopia_params;

  
    std::vector<Mat> interpolations;    // interpolation operator for multigrid in matrix form
    std::vector<Mat> projections_down;  // projections to project interates down to the hierarchy 
    std::vector<SNES> sub_sness;        // subsnes for local solves/ multilevel nonlinear solve 


    std::vector<Vec> init_guesses;              // vector of initial guesses for sub-solvers - could be potentially extracted from sub_sness
    std::vector<Vec> bc_dirichlet_marker;       // boolean vector of BC conditions


    // vectors for pointwise constraints 
    Vec upper_bound =   nullptr;             
    Vec lower_bound =   nullptr; 

    Vec fracture_marker =   nullptr; 

     /* place holder for contact variables */
    Mat Z,Zinv;     // Transformation operator as described in Dickopf/Krause 2009
    Mat O;
    Mat Dinv;
    Vec normals;    // A node*dimension vector containg the normal vectors on slave side
    Vec gap;        // the gap in normal direction from the slave side
    Vec lb;
    Vec ub;


       
    MOOSEUTOPIA_LOG MooseUtopia_log;
    PetscBool initialized;

    
} SNES_MOOSEUTOPIA;

extern PetscErrorCode parameters_to_utopia(const SNES_MOOSEUTOPIA &snes_MooseUtopia, utopia::InputParameters & utopia_params);
// extern PetscErrorCode determine_dc_bnd_var_id(const std::vector<std::string> & _variables, std::vector<int> _dc_boundary_id,                                libMesh::System * _sys, std::vector<std::vector<int> > * _dc_variables_id);
// extern PetscErrorCode splice_interpolation(Mat In, Mat * Out, int dim);
// extern PetscErrorCode reduce_interpolation(Mat, Mat * , std::vector<std::string>, libMesh::System *);
// extern PetscErrorCode dirichletnodes(Vec *, const std::vector<int> _dc_boundary_id, libMesh::System *);
// extern PetscErrorCode householder_block_transformation(Vec _normals, int _blocksize, Mat _template, Mat * _Orth);
// extern PetscErrorCode mat_truncate(Mat, PetscScalar);
// extern PetscErrorCode remove_singularities(Mat * matrix);
// extern PetscErrorCode remove_zero_rows(Mat * matrix);
// extern PetscErrorCode remove_zero_cols(Mat * matrix);

//extern PetscErrorCode constrain_blocknd(MatScalar * diag, PetscScalar * t, PetscScalar * x, PetscScalar * lb, PetscScalar * ub, int row, int bs, PetscBool * constrained);
//extern PetscErrorCode constrain_block_xl(MatScalar * diag, PetscScalar * t, PetscScalar * x, PetscScalar * xl, PetscScalar * , int row, int bs, PetscBool * constrained);
//extern PetscErrorCode constrain_box(BLK_SOLVER_CTX blk_ctx, MatScalar * diag, PetscScalar * t, PetscScalar * x, PetscScalar * lb, PetscScalar * ub, int row, int bs, PetscBool * constrained);
//
//extern PetscErrorCode  constrain_block_xl(BLK_SOLVER_CTX blk_ctx, MatScalar * diag, PetscScalar * t, PetscScalar * x,
//                                   PetscScalar * _xl, PetscScalar * ub, int row,
//                                   int bs, PetscBool * constrained);
//
//extern PetscErrorCode constrain_blocknd(BLK_SOLVER_CTX blk_ctx, MatScalar * diag, PetscScalar * t, PetscScalar * x,
//                                 PetscScalar * _xl, PetscScalar * ub, int row,
//                                 int bs, PetscBool * constrained);

#endif

