#include <MooseUtopia_Solver.h>
#include <math.h>

#include "petscis.h"
#include <petsc/private/kernels/blockinvert.h>

#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"

#undef __FUNCT__
#define __FUNCT__ "v_nan_value"
extern PetscErrorCode v_nan_value(Vec v, PetscScalar replace)
{
    PetscErrorCode ierr;
    PetscInt M,m;
    int rank;
    PetscScalar *values;
    
    PetscFunctionBegin;
    
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    
    ierr = VecGetSize(v, &M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(v, &m); CHKERRQ(ierr);
    
    ierr = VecGetArray(v,&values);CHKERRQ(ierr);
    for (int i =0;i<m;i++){
        if (std::isnan(values[i])){
            values[i] = replace;
        }
    }
    ierr = VecRestoreArray(v, &values);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// rset vector to one value
#undef __FUNCT__
#define __FUNCT__ "v_value"
extern PetscErrorCode v_value(Vec v, PetscScalar replace)
{
    PetscErrorCode ierr;
    PetscInt M,m;
    int rank;
    PetscScalar *values;
    
    PetscFunctionBegin;
    
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    ierr = VecGetSize(v, &M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(v, &m); CHKERRQ(ierr);
    ierr = VecGetArray(v,&values);CHKERRQ(ierr);
    for (int i =0;i<m;i++){
        values[i] = replace;
    }
    ierr = VecRestoreArray(v, &values); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MooseUtopiaMonitor"
PetscErrorCode MooseUtopiaMonitor(SNES snes,PetscInt k,PetscReal fnorm, void * ctx)
{
    int rank;
    
    PetscFunctionBegin;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    if (ctx){
        MOOSEUTOPIA_LOG * MooseUtopia_log;
        MooseUtopia_log = (MOOSEUTOPIA_LOG *) ctx;
        
        if (k == 0){
            PetscPrintf(PETSC_COMM_WORLD,"MooseUtopiaMonitor: \n");
            PetscPrintf(PETSC_COMM_WORLD,"it, energy, ||res||, ||step||, rate, alpha, active indices, changed indices \n");
            PetscPrintf(PETSC_COMM_WORLD,"%i, %2.12e, %2.12e, %2.12e, %2.12e, %2.12e, %i, %i \n", k, MooseUtopia_log->enorm, MooseUtopia_log->fnorm, MooseUtopia_log->snorm, MooseUtopia_log->rate, MooseUtopia_log->alpha, (int)MooseUtopia_log->active_coords, (int)MooseUtopia_log->changed_coords);

            
            if (rank == 0){
                if (MooseUtopia_log->log_to_file){
                    // datetime filename
                    char date_char[25];
                    char filename[70];
                    
                    struct tm *tt;
                    time_t rawtime;
                    time(&rawtime);
                    tt = gmtime(&rawtime);
                    strftime(date_char, sizeof(date_char), "%Y%m%d_%H_%M_%S.txt", tt);
                    
                    std::copy(MooseUtopia_log->logfile_base->begin(), MooseUtopia_log->logfile_base->end(), filename);
                    filename[MooseUtopia_log->logfile_base->size()] = '\0';
                    strcat(filename, date_char);
                    PetscFOpen(PETSC_COMM_SELF,filename,"w",&(MooseUtopia_log->logfile));



                    if (!MooseUtopia_log->logfile) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
                    PetscFPrintf(PETSC_COMM_SELF, MooseUtopia_log->logfile,"dof= %i, nprocs=%i \n", MooseUtopia_log->dof, MooseUtopia_log->procs, MooseUtopia_log->alpha);
                    PetscFPrintf(PETSC_COMM_SELF, MooseUtopia_log->logfile,"it, energy, ||res||, ||step||, rate, alpha, active indices, changed indices \n");
                    PetscFPrintf(PETSC_COMM_SELF, MooseUtopia_log->logfile, "%i, %2.12e, %2.12e, %2.12e, %2.12e , %2.12e, %i, %i \n",
                                 k, MooseUtopia_log->enorm, MooseUtopia_log->fnorm, MooseUtopia_log->snorm, MooseUtopia_log->rate, (int)MooseUtopia_log->active_coords, (int)MooseUtopia_log->changed_coords);
                }
            }
            
        }else{
            PetscPrintf(PETSC_COMM_WORLD,"%i, %2.12e, %2.12e, %2.12e, %2.12e, %2.12e, %i, %i \n", k, MooseUtopia_log->enorm, MooseUtopia_log->fnorm, MooseUtopia_log->snorm, MooseUtopia_log->rate, MooseUtopia_log->alpha, (int)MooseUtopia_log->active_coords, (int)MooseUtopia_log->changed_coords);
            if (MooseUtopia_log->logfile && rank == 0 && MooseUtopia_log->log_to_file){
                PetscFPrintf(PETSC_COMM_SELF, MooseUtopia_log->logfile, "%i, %2.12e, %2.12e, %2.12e, %2.12e , %2.12e, %i, %i  \n", k, MooseUtopia_log->enorm, MooseUtopia_log->fnorm, MooseUtopia_log->snorm, MooseUtopia_log->rate, MooseUtopia_log->alpha, (int)MooseUtopia_log->active_coords, (int)MooseUtopia_log->changed_coords);
            }
        }
    }else{
        if (k == 0){
            
            PetscPrintf(PETSC_COMM_WORLD,"MooseUtopiaMonitor: \n");
            PetscPrintf(PETSC_COMM_WORLD,"it,||res|| \n");
            PetscPrintf(PETSC_COMM_WORLD,"%i, %2.12e  \n", k, fnorm);
            
        }else{
            PetscPrintf(PETSC_COMM_WORLD,"%i, %2.12e  \n", k, fnorm);
        }
    }
    PetscFunctionReturn(0);
}

/*
 GetSolveInfo - assignning info from utopia to PETSc snes
 - required for interaction with MOOSE
 
 Input Parameters:
 . snes - the SNES context
 . info - the output from utopia solve routine
 
 Application Interface Routine: SNESSolve()
 */
PetscErrorCode GetSolveInfo(SNES &snes, const utopia::SolutionStatus &status)
{
    
    switch ((PetscInt)status.reason)
    {
            // sucess
        case utopia::ConvergenceReason::CONVERGED_FNORM_ABS :
            snes->reason = SNES_CONVERGED_FNORM_ABS;
            break;
            
        case utopia::ConvergenceReason::CONVERGED_FNORM_RELATIVE :
            snes->reason = SNES_CONVERGED_FNORM_RELATIVE;
            break;
            
        case utopia::ConvergenceReason::CONVERGED_SNORM_RELATIVE :
            snes->reason = SNES_CONVERGED_SNORM_RELATIVE;
            break;
            
        case utopia::ConvergenceReason::CONVERGED_ITS :
            snes->reason = SNES_CONVERGED_ITS;
            break;
            
        case utopia::ConvergenceReason::CONVERGED_TR_DELTA :
            snes->reason = SNES_CONVERGED_TR_DELTA;
            break;
            
            // fail
        case utopia::ConvergenceReason::DIVERGED_FUNCTION_DOMAIN :
            snes->reason = SNES_DIVERGED_FUNCTION_DOMAIN;
            break;
            
        case utopia::ConvergenceReason::DIVERGED_FUNCTION_COUNT :
            snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
            break;
            
        case utopia::ConvergenceReason::DIVERGED_LINEAR_SOLVE :
            snes->reason = SNES_DIVERGED_LINEAR_SOLVE;
            break;
            
        case utopia::ConvergenceReason::DIVERGED_FNORM_NAN :
            snes->reason = SNES_DIVERGED_FNORM_NAN;
            break;
            
        case utopia::ConvergenceReason::DIVERGED_MAX_IT :
            snes->reason = SNES_DIVERGED_MAX_IT;
            break;
            
        case utopia::ConvergenceReason::DIVERGED_LINE_SEARCH :
            snes->reason = SNES_DIVERGED_LINE_SEARCH;
            break;
            
        case utopia::ConvergenceReason::DIVERGED_INNER :
            snes->reason = SNES_DIVERGED_INNER;
            break;
            
        case utopia::ConvergenceReason::DIVERGED_LOCAL_MIN :
            snes->reason = SNES_DIVERGED_LOCAL_MIN;
            break;
            
        case utopia::ConvergenceReason::CONVERGED_ITERATING :
            snes->reason = SNES_CONVERGED_ITERATING;
            break;
            
        default :
            snes->reason = SNES_CONVERGED_FNORM_ABS;
    }
    
    snes->iter = status.iterates;
    
    return 0;
}


/*
 GetSolveInfo - assignning info from utopia to PETSc snes
 - required for interaction with MOOSE
 
 Input Parameters:
 . snes - the SNES context
 . info - the output from utopia solve routine
 
 Application Interface Routine: SNESSolve()
 */
PetscErrorCode parameters_to_utopia(const SNES_MOOSEUTOPIA &snes_MooseUtopia, utopia::InputParameters & utopia_params)
{
    
    const MOOSEUTOPIA_SOLVERS_PARAMS& MooseUtopia_params = snes_MooseUtopia.MooseUtopia_params;


    //            /*    general      */
    if(MooseUtopia_params.verbose!= 0)
        utopia_params.set("verbose", MooseUtopia_params.verbose);    
    if(MooseUtopia_params.max_it!= 0)
        utopia_params.set("max_it", MooseUtopia_params.max_it);
    
    if(MooseUtopia_params.abs_tol!= 0)
        utopia_params.set("atol", MooseUtopia_params.abs_tol);
    if(MooseUtopia_params.rtol!= 0)
        utopia_params.set("rtol", MooseUtopia_params.rtol);
    if(MooseUtopia_params.stol!= 0)
        utopia_params.set("stol", MooseUtopia_params.stol);
    
     //           /*     trust_region  */
     if(MooseUtopia_params.rho_tol!= 0)
        utopia_params.set("rho_tol", MooseUtopia_params.rho_tol);
    if(MooseUtopia_params.delta_max!= 0)
        utopia_params.set("delta_max", MooseUtopia_params.delta_max);
    if(MooseUtopia_params.delta0!= 0)
        utopia_params.set("delta0", MooseUtopia_params.delta0);
    if(MooseUtopia_params.gamma1!= 0)
        utopia_params.set("gamma1", MooseUtopia_params.gamma1);
    if(MooseUtopia_params.gamma2!= 0)
        utopia_params.set("gamma2", MooseUtopia_params.gamma2);
    if(MooseUtopia_params.eta1!= 0)
        utopia_params.set("eta1", MooseUtopia_params.eta1);
    if(MooseUtopia_params.eta2!= 0)
        utopia_params.set("eta2", MooseUtopia_params.eta2);

    //            /*   multigrid    */

    if(MooseUtopia_params.mg_presmoothing_steps!= 0)
        utopia_params.set("presmoothing_steps", MooseUtopia_params.mg_presmoothing_steps);
    if(MooseUtopia_params.mg_postsmoothing_steps!= 0)
        utopia_params.set("postsmoothing_steps", MooseUtopia_params.mg_postsmoothing_steps);
    if(MooseUtopia_params.mg_type!= 0)
        utopia_params.set("mg_type", MooseUtopia_params.mg_type);
 
    //            /*    linesearch    */
    if(MooseUtopia_params.ls_rho!= 0)
        utopia_params.set("rho", MooseUtopia_params.ls_rho);        
    if(MooseUtopia_params.c1!= 0)
        utopia_params.set("c1", MooseUtopia_params.c1);
    if(MooseUtopia_params.c2!= 0)
        utopia_params.set("c2", MooseUtopia_params.c2);


  /* if(MooseUtopia_params.verbose!= 0)
        utopia_params.verbose(MooseUtopia_params.verbose);
    if(MooseUtopia_params.time_statistics!= 0)
        utopia_params.time_statistics(MooseUtopia_params.time_statistics);
    
    if(MooseUtopia_params.max_it!= 0)
        utopia_params.max_it(MooseUtopia_params.max_it);
    
    if(MooseUtopia_params.abs_tol!= 0)
        utopia_params.atol(MooseUtopia_params.abs_tol);
    if(MooseUtopia_params.rtol!= 0)
        utopia_params.rtol(MooseUtopia_params.rtol);
    if(MooseUtopia_params.stol!= 0)
        utopia_params.stol(MooseUtopia_params.stol);
    
    
   
    // if(MooseUtopia_params.tr_alg!= 0)
    utopia_params.trust_region_alg(MooseUtopia_params.tr_alg);
    if(MooseUtopia_params.rho_tol!= 0)
        utopia_params.rho_tol(MooseUtopia_params.rho_tol);
    if(MooseUtopia_params.ST_tol!= 0)
        utopia_params.SteihaugToint_tol(MooseUtopia_params.ST_tol);
    if(MooseUtopia_params.delta_max!= 0)
        utopia_params.delta_max(MooseUtopia_params.delta_max);
    if(MooseUtopia_params.delta0!= 0)
        utopia_params.delta0(MooseUtopia_params.delta0);
    if(MooseUtopia_params.gamma1!= 0)
        utopia_params.gamma1(MooseUtopia_params.gamma1);
    if(MooseUtopia_params.gamma2!= 0)
        utopia_params.gamma2(MooseUtopia_params.gamma2);
    if(MooseUtopia_params.eta1!= 0)
        utopia_params.eta1(MooseUtopia_params.eta1);
    if(MooseUtopia_params.eta2!= 0)
        utopia_params.eta2(MooseUtopia_params.eta2);
    
    
    //    if(snes_passo.nbgs.blocksize!= 0)
    //        utopia_params.block_size(snes_passo.nbgs.blocksize);
    //    if(snes_passo.nbgs.omega!= 0!= 0)
    //        utopia_params.omega(snes_passo.nbgs.omega);
    if(MooseUtopia_params.mg_presmoothing_steps!= 0)
        utopia_params.pre_smoothing_steps(MooseUtopia_params.mg_presmoothing_steps);
    if(MooseUtopia_params.mg_postsmoothing_steps!= 0)
        utopia_params.post_smoothing_steps(MooseUtopia_params.mg_postsmoothing_steps);
    if(MooseUtopia_params.mg_type!= 0)
        utopia_params.mg_type(MooseUtopia_params.mg_type);
    if(MooseUtopia_params.smoother_type!= 0)
        utopia_params.smoother_type(MooseUtopia_params.smoother_type);
    
    
    if(MooseUtopia_params.overlap!= 0)
        utopia_params.overlap(MooseUtopia_params.overlap);
    if(MooseUtopia_params.local_max_it!= 0)
        utopia_params.local_max_it(MooseUtopia_params.local_max_it);
    
   
    utopia_params.line_search_alg(MooseUtopia_params.ls_alg);
    if(MooseUtopia_params.ls_rho!= 0)
        utopia_params.ls_rho(MooseUtopia_params.ls_rho);
    if(MooseUtopia_params.c1!= 0)
        utopia_params.c1(MooseUtopia_params.c1);
    if(MooseUtopia_params.c2!= 0)
        utopia_params.c2(MooseUtopia_params.c2);
    
    
    // TODO:: add LS parameters .....
    if(MooseUtopia_params.ls_verbose!= 0)
        utopia_params.linear_solver_verbose(MooseUtopia_params.ls_verbose);
    
    if(MooseUtopia_params.linear_solver_verbose!= 0)
        utopia_params.linear_solver_verbose(MooseUtopia_params.linear_solver_verbose);
    
    if(MooseUtopia_params.ksp_atol!= 0)
        utopia_params.ksp_atol(MooseUtopia_params.ksp_atol);
    if(MooseUtopia_params.ksp_rtol!= 0)
        utopia_params.ksp_rtol(MooseUtopia_params.ksp_rtol);
    if(MooseUtopia_params.ksp_dtol!= 0)
        utopia_params.ksp_dtol(MooseUtopia_params.ksp_dtol);
    if(MooseUtopia_params.ksp_max_it!= 0)
        utopia_params.ksp_max_it(MooseUtopia_params.ksp_max_it);
    
    utopia_params.lin_solver_type(MooseUtopia_params.ksp_type);
    
    utopia_params.preconditioner_type(MooseUtopia_params.pc_type);
    utopia_params.preconditioner_factor_mat_solver_package(MooseUtopia_params.pc_factor_mat_solver_package);
   */ 
    return 0;
}

// #undef __FUNCT__
// #define __FUNCT__ "xtAx"
// PetscErrorCode xtAx(Mat A, Vec x, PetscScalar * _enorm)
// {
//     PetscErrorCode ierr;
//     Vec tmp;
    
//     PetscFunctionBegin;
    
//     ierr = VecDuplicate(x, &tmp); CHKERRQ(ierr);
//     ierr = MatMult(A, x, tmp); CHKERRQ(ierr);
//     ierr = VecDot(x,tmp, _enorm); CHKERRQ(ierr);
    
    
//     ierr = VecDestroy(&tmp);CHKERRQ(ierr);
    
//     PetscFunctionReturn(0);
// }

// #undef __FUNCT__
// #define __FUNCT__ "lin_residual"
// PetscErrorCode lin_residual(Mat A, Vec b, Vec x, Vec res)
// {
//     PetscErrorCode ierr;
//     PetscFunctionBegin;
    
//     ierr = VecScale(x, -1); CHKERRQ(ierr);
//     ierr = MatMultAdd(A, x, b, res); CHKERRQ(ierr);
//     ierr = VecScale(x, -1);  CHKERRQ(ierr);
    
//     PetscFunctionReturn(0);
// }

// #undef __FUNCT__
// #define __FUNCT__ "mat_truncate"
// extern PetscErrorCode mat_truncate(Mat _M, PetscScalar _epsilon){
    
//     Mat M_loc;
//     PetscInt rank, size;
//     PetscErrorCode ierr;
    
//     PetscFunctionBegin;
    
//     MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
//     MPI_Comm_size(PETSC_COMM_WORLD,&size);
    
//     if(size>1){
//         ierr = MatMPIAIJGetLocalMat(_M,MAT_INITIAL_MATRIX,&M_loc); CHKERRQ(ierr);
//     }else{
//         M_loc = _M;
//     }
    
//     Mat_SeqAIJ * mat_header   = (Mat_SeqAIJ *)M_loc->data;
    
//     for (int i=0;i< mat_header->maxnz;i++){
        
//         if (fabs(mat_header->a[i]) < _epsilon){
//             mat_header->a[i] = 0;
//         };
//     }
    
//     PetscFunctionReturn(0);
// }



// #undef __FUNCT__
// #define __FUNCT__ "splice_interpolation"
// PetscErrorCode splice_interpolation(Mat In, Mat * Out, int dim){
    
//     PetscErrorCode  ierr;
//     PetscInt        rank, size;
//     Mat TMP;
    
//     PetscFunctionBegin;
//     std::chrono::high_resolution_clock::time_point _t_start = std::chrono::high_resolution_clock::now();
    
//     MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
//     MPI_Comm_size(PETSC_COMM_WORLD,&size);
    
//     // memory-leak ?!?
//     ierr = MatCreateMAIJ(In,dim,&TMP); CHKERRQ(ierr);
    
//     if(size>1){
//         ierr = MatConvert(TMP, MATMPIAIJ, MAT_INITIAL_MATRIX, Out);CHKERRQ(ierr);
//     }else{
//         ierr = MatConvert(TMP, MATAIJ, MAT_INITIAL_MATRIX, Out);CHKERRQ(ierr);
//     }
    
//     ierr = MatDestroy(&TMP); CHKERRQ(ierr);
//     std::chrono::high_resolution_clock::time_point _t_end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double, std::milli> fp_ms = _t_end - _t_start ;
//     PetscPrintf(PETSC_COMM_WORLD, "  Walltime of splice_interpolation: %f ms. \n", fp_ms.count());
    
//     PetscFunctionReturn(0);
// };


// #undef __FUNCT__
// #define __FUNCT__ "reduce_interpolation"
// extern PetscErrorCode reduce_interpolation(Mat _I, Mat * _I_red, std::vector<std::string> _variable_names, libMesh::System * _sys){
//     PetscFunctionBegin;
    
//     PetscErrorCode ierr;
//     IS is_r, is_c;
    
//     int length_r, length_c;
//     int m,n,N;
//     int first, last;
//     int first_c, last_c;
//     int * indices_r;
//     int * indices_c;
//     int dim, new_dim;
//     int nodes_h, nodes_H;
    
//     std::chrono::high_resolution_clock::time_point _t_start;
//     std::chrono::high_resolution_clock::time_point _t_end;
//     std::chrono::duration<double, std::milli> fp_ms;
    
    
//     ierr = MatGetSize(_I, NULL, &N); CHKERRQ(ierr);
//     ierr = MatGetLocalSize(_I, &m, &n); CHKERRQ(ierr);
//     ierr = MatGetOwnershipRange(_I, &first, &last); CHKERRQ(ierr);
//     ierr = MatGetOwnershipRangeColumn(_I, &first_c, &last_c); CHKERRQ(ierr);
//     dim = _sys->n_vars();
//     new_dim = _variable_names.size();
    
//     if (new_dim < dim){
        
//         PetscPrintf(MPI_COMM_WORLD, "reduce_interpolation new_dim < dim \n");
//         nodes_h = m/dim;
//         nodes_H = n/dim;
        
//         length_r = m*new_dim/dim;
//         _t_start = std::chrono::high_resolution_clock::now();
//         // this should use the number of the variable_name!!!
//         ierr = PetscMalloc1(length_r,&indices_r);CHKERRQ(ierr);
//         for (int i=0;i<nodes_h;i++){
//             for (int j=0;j<new_dim;j++){
//                 indices_r[i*new_dim + j] = first + i*dim + j;
//             }
//         }
        
        
//         _t_start = std::chrono::high_resolution_clock::now();
//         ierr = ISCreateGeneral(PETSC_COMM_WORLD,length_r,indices_r,PETSC_COPY_VALUES,&is_r); CHKERRQ(ierr);
//         _t_end = std::chrono::high_resolution_clock::now();
//         fp_ms = _t_end - _t_start ;
//         PetscPrintf(PETSC_COMM_WORLD, "\t computation of is_r %f ms. \n",fp_ms.count());
        
        
//         length_c = n*new_dim/dim;
//         ierr = PetscMalloc1(length_c,&indices_c);CHKERRQ(ierr);
//         for (int i=0;i<nodes_H;i++){
//             for (int j=0;j<new_dim;j++){
//                 indices_c[i*new_dim +j] = first_c + i*dim + j;
//             }
//         }
        
//         ierr = ISCreateGeneral(PETSC_COMM_WORLD,length_c,indices_c,PETSC_COPY_VALUES,&is_c); CHKERRQ(ierr);
//         _t_end = std::chrono::high_resolution_clock::now();
        
//         _t_end = std::chrono::high_resolution_clock::now();
//         fp_ms = _t_end - _t_start ;
//         PetscPrintf(PETSC_COMM_WORLD, "\t computation of is_c %f ms. \n",fp_ms.count());
        
//         _t_start = std::chrono::high_resolution_clock::now();
        
//         ierr = MatGetSubMatrix(_I,is_r,is_c,MAT_INITIAL_MATRIX,_I_red); CHKERRQ(ierr);
//         _t_end = std::chrono::high_resolution_clock::now();
//         _t_end = std::chrono::high_resolution_clock::now();
//         fp_ms = _t_end - _t_start ;
//         PetscPrintf(PETSC_COMM_WORLD, "\t computation of submatrix %f ms. \n",fp_ms.count());
        
//         ierr = ISDestroy(&is_r); CHKERRQ(ierr);
//         ierr = ISDestroy(&is_c); CHKERRQ(ierr);
//         ierr = PetscFree(indices_r); CHKERRQ(ierr);
//         ierr = PetscFree(indices_c); CHKERRQ(ierr);
        
//     }else{
//         PetscPrintf(MPI_COMM_WORLD, "reduce_interpolation just duplicate \n");
//         ierr = MatDuplicate(_I,MAT_COPY_VALUES, _I_red); CHKERRQ(ierr);
//     }
//     PetscFunctionReturn(0);
// };

// /**
//  * @brief      function to determine which variable has D-BC  at certain node
//  *
//  * In the input file: space stands for new BC_id
//  *                    dash stands for variable_id on given node_set
//  *
//  *      for example:
//  *          dc_boundaries = '1 2'
//  *          dc_variables  = '0-3 0-1-2-3'
//  *     -> variables 0 and 3   have prescribed BC on node_sets with ID 1
//  *     -> variables 0,1,2,3   have prescribed BC on node_sets with ID 2
//  *
//  *     In case u messed up in the input file, all variables are assumed to have D-BC on given node_set
//  *
//  * @param[in]  _variables  string as given in input file
//  * result is stored in _dc_variables_id
//  */

// #undef __FUNCT__
// #define __FUNCT__ "determine_dc_bnd_var_id"
// PetscErrorCode determine_dc_bnd_var_id(const std::vector<std::string> & _variables,
//                                        std::vector<int> _dc_boundary_id,
//                                        libMesh::System * _sys,
//                                        std::vector<std::vector<int> > * _dc_variables_id)
// {
//     PetscFunctionBegin;
//     // automatic fill-in
//     std::vector<int> vec(_sys->n_vars());
//     std::iota(vec.begin(), vec.end(), 0);
    
//     unsigned int i;
//     auto str_tmp = _variables.begin();
//     // going over all BC_ids
//     for(i = 0, str_tmp; str_tmp != _variables.end(); i++, str_tmp++)
//     {
//         std::vector<std::string> tmp = split_string(*str_tmp, '-');
        
//         // check if variable assigned in the input file exists for given simulation
//         bool var_flg = 1;
//         for(auto t = tmp.begin(); t != tmp.end(); ++t)
//         {
//             if(atoi(t->c_str()) >= _sys->n_vars())
//                 var_flg = 0;
//         }
        
//         // in case u havent put anything into input file, or u put too much
//         if(*str_tmp == "-1" || var_flg == 0)
//         {
//             _dc_variables_id->push_back(vec);
//         }
//         else
//         {
            
//             unsigned int j;
//             std::vector<int > one_BC_id;
//             auto str_in = tmp.begin();
//             for(j = 0, str_in; str_in != tmp.end(); j ++, str_in++)
//             {
//                 one_BC_id.push_back(atoi(str_in->c_str()));
//             }
//             _dc_variables_id->push_back(one_BC_id);
//         }
//     }
    
//     // check if u have same number of BC_ids in both parameters
//     if(_dc_variables_id->size() != _dc_boundary_id.size())
//     {
//         _dc_variables_id->clear();
//         for(auto i = 0; i != _dc_boundary_id.size(); i++)
//         {
//             _dc_variables_id->push_back(vec);
//         }
//     }
    
//     // print out what is considered for zero-ing
//     std::cout<<" ------ BC CONDITIONS  ------ \n";
//     unsigned int t = 0;
//     for(auto i = _dc_variables_id->begin(); i != _dc_variables_id->end();  t++, i++)
//     {
//         std::cout<<"\n BC_id:  "<< _dc_boundary_id[t] << "   var_ids:  ";
//         std::for_each(i->begin(), i->end(), [](int i){ std::cout << i << "  " ; });
//     }
//     PetscFunctionReturn(0);
// };

// #undef __FUNCT__
// #define __FUNCT__ "apply_dc_bnd_to_interpolation"
// PetscErrorCode apply_dc_bnd_to_interpolation(Mat _I, const std::string _variables , const std::vector<int> _dc_boundary_id, libMesh::System * _sys)
// {
//     PetscErrorCode ierr;
//     std::vector<int> zero_rows;
//     std::vector<std::string> _variable;
//     std::vector<std::vector<int> > _dc_variables_id;
    
//     PetscFunctionBegin;
    
//     libMesh::MeshBase * _mesh = &_sys->get_mesh();
//     unsigned int i = 0;
    
//     libMesh::MeshBase::const_node_iterator it = _mesh->nodes_begin();
//     libMesh::MeshBase::const_node_iterator it_end = _mesh->nodes_end();
    
//     for(auto boundary = _dc_boundary_id.begin(); boundary != _dc_boundary_id.end(); ++boundary, i++)
//     {
//         // iterate just over boundary nodes
//         for (; it != it_end; ++it)
//         {
//             libMesh::Node * _node =  *it;
//             // check if node is in active boundary list
//             if(_mesh->boundary_info->has_boundary_id(_node, *boundary)){
//                 // loop over all variables at this node
//                 for(auto v = 0; v < _variable.size(); v++)
//                 {
//                     const  libMesh::Variable & var  = _sys->variable(v);
//                     unsigned int var_num  = var.number();
                    
//                     _node->n_dofs(0,0);
//                     // see if this variable has any dofs at this node
//                     if(_node->n_dofs(_sys->number(), var_num) > 0)
                        
//                     {
//                         // check if given variable has BC on node
//                         if(std::find(_dc_variables_id[i].begin(), _dc_variables_id[i].end(), var_num) != _dc_variables_id[i].end())
//                         {
//                             // different components are not supported by moose at the moment...
//                             zero_rows.push_back(_node->dof_number(_sys->number(), var_num, 0));
//                         }
//                     }
//                 }
//             }
//         }
//     }
    
//     ierr = MatZeroRows(_I, zero_rows.size(), &zero_rows[0], 0, NULL, NULL); CHKERRQ(ierr);
//     zero_rows.clear();
//     return 0;
//     PetscFunctionReturn(0);
// }

// /*
 
//  Notes:
//  - If  _normals is a null-pointer this should return an identy matrix
//  - Mat _template: the Orthogonal matrix will be structurally set up as a duplicate of this matrix
 
//  */
// #undef __FUNCT__
// #define __FUNCT__ "householder_block_transformation"
// extern PetscErrorCode householder_block_transformation(Vec _normals, int _blocksize, Mat _template, Mat * _Orth){
    
//     PetscErrorCode ierr;
//     int m;
//     int first, last;
    
//     PetscFunctionBegin;
//     if (_blocksize==3){
        
//         PetscScalar * normal_values = nullptr;
//         PetscScalar * values_blk = nullptr;
//         PetscInt * indices_blk = nullptr;
        
//         ierr = PetscMalloc1(_blocksize,&indices_blk);CHKERRQ(ierr);
//         ierr = MatDuplicate(_template, MAT_DO_NOT_COPY_VALUES, _Orth); CHKERRQ(ierr);
        
//         ierr = MatShift(*_Orth, 1);CHKERRQ(ierr);
//         ierr = MatGetOwnershipRange(*_Orth,&first,&last); CHKERRQ(ierr);
//         ierr = VecGetLocalSize(_normals, &m);CHKERRQ(ierr);
//         ierr = VecGetArray(_normals, &normal_values);CHKERRQ(ierr);
        
//         Eigen::Map<Eigen::VectorXd> blk_normal(&normal_values[0], _blocksize);
//         Eigen::VectorXd blk_e1(_blocksize);
//         blk_e1 << 1,0,0;
//         Eigen::VectorXd blk_v(_blocksize);
//         Eigen::MatrixXd blk_hh(_blocksize,_blocksize);
//         Eigen::MatrixXd blk_i(_blocksize,_blocksize);
        
//         blk_i.setIdentity();
//         if (_normals != NULL){
//             for (int i=0;i<m/_blocksize;i++){
//                 new (&blk_normal) Eigen::Map<Eigen::VectorXd> (&normal_values[i*_blocksize], _blocksize,1);
                
//                 if (blk_normal.norm() > std::numeric_limits<double>::epsilon()) // hack: we should loop over nodes of \Gamma_C
//                 {
//                     for (int j=0;j<_blocksize;j++)
//                     {
//                         indices_blk[j] = first+i*_blocksize+j;
//                     }
//                     blk_v = (- blk_e1 -blk_normal);
//                     blk_v *= 0.5;
//                     blk_hh = blk_i - 2*(blk_v*blk_v.transpose())/(blk_v.transpose()*blk_v);
                    
//                     for (int j=0;j<_blocksize;j++)
//                     {
//                         indices_blk[j] = first+i*_blocksize+j;
//                     }
                    
//                     values_blk = blk_hh.data();
                    
//                     ierr = MatSetValues(*_Orth,_blocksize,indices_blk,_blocksize,indices_blk,values_blk,INSERT_VALUES); CHKERRQ(ierr);
//                 }
//             }
//         }
//         ierr = VecRestoreArray(_normals, &normal_values); CHKERRQ(ierr);
//         ierr = MatAssemblyBegin(*_Orth, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//         ierr = MatAssemblyEnd(*_Orth,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        
//         ierr = PetscFree(indices_blk);CHKERRQ(ierr);
//     }
//     PetscFunctionReturn(0);
// };
// #undef __FUNCT__
// #define __FUNCT__ "remove_zero_rows"
// extern PetscErrorCode remove_zero_rows(Mat * M){
    
//     PetscErrorCode ierr;
//     Mat tmp;
//     IS nonzero_rows;
    
//     PetscFunctionBegin;
    
//     ierr = MatFindNonzeroRows(*M,&nonzero_rows); CHKERRQ(ierr);
//     if  (nonzero_rows != NULL){
//         ierr = MatGetSubMatrix(*M,nonzero_rows,NULL,MAT_INITIAL_MATRIX,&tmp); CHKERRQ(ierr);
//         ierr = MatDestroy(M); CHKERRQ(ierr);
//         M = &tmp;
//     }
    
//     PetscFunctionReturn(0);
// };

// #undef __FUNCT__
// #define __FUNCT__ "remove_singularities"
// extern PetscErrorCode remove_singularities(Mat * matrix){
    
//     PetscErrorCode ierr;
    
//     Vec diag;
    
//     int M,N,m,n, first, last;
//     PetscScalar * vals;
//     PetscInt *idx;
//     std::vector<int> zero_rows;
    
//     std::chrono::high_resolution_clock::time_point _t_start;
//     std::chrono::high_resolution_clock::time_point _t_end;
//     std::chrono::duration<double, std::milli> fp_ms;
    
//     PetscFunctionBegin;
    
//     ierr = MatGetSize(*matrix, &M, &N); CHKERRQ(ierr);
//     ierr = MatGetLocalSize(*matrix, &m, &n); CHKERRQ(ierr);
    
//     ierr = VecCreate(PETSC_COMM_WORLD, &diag); CHKERRQ(ierr);
//     ierr = VecSetSizes(diag, n, PETSC_DECIDE); CHKERRQ(ierr);
//     ierr = VecSetFromOptions(diag); CHKERRQ(ierr);
//     ierr = VecAssemblyBegin(diag); CHKERRQ(ierr);
//     ierr = VecAssemblyEnd(diag); CHKERRQ(ierr);
    
//     _t_start = std::chrono::high_resolution_clock::now();
//     ierr = MatGetRowMaxAbs(*matrix,diag,NULL); CHKERRQ(ierr);
//     _t_end = std::chrono::high_resolution_clock::now();
//     fp_ms = _t_end - _t_start ;
//     PetscPrintf(PETSC_COMM_WORLD, "\t MatGetRowMaxAbs %f ms. \n",fp_ms.count());
    
//     _t_start = std::chrono::high_resolution_clock::now();
//     ierr = VecGetOwnershipRange(diag, &first, &last); CHKERRQ(ierr);
//     ierr = VecGetArray(diag, &vals); CHKERRQ(ierr);
//     PetscPrintf(PETSC_COMM_WORLD, "setting rows ");
//     for (int i=0;i<n;i++){
//         if (fabs(vals[i]) < PASSO_TOL){
//             // set zero row & col
//             zero_rows.push_back((i+first));
//             PetscPrintf(PETSC_COMM_WORLD, "%i ", i+first);
//             vals[i] = 1;
//         }else{
//             vals[i] = 0;
//         };
//     }
//     PetscPrintf(PETSC_COMM_WORLD, " to one. \n");
//     _t_end = std::chrono::high_resolution_clock::now();
//     fp_ms = _t_end - _t_start ;
//     PetscPrintf(PETSC_COMM_WORLD, "\t finding zero rows: %f ms. \n",fp_ms.count());
    
    
//     ierr = VecRestoreArray(diag, &vals); CHKERRQ(ierr);
//     ierr = VecAssemblyBegin(diag); CHKERRQ(ierr);
//     ierr = VecAssemblyEnd(diag); CHKERRQ(ierr);
    
//     ierr = MatZeroRows(*matrix, zero_rows.size(), &zero_rows[0], 0, NULL, NULL); CHKERRQ(ierr);
//     ierr = MatTranspose(*matrix, MAT_REUSE_MATRIX, matrix);
//     ierr = MatZeroRows(*matrix, zero_rows.size(), &zero_rows[0], 0, NULL, NULL); CHKERRQ(ierr);
//     ierr = MatTranspose(*matrix, MAT_REUSE_MATRIX, matrix);
    
//     PetscPrintf(PETSC_COMM_WORLD, "Start adding diagonal ones ... ");
//     _t_start = std::chrono::high_resolution_clock::now();
//     ierr = MatSetOption(*matrix,MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);CHKERRQ(ierr);
//     ierr = MatDiagonalSet(*matrix,diag,ADD_VALUES); CHKERRQ(ierr);
//     _t_end = std::chrono::high_resolution_clock::now();
//     fp_ms = _t_end - _t_start ;
//     PetscPrintf(PETSC_COMM_WORLD, "done %f ms. \n",fp_ms.count());
    
//     ierr = VecDestroy(&diag);CHKERRQ(ierr);
    
//     PetscFunctionReturn(0);
// };





// #undef __FUNCT__
// #define __FUNCT__ "remove_zero_cols"
// extern PetscErrorCode remove_zero_cols(Mat * M){
    
//     PetscErrorCode ierr;
//     Mat tmp, tmp2;
//     IS nonzero_rows;
    
//     PetscFunctionBegin;
    
//     //    ierr = ISCreate(PETSC_COMM_WORLD, &nonzero_rows); CHKERRQ(ierr);
//     ierr = MatTranspose(*M, MAT_INITIAL_MATRIX,&tmp); CHKERRQ(ierr);
//     ierr = MatDestroy(M);CHKERRQ(ierr);
//     ierr = MatFindNonzeroRows(tmp,&nonzero_rows); CHKERRQ(ierr);
//     if  (nonzero_rows != NULL){
//         ierr = MatGetSubMatrix(tmp,nonzero_rows,NULL,MAT_INITIAL_MATRIX,&tmp2); CHKERRQ(ierr);
//         ierr = MatTranspose(tmp2, MAT_INITIAL_MATRIX,M); CHKERRQ(ierr);
        
//         ierr = MatDestroy(&tmp2); CHKERRQ(ierr);
//     }else{
//         ierr = MatTranspose(tmp, MAT_INITIAL_MATRIX,M); CHKERRQ(ierr);
//     }
    
//     ierr = MatDestroy(&tmp); CHKERRQ(ierr);
    
//     PetscFunctionReturn(0);
// };

// /*  computes enorm  =  <c,c>_A - <f,c> */
// #undef __FUNCT__
// #define __FUNCT__ "energy_norm"
// extern PetscErrorCode energy_norm(Mat A, Vec f, Vec c, PetscScalar * enorm){
//     PetscFunctionBegin;
//     PetscErrorCode ierr;
//     PetscScalar ctAc=0, ftc=0, norm_c=0, norm_f = 0;

//     xtAx(A, c, &ctAc);
//     ierr = VecNorm(c, NORM_2, &norm_c); CHKERRQ(ierr);
//     ierr = VecNorm(f, NORM_2, &norm_f); CHKERRQ(ierr);
//     ierr = VecDot(f,c, &ftc); CHKERRQ(ierr);
//     *enorm = (0.5*ctAc - ftc);
    
//     //    PetscPrintf(PETSC_COMM_WORLD, "energy_norm(...): ctAc = %f, ftc=%f, |c|= %f, energy= %f \n ", ctAc, ftc, norm_c, *enorm);
//     PetscFunctionReturn(0);
// };

// /*  residual of the inactive set (normally nodes not in contact) */
// #undef __FUNCT__
// #define __FUNCT__ "inactive_set_residual_norm"
// extern PetscErrorCode inactive_set_residual_norm(Mat A, Vec f, Vec c, Vec active_set, PetscScalar * norm){
//     PetscFunctionBegin;

//     PetscErrorCode ierr;
//     PetscInt n,m;
//     Vec mask, c_tmp, res_tmp;

//     ierr = MatGetLocalSize(A, &n, &m); CHKERRQ(ierr);
//     {   // determine inactive nodes

//         ierr = VecDuplicate(active_set, &mask); CHKERRQ(ierr);
//         ierr = VecScale(mask,0); CHKERRQ(ierr);

//         PetscScalar * active_coord_vals, * mask_vals;
//         ierr = VecGetArray(active_set, &active_coord_vals); CHKERRQ(ierr);
//         ierr = VecGetArray(mask, &mask_vals); CHKERRQ(ierr);
//         for (int i=0;i<n;i++){
//             if (active_coord_vals[i] == 0){
//                 mask_vals[i] = 1;
//             }else{
//                 mask_vals[i] = 1;
//             }
//         }
//         ierr = VecRestoreArray(active_set, &active_coord_vals); CHKERRQ(ierr);
//         ierr = VecRestoreArray(mask, &mask_vals); CHKERRQ(ierr);

//     }

//     ierr = VecDuplicate(c, &c_tmp); CHKERRQ(ierr);
//     ierr = VecDuplicate(c, &res_tmp); CHKERRQ(ierr);
//     ierr = VecPointwiseMult(c_tmp, c, mask); CHKERRQ(ierr);
//     lin_residual(A, f, c_tmp, res_tmp);
//     ierr = VecNorm(res_tmp, NORM_2, norm);

//     ierr = VecDestroy(&mask);
//     ierr = VecDestroy(&c_tmp);
//     ierr = VecDestroy(&res_tmp);

//     PetscFunctionReturn(0);
// };


// /*  computes alpha = <g,c>/<c,c>_A */
// #undef __FUNCT__
// #define __FUNCT__ "passo_linesearch"
// extern PetscErrorCode passo_linesearch(Vec g, Mat A, Vec c, PetscScalar * alpha){
//     PetscFunctionBegin;
    
//     PetscErrorCode ierr;
//     PetscScalar nom, denom =1;
    
//     xtAx(A, c, &denom);
//     ierr = VecDot(g, c, &nom); CHKERRQ(ierr);
//     *alpha = nom/denom;
    
//     PetscFunctionReturn(0);
// }

/*@C
 SNESConvergedPasso - Convergence test of the solvers for
 systems of nonlinear equations (copied from default).
 
 Collective on SNES
 
 Input Parameters:
 +  snes - the SNES context
 .  it - the iteration (0 indicates before any Newton steps)
 .  xnorm - 2-norm of current iterate
 .  snorm - 2-norm of current step
 .  fnorm - 2-norm of function at current iterate
 -  dummy - unused context
 
 Output Parameter:
 .   reason  - one of
 $  SNES_CONVERGED_FNORM_ABS       - (fnorm < abstol),
 $  SNES_CONVERGED_SNORM_RELATIVE  - (snorm < stol*xnorm),
 $  SNES_CONVERGED_FNORM_RELATIVE  - (fnorm < rtol*fnorm0),
 $  SNES_DIVERGED_FUNCTION_COUNT   - (nfct > maxf),
 $  SNES_DIVERGED_FNORM_NAN        - (fnorm == NaN),
 $  SNES_CONVERGED_ITERATING       - (otherwise),
 
 where
 +    maxf - maximum number of function evaluations,
 set with SNESSetTolerances()
 .    nfct - number of function evaluations,
 .    abstol - absolute function norm tolerance,
 set with SNESSetTolerances()
 -    rtol - relative function norm tolerance, set with SNESSetTolerances()
 
 Level: intermediate
 
 .keywords: SNES, nonlinear, default, converged, convergence
 
 .seealso: SNESSetConvergenceTest()
 @*/
#undef __FUNCT__
#define __FUNCT__ "SNESConvergenceMooseUtopia"
PetscErrorCode
SNESConvergenceMooseUtopia(SNES snes,PetscInt it,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason,void * /*dummy */)
{
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
    PetscValidPointer(reason,6);
    
    *reason = SNES_CONVERGED_ITERATING;
    //    PetscPrintf(PETSC_COMM_SELF, "Check convergence, rtol = %2.12f \n", snes->rtol);
    
    if (!it) {
        /* set parameter for default relative tolerance convergence test */
        //        snes->ttol = fnorm*snes->rtol;
        snes->ttol = snes->rtol;
    }
    if (PetscIsInfOrNanReal(fnorm)) {
        ierr    = PetscPrintf(PETSC_COMM_SELF, "Failed to converged, function norm is NaN\n");CHKERRQ(ierr);
        *reason = SNES_DIVERGED_FNORM_NAN;
    } else if (fnorm < snes->abstol) {
        ierr    = PetscPrintf(PETSC_COMM_SELF, "Converged due to function norm %4.12e < %4.12e\n",(double)fnorm,(double)snes->abstol);CHKERRQ(ierr);
        *reason = SNES_CONVERGED_FNORM_ABS;
    } else if (snes->nfuncs >= snes->max_funcs) {
        ierr    =PetscPrintf(PETSC_COMM_SELF, "Exceeded maximum number of function evaluations: %D > %D\n",snes->nfuncs,snes->max_funcs);CHKERRQ(ierr);
        *reason = SNES_DIVERGED_FUNCTION_COUNT;
    }
    
    if (it && !*reason) {
        if (fnorm <= snes->ttol) {
            ierr    = PetscPrintf(PETSC_COMM_SELF,"Converged due to function norm %4.12e < %4.12e \n",(double)fnorm,(double)snes->ttol);CHKERRQ(ierr);
            *reason = SNES_CONVERGED_FNORM_RELATIVE;
        } else if (snorm < snes->stol*xnorm) {
            ierr    =  PetscPrintf(PETSC_COMM_SELF, "Converged due to small update length: %4.12e < %4.12e * %4.12e\n",(double)snorm,(double)snes->stol,(double)xnorm);CHKERRQ(ierr);
            *reason = SNES_CONVERGED_SNORM_RELATIVE;
        }
    }
    PetscFunctionReturn(0);
}
#undef __FUNCT__


#define __FUNCT__ "SNESSetFromOptions_MooseUtopia_options"
PetscErrorCode
SNESSetFromOptions_MooseUtopia_options(PetscOptionItems * PetscOptionsObject, SNES snes)
{
    SNES_MOOSEUTOPIA * ctx = (SNES_MOOSEUTOPIA *)snes->data;
    
    PetscErrorCode ierr;
    PetscFunctionBegin;
    ierr = PetscOptionsHead(PetscOptionsObject, "SNES MOOSEUTOPIA PARAMETERS options");
    CHKERRQ(ierr);
    
    PetscOptionsBool(
                     "-MooseUtopia_verbose", "alg", "None", ctx->MooseUtopia_params.verbose, &ctx->MooseUtopia_params.verbose, NULL);
    PetscOptionsBool("-MooseUtopia_time_statistics",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.time_statistics,
                     &ctx->MooseUtopia_params.time_statistics,
                     NULL);
    PetscOptionsInt(
                    "-MooseUtopia_max_it", "alg", "None", ctx->MooseUtopia_params.max_it, &ctx->MooseUtopia_params.max_it, NULL);
    
    PetscOptionsReal(
                     "-MooseUtopia_atol", "alg", "None", ctx->MooseUtopia_params.abs_tol, &ctx->MooseUtopia_params.abs_tol, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_stol", "alg", "None", ctx->MooseUtopia_params.stol, &ctx->MooseUtopia_params.stol, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_rtol", "alg", "None", ctx->MooseUtopia_params.rtol, &ctx->MooseUtopia_params.rtol, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_tol", "alg", "None", ctx->MooseUtopia_params.tol, &ctx->MooseUtopia_params.tol, NULL);
    
    //------------------------------------------- TR -------------------------------------------//
    PetscStrcpy(ctx->MooseUtopia_params.tr_alg, "DOGLEG");
    PetscOptionsString(
                       "-MooseUtopia_tr_alg", "", "", ctx->MooseUtopia_params.tr_alg, ctx->MooseUtopia_params.tr_alg, 256, NULL);
    PetscStrcpy(ctx->MooseUtopia_params.ls_alg, "BACKTRACKING");
    PetscOptionsString(
                       "-MooseUtopia_ls_alg", "", "", ctx->MooseUtopia_params.ls_alg, ctx->MooseUtopia_params.ls_alg, 256, NULL);
    
    PetscOptionsReal("-MooseUtopia_delta_max",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.delta_max,
                     &ctx->MooseUtopia_params.delta_max,
                     NULL);
    PetscOptionsReal(
                     "-MooseUtopia_delta0", "alg", "None", ctx->MooseUtopia_params.delta0, &ctx->MooseUtopia_params.delta0, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_gamma1", "alg", "None", ctx->MooseUtopia_params.gamma1, &ctx->MooseUtopia_params.gamma1, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_gamma2", "alg", "None", ctx->MooseUtopia_params.gamma1, &ctx->MooseUtopia_params.gamma2, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_eta1", "alg", "None", ctx->MooseUtopia_params.eta1, &ctx->MooseUtopia_params.eta1, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_eta2", "alg", "None", ctx->MooseUtopia_params.eta2, &ctx->MooseUtopia_params.eta2, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_rho_tol", "alg", "None", ctx->MooseUtopia_params.rho_tol, &ctx->MooseUtopia_params.rho_tol, NULL);
    PetscOptionsReal(
                     "-MooseUtopia_ST_tol", "alg", "None", ctx->MooseUtopia_params.ST_tol, &ctx->MooseUtopia_params.ST_tol, NULL);
    
    PetscOptionsReal(
                     "-MooseUtopia_nl_atol", "alg", "None", ctx->MooseUtopia_params.nl_atol, &ctx->MooseUtopia_params.nl_atol, NULL);
    PetscOptionsInt("-MooseUtopia_local_max_it",
                    "alg",
                    "None",
                    ctx->MooseUtopia_params.local_max_it,
                    &ctx->MooseUtopia_params.local_max_it,
                    NULL);
    PetscOptionsInt(
                    "-MooseUtopia_overlap", "alg", "None", ctx->MooseUtopia_params.overlap, &ctx->MooseUtopia_params.overlap, NULL);
    
    //------------------------------------------- MG -------------------------------------------//
    PetscOptionsInt("-MooseUtopia_presmoothing",
                    "alg",
                    "None",
                    ctx->MooseUtopia_params.mg_presmoothing_steps,
                    &ctx->MooseUtopia_params.mg_presmoothing_steps,
                    NULL);
    PetscOptionsInt("-MooseUtopia_postsmoothing",
                    "alg",
                    "None",
                    ctx->MooseUtopia_params.mg_postsmoothing_steps,
                    &ctx->MooseUtopia_params.mg_postsmoothing_steps,
                    NULL);
    PetscOptionsInt(
                    "-MooseUtopia_MG_type", "alg", "None", ctx->MooseUtopia_params.mg_type, &ctx->MooseUtopia_params.mg_type, NULL);
    PetscOptionsInt("-MooseUtopia_blocksize", "alg", "None", ctx->MooseUtopia_params.blocksize, &ctx->MooseUtopia_params.blocksize, NULL);
    PetscOptionsReal("-MooseUtopia_omega", "alg", "None", ctx->MooseUtopia_params.omega, &ctx->MooseUtopia_params.omega, NULL);
    
    PetscStrcpy(ctx->MooseUtopia_params.smoother_type, "jacobi");
    PetscOptionsString("-MooseUtopia_smoother_type",
                       "",
                       "",
                       ctx->MooseUtopia_params.smoother_type,
                       ctx->MooseUtopia_params.smoother_type,
                       256,
                       NULL);
    
    PetscStrcpy(ctx->MooseUtopia_params.coarse_solve_type, "jacobi");
    PetscOptionsString("-MooseUtopia_coarse_solve_type",
                       "",
                       "",
                       ctx->MooseUtopia_params.coarse_solve_type,
                       ctx->MooseUtopia_params.coarse_solve_type,
                       256,
                       NULL);
    
    PetscOptionsBool("-MooseUtopia_static_time_step",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.static_time_step,
                     &ctx->MooseUtopia_params.static_time_step,
                     NULL);
    
    PetscOptionsBool("-MooseUtopia_mg_native",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.mg_native,
                     &ctx->MooseUtopia_params.mg_native,
                     NULL);
    
    //------------------------------------------- LS -------------------------------------------//
    PetscOptionsReal(
                     "-MooseUtopia_ls_rho", "alg", "None", ctx->MooseUtopia_params.ls_rho, &ctx->MooseUtopia_params.ls_rho, NULL);
    PetscOptionsReal("-MooseUtopia_c1", "alg", "None", ctx->MooseUtopia_params.c1, &ctx->MooseUtopia_params.c1, NULL);
    PetscOptionsReal("-MooseUtopia_c2", "alg", "None", ctx->MooseUtopia_params.c2, &ctx->MooseUtopia_params.c2, NULL);
    
    //------------------------------------------- Linear solver
    //-------------------------------------------//
    PetscStrcpy(ctx->MooseUtopia_params.ksp_type, "preonly");
    PetscOptionsString(
                       "-MooseUtopia_ksp_type", "", "", ctx->MooseUtopia_params.ksp_type, ctx->MooseUtopia_params.ksp_type, 256, NULL);
    PetscStrcpy(ctx->MooseUtopia_params.pc_type, "lu");
    PetscOptionsString(
                       "-MooseUtopia_pc_type", "", "", ctx->MooseUtopia_params.pc_type, ctx->MooseUtopia_params.pc_type, 256, NULL);
    PetscStrcpy(ctx->MooseUtopia_params.pc_factor_mat_solver_package, "mumps");
    PetscOptionsString("-MooseUtopia_pc_factor_mat_solver_package",
                       "",
                       "",
                       ctx->MooseUtopia_params.pc_factor_mat_solver_package,
                       ctx->MooseUtopia_params.pc_factor_mat_solver_package,
                       256,
                       NULL);
    PetscOptionsReal("-MooseUtopia_ksp_atol",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.ksp_atol,
                     &ctx->MooseUtopia_params.ksp_atol,
                     NULL);
    PetscOptionsReal("-MooseUtopia_ksp_rtol",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.ksp_rtol,
                     &ctx->MooseUtopia_params.ksp_rtol,
                     NULL);
    PetscOptionsReal("-MooseUtopia_ksp_dtol",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.ksp_dtol,
                     &ctx->MooseUtopia_params.ksp_dtol,
                     NULL);
    PetscOptionsReal("-MooseUtopia_ksp_max_it",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.ksp_max_it,
                     &ctx->MooseUtopia_params.ksp_max_it,
                     NULL);
    PetscOptionsBool("-MooseUtopia_linear_solver_verbose",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.linear_solver_verbose,
                     &ctx->MooseUtopia_params.linear_solver_verbose,
                     NULL);
    
    // ----------------------------------------- Hessian update
    // -----------------------------------------------//
    PetscOptionsInt("-MooseUtopia_hessian_update",
                    "alg",
                    "None",
                    ctx->MooseUtopia_params.hessian_update,
                    &ctx->MooseUtopia_params.hessian_update,
                    NULL);
    
    //------------------------------------------- general logging and debug options
    //------------------------------------//
    PetscOptionsBool("-MooseUtopia_log_iterates",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.MooseUtopia_log_iterates,
                     &ctx->MooseUtopia_params.MooseUtopia_log_iterates,
                     NULL);
    PetscOptionsBool("-MooseUtopia_log_system",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.MosseUtopia_log_system,
                     &ctx->MooseUtopia_params.MosseUtopia_log_system,
                     NULL);
    PetscOptionsBool("-MooseUtopia_log_norms",
                     "alg",
                     "None",
                     ctx->MooseUtopia_params.MooseUtopia_log_norms,
                     &ctx->MooseUtopia_params.MooseUtopia_log_norms,
                     NULL);
    
    PetscOptionsString("-MooseUtopia_path_to_log", "", "", ctx->MooseUtopia_params.path_to_log, ctx->MooseUtopia_params.path_to_log, 256, NULL);

    PetscOptionsString("-MooseUtopia_consistency_level", "", "", ctx->MooseUtopia_params.consistency_level, ctx->MooseUtopia_params.consistency_level, 256, NULL);


    ierr = PetscOptionsTail();
    CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}
#undef __FUNCT__

