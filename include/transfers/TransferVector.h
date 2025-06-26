// #ifndef VectorTransfer_H
// #define VectorTransfer_H

// #include "MultiAppTransfer.h"
// #include <utopia_fe.hpp>
// #include <utopia_LibMeshBackend.hpp>
// #include "utopia.hpp"
// #include <memory>
// #include "UserObjectInterface.h"
// #include "StoreTransferOperators.h"
// class VectorTransfer;


// template <>
// InputParameters validParams<VectorTransfer>();

// /* It projects values from a master mesh to a slave mesh */

// class VectorTransfer : public MultiAppTransfer,
// public UserObjectInterface
// {
// public: VectorTransfer(const InputParameters & parameters);

// virtual void initialSetup();
  
// virtual void execute();

// typedef utopia::DSMatrixd SparseMatT;

// typedef utopia::DVectord VecT;

// protected:

// AuxVariableName _to_var_name;

// VariableName _from_var_name;

// bool _impact;

// //userobject to store our operators
// const StoreTransferOperators & _operator_storage;

// std::string _operator_type;

// void buildTransferOperators();
  
// virtual void get_solution_vec(FEProblemBase * problem, std::string var_name,Vec &vector);
    
// virtual void apply_transpose_operator();
   
// virtual void apply_operator();
    
// int _n_var;
    
// void  Permutation(FEProblemBase *problem, std::shared_ptr<SparseMatT> P, std::string main_var,
//                   std::string aux_var, unsigned int number);

// void copy_sol(VecT sol_to_copy);
    
// std::shared_ptr<SparseMatT> _P= NULL;
    
// std::shared_ptr<SparseMatT> _P_reverse= NULL;

// unsigned int _var_1_from, _var_1_to;
// unsigned int _var_2_from, _var_2_to;
// unsigned int _var_3_from, _var_3_to;
    
    

// };
// //
// #endif /* VectorTrasfer */
