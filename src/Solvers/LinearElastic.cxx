/**
 * @file LinearElastic.cxx
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A wrapper for PETSc linear solver `ksp` mainly for linear elasticity.
 * @date 2024-05-28
 * 
 * @copyright Copyright (c) 2024
 * 
 * Updates (when, what and who)
 * 
 */

#include "Solvers/LinearElastic.h"
#include<iostream>

LinearElastic::LinearElastic(Mat &A){

    // Initialize the solver. Default GMRES. We could use direct solver.
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetType(ksp, KSPPREONLY); // Direct solver for linear systems
    KSPSetFromOptions(ksp);

    // Set preconditioner
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU); // LU preconditioner for direct solver
    PCSetFromOptions(pc);
}

LinearElastic::~LinearElastic(){

    // Exit message
    std::cout << "LinearElastic solver exited correctly" << "\n";
}

void LinearElastic::Solve(Vec &x, Vec &b){

    KSPSolve(ksp, b, x);
    VecCopy(x, b);
}
