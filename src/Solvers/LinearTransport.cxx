#include "Solvers/LinearTransport.h"
#include<iostream>

LinearTransport::LinearTransport(Mat &A){

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

LinearTransport::~LinearTransport(){

    // Exit message
    std::cout << "LinearTransport solver exited correctly" << "\n";
}

void LinearTransport::UpdateKSP(Mat &A){

    KSPSetOperators(ksp, A, A);
    KSPSetFromOptions(ksp);

    // Set preconditioner
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU); // LU preconditioner for direct solver
    PCSetFromOptions(pc);
}

void LinearTransport::SolveTransport(Vec &x, Vec &F){

    KSPSolve(ksp, F, x);
}