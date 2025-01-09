#include "Solvers/LinearMech.h"
#include<iostream>

LinearMech::LinearMech(Mat &A){

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

LinearMech::~LinearMech(){

    // Exit message
    std::cout << "LinearMech solver exited correctly" << "\n";
}

void LinearMech::Solve(Vec &x, Vec &b){

    KSPSolve(ksp, b, x);
}
