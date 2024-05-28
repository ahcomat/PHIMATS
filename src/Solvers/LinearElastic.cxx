#include "Solvers/LinearElastic.h"
#include<iostream>

LinearElastic::LinearElastic(Mat &A){

    // Initialize the solver. Default GMRES. We could use direct solver.
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetType(ksp, KSPPREONLY); // Direct solver for linear systems
    KSPSetFromOptions(ksp);

    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU); // LU preconditioner for direct solver
    PCSetFromOptions(pc);
}

LinearElastic::~LinearElastic(){

    KSPDestroy(&ksp);

    // Exit message
    std::cout << "LinearElastic solver exited correctly" << "\n";

}

void LinearElastic::Solve(Vec &x, Vec &b){

    // KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);

    KSPSolve(ksp, b, x);
    VecCopy(x, b);
}
