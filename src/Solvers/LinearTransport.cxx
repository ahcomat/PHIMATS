#include "Solvers/LinearTransport.h"
#include<iostream>

LinearTransport::LinearTransport(Mat &A, Logger& logger, string solverType){

    // Initialize the solver.
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);

    try{
        if (solverType=="DIRECT"){

            // Direct solver
            KSPSetType(ksp, KSPPREONLY); // Direct solver for linear systems
            KSPSetFromOptions(ksp);

            // Set preconditioner
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCLU); // LU preconditioner for direct solver
            PCSetFromOptions(pc);

            logger.log("    Using < "+solverType+" > solver" , "INFO");
            logger.log("", "", false);

        } else if (solverType=="GMRES"){

            // GMRES solver
            KSPSetType(ksp, KSPGMRES);
            KSPGMRESSetRestart(ksp, 50); // Optional: Set GMRES restart value
            KSPSetTolerances(ksp, 1e-12, 1e-12, PETSC_DEFAULT, 1000);
            KSPSetFromOptions(ksp);

            // Set preconditioner
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCILU);
            PCFactorSetLevels(pc, 2); // ILU with limited fill
            PCSetFromOptions(pc);

            logger.log("    Using solver < "+solverType+" >" , "INFO");
            logger.log("", "", false);

        } else {

            logger.log("Undefined solver type < " + solverType + " >", "ERROR");
            throw std::runtime_error("Undefined solver type < " + solverType + " >. Supported options are < DIRECT, GMRES >");

        }
    } catch (const std::runtime_error& e) {
        cerr << "\nException caught:: " << e.what() << endl;
        cerr << "Terminating!\n" << endl;
        exit(EXIT_FAILURE);
    }
}

LinearTransport::~LinearTransport(){

    // Exit message
    std::cout << "LinearTransport solver exited correctly" << "\n";
}

void LinearTransport::UpdateKSP(Mat &A){

    KSPSetOperators(ksp, A, A);
}

void LinearTransport::SolveTransport(Vec &x, Vec &F){

    KSPSolve(ksp, F, x);
}