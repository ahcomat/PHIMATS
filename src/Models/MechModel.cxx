#include<iostream>
#include<variant>
#include <unordered_set>

#include "Models/MechModel.h"
#include "Materials/Mechanics/IsoHard.h"

using namespace std;

MechModel::MechModel(vector<BaseElemMech*> elements, H5IO& H5File_in, Logger& logger, const int NR_update)
    : logger(logger), NR_freq(NR_update) {

    string dsetName;
    dsetName = "SimulationParameters/nSteps";
    nSteps = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/nDim";
    nDim = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/nTotElements";
    nTotElements = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/nElementSets";
    nElementSets = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/nTotDofs";
    nTotDofs = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/nTotNodes";
    nTotNodes = H5File_in.ReadScalar(dsetName);

    // Set the type and size
    nodCount.resize(nTotNodes);
    nodStran_eq.resize(nTotNodes);
    nodStres_eq.resize(nTotNodes);
    nodStres_h.resize(nTotNodes);
    nodRho.resize(nTotNodes);
    if (nDim == 2){ // Case 2D model

        nodStres = vector<ColVecd3>(nTotNodes);
        nodStran = vector<ColVecd3>(nTotNodes);
        nodStran_e = vector<ColVecd3>(nTotNodes);
        nodStran_p = vector<ColVecd3>(nTotNodes);

    } else if (nDim == 3){ // Case 3D

        nodStres = vector<ColVecd6>(nTotNodes);
        nodStran = vector<ColVecd6>(nTotNodes);
        nodStran_e = vector<ColVecd6>(nTotNodes);
        nodStran_p = vector<ColVecd6>(nTotNodes);

    }

    // Allocate memory for `Fint` and `indices`.
    PetscMalloc1(nTotDofs, &Fint);
    PetscMalloc1(nTotDofs, &indices);

    // Initialize to zero.
    for (int iDof=0; iDof<nTotDofs; iDof++){
        Fint[iDof] = 0;
    }

    // Initialize to zeros, otherwise will get garbage memory values.
    setZeroNodVals();

    // Set indices
    for (int iDof=0; iDof<nTotDofs; iDof++){
        indices[iDof] = iDof;
    }

    InitializePETSc(elements);

}

MechModel::~MechModel(){

    // Deallocate memory.
    PetscFree(presDofs); PetscFree(presVals); PetscFree(Fint); PetscFree(indices);
    VecDestroy(&vecFext); VecDestroy(&vecDisp); VecDestroy(&vecFint); 
    VecDestroy(&vecR); MatDestroy(&matA);
    SNESDestroy(&snes);
    // Finalize PETSc
    PetscFinalize();
    // Exit message
    cout << "MechModel elements exited correctly" << "\n";
}

void MechModel::setZeroNodVals(){

    if (nDim == 2){ // Case 2D model

        for(int iNod=0; iNod<nTotNodes; iNod++){

            std::get<std::vector<ColVecd3>>(nodStran).at(iNod).setZero();
            std::get<std::vector<ColVecd3>>(nodStran_e).at(iNod).setZero();
            std::get<std::vector<ColVecd3>>(nodStran_p).at(iNod).setZero();
            std::get<std::vector<ColVecd3>>(nodStres).at(iNod).setZero();
            nodStran_eq.at(iNod) = 0;
            nodStres_eq.at(iNod) = 0;
            nodStres_h.at(iNod) = 0;
            nodRho.at(iNod) = 0;
            nodCount.at(iNod) = 0;
            
        }

    } else if (nDim == 3){ // Case 3D

        for(int iNod=0; iNod<nTotNodes; iNod++){
            
            std::get<std::vector<ColVecd6>>(nodStran).at(iNod).setZero();
            std::get<std::vector<ColVecd6>>(nodStran_e).at(iNod).setZero();
            std::get<std::vector<ColVecd6>>(nodStran_p).at(iNod).setZero();
            std::get<std::vector<ColVecd6>>(nodStres).at(iNod).setZero();
            nodStran_eq.at(iNod) = 0;
            nodStres_eq.at(iNod) = 0;
            nodStres_h.at(iNod) = 0;
            nodRho.at(iNod) = 0;
            nodCount.at(iNod) = 0;

        }
    }
}

void MechModel::InitializeDirichBC(H5IO& H5File_in){

    // Read Dirichlet BCs
    string dsetName;
    dsetName = "SimulationParameters/nPresDofs";
    nPresDofs = H5File_in.ReadScalar(dsetName);
    PetscMalloc1(nPresDofs, &presDofs);
    PetscMalloc1(nPresDofs, &presVals); 
    PetscMalloc1(nPresDofs, &presZeros);

    vector<double> dummy(3);
    for (int iPresDof=0; iPresDof<nPresDofs; iPresDof++){
        // Read values
        dsetName = "PrescribedDOFs/Prescribed_"+to_string(iPresDof);
        H5File_in.ReadField1D(dsetName, dummy);
        // Assign values
        presDofs[iPresDof] = nDim*dummy.at(0)+dummy.at(1); // nDim*iNode+dof
        presVals[iPresDof] = dummy.at(2)/nSteps;
        // // TODO: For debug!
        // cout << presDofs[iPresDof] << " --> " << presVals[iPresDof] << "\n";
        presZeros[iPresDof] = 0;
    }
}

void MechModel::setDirichBC(){

    VecSetValues(vecFext, nPresDofs, presDofs, presVals, INSERT_VALUES); 
    VecAssemblyBegin(vecFext); VecAssemblyEnd(vecFext);
}

int MechModel::get_nSteps() const{
    
    return nSteps;
}

void MechModel::InitializePETSc(vector<BaseElemMech*> elements){

    // TODO: For debug!
    // for(auto* elem : elements)
    //     for (int dispDof : elem->get_elemDispDof(0))
    //         cout << dispDof << "\n";

    // Initialize the vectors
    VecCreate(PETSC_COMM_WORLD, &vecFext);
    VecSetSizes(vecFext, PETSC_DECIDE, nTotDofs);

    // Since we are interested in sequential implementation for now.
    VecSetType(vecFext, VECSEQ);
    VecSet(vecFext, 0.0); // Set all values to zero.

    // Initialize the displacement increment (solution) vector
    VecDuplicate(vecFext, &vecDeltaDisp);      
    VecSet(vecDeltaDisp, 0.0); 
    
    // Initialize the displacement vector
    VecDuplicate(vecFext, &vecDisp);      
    VecSet(vecDisp, 0.0); 
    
    // Initialize the Fint vector
    VecDuplicate(vecFext, &vecFint);      
    VecSet(vecFint, 0.0); 
    
    // Initialize the residual vector
    VecDuplicate(vecFext, &vecR);         
    VecSet(vecR, 0.0); 

    // Initialize the coefficient matrix.
    MatCreate(PETSC_COMM_WORLD, &matA);
    MatSetSizes(matA, PETSC_DECIDE, PETSC_DECIDE, nTotDofs, nTotDofs);
    
    // Since we are interested in only sequential in this implementation. 
    MatSetType(matA, MATSEQAIJ);
    // MatSetFromOptions(matA);  // for command line options, but we dont do it here.

    // Preallocate the coefficient matrix.
    vector<unordered_set<int>> gDofs(nTotDofs); // vector to store dofs per row.
    int itotv, jtotv; // for global row and colum dof.
    PetscInt *nnz; // Array for the number of zeros per row
    PetscMalloc1(nTotDofs, &nnz); // Allocates the size of nnz
    
    // Find the number of non zeros (columns) per row for preallocation.
    // TODO: Make a check for elements.size()==nElementSets?
    for (auto* elem : elements){  // Loop through element sets

        nElDispDofs = elem->get_nElDispDofs();
        nElements = elem->get_nElements();
        const vector<vector<int>>& elemDispDof_ptr = elem->get_elemDispDof();
    
        for (int iElem=0; iElem<nElements; iElem++){ // Loop through all elements per element set

            for (int idof=0; idof<nElDispDofs; idof++){ // Row

                itotv = elemDispDof_ptr.at(iElem).at(idof); // global row number
                
                for (int jdof=0; jdof<nElDispDofs; jdof++){// Column.

                    jtotv = elemDispDof_ptr.at(iElem).at(jdof); // global column number

                    // If jtotv is not in gDofs.at(itotv)
                    if (find(gDofs.at(itotv).begin(), gDofs.at(itotv).end(), jtotv) == gDofs.at(itotv).end()){

                        gDofs[itotv].insert(jtotv); // Insert column index (automatically handles duplicates)

                    }
                }
            }
        }
    }   

    // Add the number of non zero columns to nnz
    for (int iDof=0; iDof<nTotDofs; iDof++){
        nnz[iDof] = gDofs.at(iDof).size();
    }


    // // To account for Dirichlet BCs in preallocation. Doesn't work very great though. 
    // for (int iPresDof=0; iPresDof<nPresDofs; iPresDof++){
    //     PetscInt row = presDofs[iPresDof];

    //     // Ensure the row has at least one nonzero (for boundary conditions)
    //     nnz[row] = max(nnz[row], static_cast<PetscInt>(1));
    // }

    // Preallocate the stiffness matrix.
    MatSeqAIJSetPreallocation(matA, PETSC_DEFAULT, nnz); 
    PetscFree(nnz);

    // No new memory is allocated
    MatSetOption(matA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    MatSetOption(matA,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

    // // TODO: For debugging. Check the preallocation indirectly by querying the number of nonzeros
    // MatInfo info;
    // MatGetInfo(matA, MAT_LOCAL, &info);
    // if (info.nz_allocated > 0) {
    //     PetscPrintf(PETSC_COMM_WORLD, "Matrix is preallocated.\n");
    // } else {
    //     PetscPrintf(PETSC_COMM_WORLD, "Matrix is not preallocated.\n");
    // }

    // Initialize the `SNES` solver.
    SNESCreate(PETSC_COMM_WORLD, &snes);

    // Parameters for SNES convergence. Details can be found in 
    // https://petsc.org/main/manualpages/SNES/SNESSetTolerances/
    SNESSetTolerances(snes, 1e-50, 1e-8, 1e-8, 100, 1000);

    // Get KSP from SNES
    SNESGetKSP(snes, &ksp);

    // Get PC from KSP
    KSPGetPC(ksp, &pc);

    // PCSetType(pc, PCLU);
    // PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

    KSPSetType(ksp, KSPGMRES);
    KSPGMRESSetRestart(ksp, 50); // Optional: Set GMRES restart value
    KSPSetTolerances(ksp, 1e-12, 1e-12, PETSC_DEFAULT, 1000);
    PCSetType(pc, PCILU);
    PCFactorSetLevels(pc, 2); // ILU with limited fill

    KSPSetFromOptions(ksp);
    PCSetFromOptions(pc);
    SNESSetFromOptions(snes);
}

void MechModel::CalcElemStiffMatx(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats){

    try{

        for (int iSet=0; iSet<nElementSets; iSet++){
            
            if (typeid(*mats[iSet]) == typeid(LinearElastic)){ // Becuase some material models inherit from `Elastic`

                LinearElastic* elasticMat = dynamic_cast<LinearElastic*>(mats[iSet]);
                elements[iSet]->CalcElemStiffMatx(elasticMat->getDMatx());
    
            } else if (typeid(*mats[iSet]) == typeid(IsoHard)){

                IsoHard*plasticMat = dynamic_cast<IsoHard*>(mats[iSet]);
                elements[iSet]->CalcElemStiffMatx(plasticMat->getDMatx());

            } else {

                throw std::runtime_error("Undefined material model < " + std::string(typeid(*mats[iSet]).name()) + " >");
            }

        }

    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}

PetscErrorCode MechModel::AssembleElementMatrix(const auto* elStiffMatx_ptr,
                                                const std::vector<std::vector<int>>& elemDispDof_ptr,
                                                PetscInt nElDispDofs,
                                                PetscInt nElements) {
    const auto& elStiffMatx_ref = *elStiffMatx_ptr; // Dereference the pointer

    // Allocate memory for indices
    PetscInt* i1;
    PetscInt* j1;
    PetscMalloc1(nElDispDofs, &i1);
    PetscMalloc1(nElDispDofs, &j1);

    for (int iElem = 0; iElem < nElements; iElem++) {
        // Get the DOFs associated with the element
        for (int iElDof = 0; iElDof < nElDispDofs; iElDof++) {
            i1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
            j1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
        }

        // Add the local stiffness matrix to the global matrix
        MatSetValues(matA, nElDispDofs, i1, nElDispDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
    }

    // Free memory
    CHKERRQ(PetscFree(i1)); // Check return value of PetscFree
    CHKERRQ(PetscFree(j1)); // Check return value of PetscFree

    return 0;
}

PetscErrorCode MechModel::Assemble(std::vector<BaseElemMech*> elements) {
    
    // Has to be zeroed for iterative solver. 
    MatZeroEntries(matA);

    for (auto* elem : elements) {  // Loop through element sets
        // Element properties
        const auto& elemDispDof_ptr = elem->get_elemDispDof();
        PetscInt nElDispDofs = elem->get_nElDispDofs();
        PetscInt nElements = elem->get_nElements();

        // Get the stiffness matrix variant
        const T_ElStiffMatx& T_elStiffMatx_ref = elem->getElStiffMatx();

        // Use std::visit to handle the variant
        std::visit(
            [&](auto&& elStiffMatx_ptr) {
                if (elStiffMatx_ptr) {
                    AssembleElementMatrix(elStiffMatx_ptr, elemDispDof_ptr, nElDispDofs, nElements);
                } else {
                    throw std::runtime_error("Unsupported element type in T_elStiffMatx_ref");
                }
            },
            T_elStiffMatx_ref);
    }

    // Final assembly
    MatAssemblyBegin(matA, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(matA, MAT_FINAL_ASSEMBLY);

    // Apply Dirichlet boundary conditions (requires MAT_FINAL_ASSEMBLY)
    MatZeroRows(matA, nPresDofs, presDofs, 1.0, NULL, NULL);

    // Save the matrix to MatrixMarket format
    // PetscViewer viewer;
    // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix_output.mtx", &viewer);
    // PetscViewerPushFormat(viewer,  PETSC_VIEWER_ASCII_MATRIXMARKET);
    // PetscViewerPushFormat(viewer,  PETSC_VIEWER_ASCII_MATLAB);

    // MatView(matA, viewer);
    // PetscViewerDestroy(&viewer);

    return 0;
}

void MechModel::SolveSNES(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats, int iStep){

    // Set counter to zero.
    iterCounter = 0; 

    // Set to zero.
    VecSet(vecDeltaDisp, 0.0); 

    KSPSetOperators(ksp, matA, matA);
    KSPSolve(ksp, vecFext, vecDeltaDisp);

    // Create a context for PETSc
    AppCtx *user = new AppCtx{elements, mats, iStep, this};

    // Set the residual function and context
    SNESSetFunction(snes, vecR, ResidualCallback, user); 

    // Set the Jacobian
    SNESSetJacobian(snes, matA, matA, JacobianCallback, user); 

    // Update Jacobian and preconditioner every `NR_freq` iterations
    SNESSetLagJacobian(snes, NR_freq); 
    SNESSetLagPreconditioner(snes, NR_freq); 
    
    // Solve
    SNESSolve(snes, NULL, vecDeltaDisp);

    // Check convergence
    ierr = SNESGetConvergedReason(snes, &reason);
    if (reason < 0) { // Convergence failure
        PetscPrintf(PETSC_COMM_WORLD, "SNES failed to converge: reason %d\n", reason);
        throw std::runtime_error("SNES failed to converge");
    }

    // Update total displacement
    VecAXPY(vecDisp, 1.0, vecDeltaDisp);

    // Calculate total strain
    VecGetArrayRead(vecDisp, &globalBuffer);
    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->CalcElStran(globalBuffer);
    }
    VecRestoreArrayRead(vecDisp, &globalBuffer);

    // Assign values to old
    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->getNew();
    }
}

PetscErrorCode MechModel::ResidualCallback(SNES snes, Vec deltaU, Vec R, void *ctx){
    
    // Cast the context to AppCtx
    AppCtx *user = static_cast<AppCtx*>(ctx);

    // Compute the internal forces using your existing CalcResidual method
    PetscErrorCode ierr = user->mechModel->CalcResidual(deltaU, R, user->elements, user->mats, user->iStep);
    CHKERRQ(ierr);

    // // TODO: For debugging !
    // VecView(R, PETSC_VIEWER_STDOUT_WORLD);

    return 0;
}

PetscErrorCode MechModel::JacobianCallback(SNES snes, Vec deltaU, Mat J, Mat P, void *ctx){
    
    // Cast the context to AppCtx
    AppCtx *user = static_cast<AppCtx*>(ctx);

    // Assemble the Jacobian (stiffness) matrix.
    PetscErrorCode ierr = user->mechModel->Assemble(user->elements);
    CHKERRQ(ierr);

    // // TODO: For debugging. 
    // MatView(J, PETSC_VIEWER_STDOUT_WORLD);

    return 0;
}

PetscErrorCode MechModel::CalcResidual(Vec deltaU, Vec R, vector<BaseElemMech*> elements, vector<BaseMechanics*> mats, int iStep){

        try{
        for (int iSet=0; iSet<nElementSets; iSet++){
            
            if (typeid(*mats[iSet]) == typeid(LinearElastic)){ // Becuase some material models inherit from `Elastic`

    
            } else if (typeid(*mats[iSet]) == typeid(IsoHard)){

                // // TODO: For debug!
                // cout << std::string(typeid(*mats[iSet]).name()) << "\n"; 

                updateStiffMat = iterCounter % NR_freq == 0;

                // Calculate strain increment
                VecGetArrayRead(deltaU, &globalBuffer);
                elements[iSet]->CalcElDStran(globalBuffer);
                VecRestoreArrayRead(deltaU, &globalBuffer);

                // Retrun mapping
                elements[iSet]->CalcRetrunMapping(mats[iSet], updateStiffMat, iStep);

                // Has to be set to zero before any calculation, otherwise it accumulates.
                for (int iDof=0; iDof<nTotDofs; iDof++){
                    Fint[iDof] = 0;
                }
                // Calculate internal force
                elements[iSet]->CalcFint(Fint);
                VecSetValues(vecFint, nTotDofs, indices, Fint, INSERT_VALUES); 
                VecAssemblyBegin(vecFint); VecAssemblyEnd(vecFint);

                // VecSetValues(vecFext, nPresDofs, indices, presZeros, INSERT_VALUES); 
                // VecAssemblyBegin(vecFext); VecAssemblyEnd(vecFext);
                /* 
                << NOTE >> We are solving the system Fext-Fint = R, while PETSc solves Ju = -R. 
                It will multiply R with -1, so we keep this in mind while providing R. 
                */
                // VecWAXPY(R, -1.0, vecFext, vecFint);
                VecAXPBY(R, 1.0, 0.0, vecFint);
                VecSetValues(R, nPresDofs, presDofs, presZeros, INSERT_VALUES); 
                VecAssemblyBegin(R); VecAssemblyEnd(R);

                // Update iteration counter.
                iterCounter++;

                // // TODO: For debugging 
                // VecView(deltaU, PETSC_VIEWER_STDOUT_WORLD);

                // PetscReal l2norm;
                // VecNorm(R, NORM_2, &l2norm);
                // cout << "-- L2Norm: " << l2norm << "\n" ;


            } else {

                throw std::runtime_error("Undefined material model < " + std::string(typeid(*mats[iSet]).name()) + " >");
            }
        }

    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

Vec& MechModel::getB(){

    return vecFext;
}

Vec& MechModel::getX(){

    return vecDisp;
}

Mat& MechModel::getA(){

    return matA;
}

void MechModel::CalcStres(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats){

    VecGetArrayRead(vecDisp, &globalBuffer);

    for (int iSet=0; iSet<nElementSets; iSet++){
        
        elements[iSet]->CalcStres(mats[iSet]->getDMatx(), globalBuffer, Fint, nodStres, nodStran, nodCount);
    }

    VecRestoreArrayRead(vecDisp, &globalBuffer);

    // Number averaging the nodal values
    if (nDim==2){

        for(int iNod=0; iNod<nTotNodes; iNod++){
            
            std::get<std::vector<ColVecd3>>(nodStran).at(iNod) = std::get<std::vector<ColVecd3>>(nodStran).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd3>>(nodStres).at(iNod) = std::get<std::vector<ColVecd3>>(nodStres).at(iNod)/nodCount.at(iNod);
        }

    } else if (nDim==3){

        for(int iNod=0; iNod<nTotNodes; iNod++){
            
            std::get<std::vector<ColVecd6>>(nodStran).at(iNod) = std::get<std::vector<ColVecd6>>(nodStran).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd6>>(nodStres).at(iNod) = std::get<std::vector<ColVecd6>>(nodStres).at(iNod)/nodCount.at(iNod);
        }
    }

    // TODO: For debug!
    // cout << std::get<std::vector<ColVecd3>>(nodStran).at(0) << "\n";

}

void MechModel::CalcNodVals(vector<BaseElemMech*> elements){

    for (int iSet=0; iSet<nElementSets; iSet++){
        
        elements[iSet]->CalcNodVals(nodStres, nodStran, nodStran_e, nodStran_p, nodStran_eq, nodStres_eq, nodStres_h, nodRho, nodCount);
    }

    // Number averaging the nodal values
    if (nDim==2){

        for(int iNod=0; iNod<nTotNodes; iNod++){
            
            std::get<std::vector<ColVecd3>>(nodStran).at(iNod) = std::get<std::vector<ColVecd3>>(nodStran).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd3>>(nodStran_e).at(iNod) = std::get<std::vector<ColVecd3>>(nodStran_e).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd3>>(nodStran_p).at(iNod) = std::get<std::vector<ColVecd3>>(nodStran_p).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd3>>(nodStres).at(iNod) = std::get<std::vector<ColVecd3>>(nodStres).at(iNod)/nodCount.at(iNod);
            nodStran_eq.at(iNod) = nodStran_eq.at(iNod)/nodCount.at(iNod);
            nodStres_eq.at(iNod) = nodStres_eq.at(iNod)/nodCount.at(iNod);
            nodStres_h.at(iNod) = nodStres_h.at(iNod)/nodCount.at(iNod);
            nodRho.at(iNod) = nodRho.at(iNod)/nodCount.at(iNod);
        }

    } else if (nDim==3){

        for(int iNod=0; iNod<nTotNodes; iNod++){
            
            std::get<std::vector<ColVecd6>>(nodStran).at(iNod) = std::get<std::vector<ColVecd6>>(nodStran).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd6>>(nodStran_e).at(iNod) = std::get<std::vector<ColVecd6>>(nodStran_e).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd6>>(nodStran_p).at(iNod) = std::get<std::vector<ColVecd6>>(nodStran_p).at(iNod)/nodCount.at(iNod);
            std::get<std::vector<ColVecd6>>(nodStres).at(iNod) = std::get<std::vector<ColVecd6>>(nodStres).at(iNod)/nodCount.at(iNod);
            nodStran_eq.at(iNod) = nodStran_eq.at(iNod)/nodCount.at(iNod);
            nodStres_eq.at(iNod) = nodStres_eq.at(iNod)/nodCount.at(iNod);
            nodStres_h.at(iNod) = nodStres_h.at(iNod)/nodCount.at(iNod);
            nodRho.at(iNod) = nodRho.at(iNod)/nodCount.at(iNod);
        }
    }

    // TODO: For debug!
    // cout << std::get<std::vector<ColVecd3>>(nodStran).at(0) << "\n";

}

void MechModel::WriteOut(vector<BaseElemMech*> elements, H5IO &H5File_out, const string iStep){

    // Displacement
    VecGetArrayRead(vecDisp, &globalBuffer);
    H5File_out.WriteArray1D("Disp/Step_"+iStep, nTotDofs, globalBuffer);
    VecRestoreArrayRead(vecDisp, &globalBuffer);
    
    // Force
    H5File_out.WriteArray1D("Force/Step_"+iStep, nTotDofs, Fint);

    // Equivalent stress and strain
    H5File_out.WriteArray1D("Stress_eq/Step_"+iStep, nTotNodes, nodStres_eq.data());
    H5File_out.WriteArray1D("Strain_eq/Step_"+iStep, nTotNodes, nodStran_eq.data());
    H5File_out.WriteArray1D("Stress_h/Step_"+iStep, nTotNodes, nodStres_h.data());
    
    H5File_out.WriteArray1D("Rho/Step_"+iStep, nTotNodes, nodRho.data()); // Will write it as output anyways for now


    // Stresses and strains
    if (nDim==2){
        H5File_out.WriteTensor("Strain/Step_"+iStep, nTotNodes, 3, nodStran);
        H5File_out.WriteTensor("Strain_e/Step_"+iStep, nTotNodes, 3, nodStran_e);
        H5File_out.WriteTensor("Strain_p/Step_"+iStep, nTotNodes, 3, nodStran_p);
        H5File_out.WriteTensor("Stress/Step_"+iStep, nTotNodes, 3, nodStres);
    } else if (nDim==3) {
        H5File_out.WriteTensor("Strain/Step_"+iStep, nTotNodes, 6, nodStran);
        H5File_out.WriteTensor("Strain_e/Step_"+iStep, nTotNodes, 6, nodStran_e);
        H5File_out.WriteTensor("Strain_p/Step_"+iStep, nTotNodes, 6, nodStran_p);
        H5File_out.WriteTensor("Stress/Step_"+iStep, nTotNodes, 6, nodStres);
    }

    // Set nodal vlaues to zero
    setZeroNodVals();
}