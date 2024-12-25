#include<iostream>
#include<variant>
#include <unordered_set>

#include "Models/MechModel.h"
#include "Materials/Mechanics/IsoHard.h"

using namespace std;

MechModel::MechModel(H5IO& H5File_in){

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
    if (nDim == 2){ // Case 2D model

        nodStres = vector<ColVecd3>(nTotNodes);
        nodStran = vector<ColVecd3>(nTotNodes);

    } else if (nDim == 3){ // Case 3D

        nodStres = vector<ColVecd6>(nTotNodes);
        nodStran = vector<ColVecd6>(nTotNodes);

    }

    // Allocate memory for `Fint` and `indices`.
    PetscMalloc1(nTotDofs, &Fint);
    PetscMalloc1(nTotDofs, &indices);

    // Initialize to zeros, otherwise will get garbage memory values.
    setZeroNodVals();

    // Set indices
    for (int iDof=0; iDof<nTotDofs; iDof++){
        indices[iDof] = iDof;
    }

}

MechModel::~MechModel(){

    // Deallocate memory.
    PetscFree(presDofs); PetscFree(presVals); PetscFree(Fint); PetscFree(indices);
    VecDestroy(&vecFext); VecDestroy(&vecDisp); VecDestroy(&vecDeltaDisp); VecDestroy(&vecFint); 
    VecDestroy(&vecR); MatDestroy(&matA);
    SNESDestroy(&snes);
    // Finalize PETSc
    PetscFinalize();
    // Exit message
    cout << "MechModel elements exited correctly" << "\n";
}

void MechModel::setZeroNodVals(){

    for (int iDof=0; iDof<nTotDofs; iDof++){
        Fint[iDof] = 0;
    }

    if (nDim == 2){ // Case 2D model

        for(int iNod=0; iNod<nTotNodes; iNod++){

            std::get<std::vector<ColVecd3>>(nodStran).at(iNod).setZero();
            std::get<std::vector<ColVecd3>>(nodStres).at(iNod).setZero();
            nodCount.at(iNod) = 0;
            
        }

    } else if (nDim == 3){ // Case 3D

        for(int iNod=0; iNod<nTotNodes; iNod++){

            std::get<std::vector<ColVecd6>>(nodStran).at(iNod).setZero();
            std::get<std::vector<ColVecd6>>(nodStres).at(iNod).setZero();
            nodCount.at(iNod) = 0;

        }
    }
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
    // VecSetFromOptions(vecFext); =// for the general case.
    VecSet(vecFext, 0.0); // Set all values to zero.
    VecAssemblyBegin(vecFext); VecAssemblyEnd(vecFext);
    
    // Initialize the solution vector
    VecDuplicate(vecFext, &vecDisp);      
    VecSet(vecDisp, 0.0); 
    VecAssemblyBegin(vecDisp); VecAssemblyEnd(vecDisp);

    // Initialize the displacement increment vector
    VecDuplicate(vecFext, &vecDeltaDisp);      
    VecSet(vecDisp, 0.0); 
    VecAssemblyBegin(vecDisp); VecAssemblyEnd(vecDisp);
    
    // Initialize the Fint vector
    VecDuplicate(vecFext, &vecFint);      
    VecSet(vecFint, 0.0); 
    VecAssemblyBegin(vecFint); VecAssemblyEnd(vecFint);
    
    // Initialize the residual vector
    VecDuplicate(vecFext, &vecR);         
    VecSet(vecR, 0.0); 
    VecAssemblyBegin(vecR); VecAssemblyEnd(vecR);

    // Initialize the coefficient matrix.
    MatCreate(PETSC_COMM_WORLD, &matA);
    MatSetSizes(matA, PETSC_DECIDE, PETSC_DECIDE, nTotDofs, nTotDofs);
    
    // Since we are interested in only sequential in this implementation. 
    MatSetType(matA, MATSEQAIJ);
    // MatSetFromOptions(matA);  // for command line options, but we dont do it here.

    // // MAT_SYMMETRIC: symmetric in terms of both structure and value
    // MatSetOption(matA, MAT_SYMMETRIC, PETSC_TRUE);

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

    // To account for Dirichlet BCs in preallocation. Doesn't work very great though. 
    for (int iPresDof=0; iPresDof<nPresDofs; iPresDof++){
        PetscInt row = presDofs[iPresDof];

        // Ensure the row has at least one nonzero (for boundary conditions)
        nnz[row] = max(nnz[row], static_cast<PetscInt>(1));
    }

    // Preallocate the stiffness matrix.
    MatSeqAIJSetPreallocation(matA, PETSC_DEFAULT, nnz); 
    PetscFree(nnz);

    // // TODO: For debugging. Check the preallocation indirectly by querying the number of nonzeros
    // MatInfo info;
    // MatGetInfo(matA, MAT_LOCAL, &info);
    // if (info.nz_allocated > 0) {
    //     PetscPrintf(PETSC_COMM_WORLD, "Matrix is preallocated.\n");
    // } else {
    //     PetscPrintf(PETSC_COMM_WORLD, "Matrix is not preallocated.\n");
    // }

    // Throw error if unallocated entry is accessed if "PETSC_TRUE".
    // Should be added after the matrix is assembled https://lists.mcs.anl.gov/pipermail/petsc-users/2019-October/039608.html
    // MatSetOption(matA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    MatSetOption(matA, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

    // Initialize the `SNES` solver.
    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFromOptions(snes);
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
        H5File_in.ReadFieldDoub1D(dsetName, dummy);
        // Assign values
        presDofs[iPresDof] = nDim*dummy.at(0)+dummy.at(1); // nDim*iNode+dof
        presVals[iPresDof] = dummy.at(2)/nSteps;
        // // TODO: For debug!
        // cout << presDofs[iPresDof] << " --> " << presVals[iPresDof] << "\n";
        presZeros[iPresDof] = 0;
    }
}

void MechModel::setDirichBC(){

    VecSetValues(vecFext, nPresDofs, presDofs, presVals, ADD_VALUES); 
    VecAssemblyBegin(vecFext); VecAssemblyEnd(vecFext);
}

void MechModel::UpdateDisp(){

    VecSetValues(vecDisp, nPresDofs, presDofs, presVals, ADD_VALUES); 
    VecAssemblyBegin(vecDisp); VecAssemblyEnd(vecDisp);
}

void MechModel::CalcElemStiffMatx(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats){

    try{

        for (int iSet=0; iSet<nElementSets; iSet++){
            
            if (typeid(*mats[iSet]) == typeid(LinearElastic)){ // Becuase some material models inherit from `Elastic`

                LinearElastic* elasticMat = dynamic_cast<LinearElastic*>(mats[iSet]);
                elements[iSet]->CalcElemStiffMatx(elasticMat->getDMatx());
    
            } else if (typeid(*mats[iSet]) == typeid(IsoHard)){

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

PetscErrorCode MechModel::Assemble(vector<BaseElemMech*> elements){

    // Has to be zeroed for iterative solver. 
    MatZeroEntries(matA);

    for (auto* elem : elements){  // Loop through element sets

        // Pointer to the element dip DOFs. Will vanish out the "for" element set loop.
        const vector<vector<int>>& elemDispDof_ptr = elem->get_elemDispDof();
        // Get the number of element dips DOFs. 
        nElDispDofs = elem->get_nElDispDofs();
        // Get the number of elements in the set. 
        nElements = elem->get_nElements();

        //  Assemble the coefficient matrix.        
        PetscInt* i1; PetscInt* j1; // Indices for row and columns to insert.
        PetscMalloc1(nElDispDofs, &i1);
        PetscMalloc1(nElDispDofs, &j1);

        /* matA constatn reference to the std::variant object returned by `elem->getElStiffMatx()`, 
        in this case `T_ElStiffMatx` variant that holds a pointer to a vector.
        You cannot modify the object through T_elStiffMatx_ref (because it's a constant reference).
        */
        const T_ElStiffMatx& T_elStiffMatx_ref = elem->getElStiffMatx();

        // // TODO: For debugging.
        // std::cout << "Active index: " << T_elStiffMatx_ref.index() << std::endl;

        if (std::holds_alternative<vector<Matd8x8>*>(T_elStiffMatx_ref)){  // Quad4 elements.

            /* The "*" operator dereferences the pointer stored in the std::variant, giving
            you a reference to the actual vector. 
            This way avoids copying the vector and can safely read from it without modifying it.
            */ 
            const vector<Matd8x8>& elStiffMatx_ref = *std::get<vector<Matd8x8>*>(T_elStiffMatx_ref);

            for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                // Get the disp dofs associated with the element
                for(int iElDof=0; iElDof<nElDispDofs; iElDof++){

                    i1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                    j1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);

                }

                MatSetValues(matA, nElDispDofs, i1, nElDispDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
            }

        } else if (std::holds_alternative<vector<Matd6x6>*>(T_elStiffMatx_ref)){  // Tri3 elements.
 
            const vector<Matd6x6>& elStiffMatx_ref = *std::get<vector<Matd6x6>*>(T_elStiffMatx_ref);

            for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                // Get the disp dofs associated with the element
                for(int iElDof=0; iElDof<nElDispDofs; iElDof++){

                    i1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                    j1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                }

                MatSetValues(matA, nElDispDofs, i1, nElDispDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
            }

        } else if (std::holds_alternative<vector<Matd24x24>*>(T_elStiffMatx_ref)){  // Hex8 elements.
 
            const vector<Matd24x24>& elStiffMatx_ref = *std::get<vector<Matd24x24>*>(T_elStiffMatx_ref);

            for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                // Get the disp dofs associated with the element
                for(int iElDof=0; iElDof<nElDispDofs; iElDof++){

                    i1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                    j1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                }

                MatSetValues(matA, nElDispDofs, i1, nElDispDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
            }
        }

        // Free memory
        PetscFree(i1);
        PetscFree(j1);
    }

    // Final assembly
    MatAssemblyBegin(matA, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(matA, MAT_FINAL_ASSEMBLY);

    // For Dirichlet boundary conditions (Requires `MAT_FINAL_ASSEMBLY`)
    MatZeroRows(matA, nPresDofs, presDofs, 1.0, NULL, NULL);

    // Sets the final sparsity structure  
    MatSetOption(matA, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);

    return 0;
}

void MechModel::SolveSNES(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats, int iStep){

    // Set counter to zero.
    iterCounter = 0;

    // Create a context for PETSc
    AppCtx *user = new AppCtx{elements, mats, iStep, this};

    // Set the residual function and context
    SNESSetFunction(snes, vecR, ResidualCallback, user); 

    // Set the Jacobian
    SNESSetJacobian(snes, matA, matA, JacobianCallback, user); 
    
    // Solve
    SNESSolve(snes, NULL, vecDeltaDisp);

}

PetscErrorCode MechModel::ResidualCallback(SNES snes, Vec deltaU, Vec R, void *ctx){
    
    // Cast the context to AppCtx
    AppCtx *user = static_cast<AppCtx*>(ctx);

    // Compute the internal forces using your existing CalcResidual method
    PetscErrorCode ierr = user->mechModel->CalcResidual(deltaU, user->elements, user->mats, user->iStep);
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

PetscErrorCode MechModel::CalcResidual(Vec deltaU, vector<BaseElemMech*> elements, vector<BaseMechanics*> mats, int iStep){

        try{
        for (int iSet=0; iSet<nElementSets; iSet++){
            
            if (typeid(*mats[iSet]) == typeid(LinearElastic)){ // Becuase some material models inherit from `Elastic`

    
            } else if (typeid(*mats[iSet]) == typeid(IsoHard)){

                // // TODO: For debug!
                // cout << std::string(typeid(*mats[iSet]).name()) << "\n"; 

                // Update iteration counter.
                iterCounter += 1;

                updateStiffMat = iterCounter % NR_freq == 0;

                // Update total displacement vector.
                VecAXPY(vecDisp, +1.0, deltaU);
                VecAssemblyBegin(vecDisp); VecAssemblyEnd(vecDisp);

                // Calculate total strain
                VecGetArrayRead(vecDisp, &globalBuffer);
                elements[iSet]->CalcElStran(globalBuffer);
                VecRestoreArrayRead(vecDisp, &globalBuffer);

                // Retrun mapping
                elements[iSet]->CalcRetrunMapping(mats[iSet], updateStiffMat, iStep);

                // Calculate internal force
                elements[iSet]->CalcFint(Fint);
                VecSetValues(vecFint, nTotDofs, indices, Fint, INSERT_VALUES); 
                VecAssemblyBegin(vecFint); VecAssemblyEnd(vecFint);

                /* 
                << NOTE >> We are solving the system Fext-Fint = R, while PETSc solves Ju = -R. 
                It will multiply R with -1, so we keep this in mind while providing R. 
                */
                VecWAXPY(vecR, -1.0, vecFext, vecFint);
                VecSetValues(vecR, nPresDofs, presDofs, presZeros, INSERT_VALUES); 
                VecAssemblyBegin(vecR); VecAssemblyEnd(vecR);

                // // TODO: For debugging !
                // VecView(vecFint, PETSC_VIEWER_STDOUT_WORLD);

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

void MechModel::WriteOut(vector<BaseElemMech*> elements, H5IO &H5File_out, const string iStep){

    // Displacements
    VecGetArrayRead(vecDisp, &globalBuffer);
    H5File_out.WriteArray_1D("Disp/Step_"+iStep, nTotDofs, globalBuffer);
    VecRestoreArrayRead(vecDisp, &globalBuffer);
    H5File_out.WriteArray_1D("Force/Step_"+iStep, nTotDofs, Fint);

    // Stresses and strains
    if (nDim==2){
        H5File_out.WriteStres("Strain/Step_"+iStep, nTotNodes, 3, nodStran);
        H5File_out.WriteStres("Stress/Step_"+iStep, nTotNodes, 3, nodStres);
    } else if (nDim==3) {
        H5File_out.WriteStres("Strain/Step_"+iStep, nTotNodes, 6, nodStran);
        H5File_out.WriteStres("Stress/Step_"+iStep, nTotNodes, 6, nodStres);
    }

    // Set nodal vlaues to zero
    setZeroNodVals();
}