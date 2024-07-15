#include<iostream>
#include<variant>

#include "Models/TrappingModel.h"

using namespace std;

TrappingModel::TrappingModel(vector<BaseElemTrap*> elements, H5IO& H5File_in){

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
    dsetName = "SimulationParameters/dt";
    dt = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/T";
    T = H5File_in.ReadScalar(dsetName);

    // Set the type and size
    nodCount.resize(nTotNodes);
    if (nDim == 2){ // Case 2D model

        nodFlux = vector<ColVecd2>(nTotNodes);

    } else if (nDim == 3){ // Case 3D

        nodFlux = vector<ColVecd3>(nTotNodes);

    }

    // Initialize to zeros
    setZero_nodFlux();

    InitializePETSc(elements);
}

TrappingModel::~TrappingModel(){

    // Deallocate memory.
    PetscFree(presDofs); PetscFree(presVals);
    VecDestroy(&F); VecDestroy(&x); MatDestroy(&K); MatDestroy(&M);
    // Finalize PETSc
    PetscFinalize();
    // Exit message
    cout << "TrappingModel model exited correctly" << "\n";
}

void TrappingModel::setZero_nodFlux(){

    if (nDim == 2){ // Case 2D model

        for(int iNod=0; iNod<nTotNodes; iNod++){

            std::get<std::vector<ColVecd2>>(nodFlux).at(iNod).setZero();
            nodCount.at(iNod) = 0;

        }

    } else if (nDim == 3){ // Case 3D

        for(int iNod=0; iNod<nTotNodes; iNod++){

            std::get<std::vector<ColVecd3>>(nodFlux).at(iNod).setZero();
            nodCount.at(iNod) = 0;
            
        }
    }
}

void TrappingModel::InitializePETSc(vector<BaseElemTrap*> elements){

    // TODO: For debug!
    // for(auto* elem : elements)
    //     for (int dispDof : elem->get_elemDispDof(0))
    //         cout << dispDof << "\n";

    // Initialize the vectors
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, PETSC_DECIDE, nTotDofs);

    VecSetType(F, VECSEQ);
    VecDuplicate(F, &x);      

    VecSet(F, 0.0); // Set all values to zero.
    VecAssemblyBegin(F); VecAssemblyEnd(F);

    VecSet(x, 0.0); // Set all values to zero.
    VecAssemblyBegin(x); VecAssemblyEnd(x);

    // Initialize the coefficient matrices.
    MatCreate(PETSC_COMM_WORLD, &K);
    MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, nTotDofs, nTotDofs);
    
    MatSetType(K, MATSEQAIJ); 

    MatDuplicate(K, MAT_DO_NOT_COPY_VALUES, &M);

    // Preallocate the coefficient matrix.
    vector<vector<int>> gDofs(nTotDofs); // vector to store dofs per row.
    int itotv, jtotv; // for global row and colum dof.
    PetscInt *nnz; // Array for the number of zeros per row
    PetscMalloc1(nTotDofs, &nnz); // Allocates the size of nnz
    
    // Find the number of non zeros (columns) per row for preallocation.
    for (auto* elem : elements){  // Loop through element sets

        nElConDofs = elem->get_nElConDofs();
        nElements = elem->get_nElements();
        const vector<vector<int>>& elemConDof_ptr = elem->get_elemConDof();
    
        for (int iElem=0; iElem<nElements; iElem++){ // Loop through all elements per element set

            for (int idof=0; idof<nElConDofs; idof++){ // Row

                itotv = elemConDof_ptr.at(iElem).at(idof); // global row number
                
                for (int jdof=0; jdof<nElConDofs; jdof++){// Column.

                    jtotv = elemConDof_ptr.at(iElem).at(jdof); // global column number

                    // If jtotv is not in gDofs.at(itotv)
                    if (find(gDofs.at(itotv).begin(), gDofs.at(itotv).end(), jtotv) == gDofs.at(itotv).end()){

                        gDofs.at(itotv).push_back(jtotv);
                    }
                }
            }
        }
    }   

    // Add the number of non zero columns to nnz
    for (int iDof=0; iDof<nTotDofs; iDof++){
        nnz[iDof] = gDofs.at(iDof).size();
    }

    // Preallocate the stiffness matrix.
    MatSeqAIJSetPreallocation(K, PETSC_DEFAULT, nnz); 
    MatSeqAIJSetPreallocation(M, PETSC_DEFAULT, nnz); 
    PetscFree(nnz);

}

void TrappingModel::WriteGradPhi(vector<BaseElemTrap*> elements, H5IO& H5File_out){

    // Nodal laplacian of phi
    double* nodLapPhi = new double[nTotNodes];
    for (int iNod=0; iNod<nTotDofs; iNod++){
        nodLapPhi[iNod] = 0;}
    // Nodal gradient of phi.
    T_nodStres nodGradPhi;     
    // Counter for integration points surrounding nodes.
    vector<double> nodPhiCount(nTotNodes);     
     
    if (nDim == 2){ // Case 2D model
        nodGradPhi = vector<ColVecd2>(nTotNodes);
    } else if (nDim == 3){ // Case 3D
        nodGradPhi = vector<ColVecd3>(nTotNodes);
    }

    // Set zeros
    if (nDim == 2){ // Case 2D model
        for(int iNod=0; iNod<nTotNodes; iNod++){
            std::get<std::vector<ColVecd2>>(nodGradPhi).at(iNod).setZero();
            nodPhiCount.at(iNod) = 0;
        }
    } else if (nDim == 3){ // Case 3D
        for(int iNod=0; iNod<nTotNodes; iNod++){
            std::get<std::vector<ColVecd3>>(nodGradPhi).at(iNod).setZero();
            nodPhiCount.at(iNod) = 0;
        }
    }

    elements[0]->CalcGrad(nodGradPhi, nodPhiCount, nodLapPhi);

    // Number averaging the nodal values
    if (nDim==2){
        for(int iNod=0; iNod<nTotNodes; iNod++){
            std::get<std::vector<ColVecd2>>(nodGradPhi).at(iNod) = std::get<std::vector<ColVecd2>>(nodGradPhi).at(iNod)/nodPhiCount.at(iNod);
        }
    } else if (nDim==3){
        for(int iNod=0; iNod<nTotNodes; iNod++){
            std::get<std::vector<ColVecd3>>(nodGradPhi).at(iNod) = std::get<std::vector<ColVecd3>>(nodGradPhi).at(iNod)/nodPhiCount.at(iNod);
        }
    }

    // Flux
    if (nDim==2){
        H5File_out.WriteStres("GradPhi", nTotNodes, 2, nodGradPhi);
    } else if (nDim==3) {
        H5File_out.WriteStres("GradPhi", nTotNodes, 3, nodGradPhi);
    }

    // laplacian
    H5File_out.WriteArray_1D("LapPhi", nTotDofs, nodLapPhi);
    delete [] nodLapPhi;
}

void TrappingModel::CalcElemStiffMatx(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats){

    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->CalcElemStiffMatx(mats[iSet], T, globalBuffer);
    }
}

void TrappingModel::Assemble(vector<BaseElemTrap*> elements){

    for (auto* elem : elements){  // Loop through element sets

        const vector<vector<int>>& elemConDof_ptr = elem->get_elemConDof();
        nElConDofs = elem->get_nElConDofs();
        nElements = elem->get_nElements();

        //  Assemble the coefficient matrix.        
        PetscInt* i1; PetscInt* j1; // Indices for row and columns to insert.
        PetscMalloc1(nElConDofs, &i1);
        PetscMalloc1(nElConDofs, &j1);

        const T_ElStiffMatx& T_elStiffMatx_ref = elem->getElStiffMatx();
        const T_ElStiffMatx& T_elCapMatx_ref = elem->getElCapMatx();

        if (std::holds_alternative<vector<Matd4x4>*>(T_elStiffMatx_ref)){  // Quad4 elements.
 
            const vector<Matd4x4>& elStiffMatx_ref = *std::get<vector<Matd4x4>*>(T_elStiffMatx_ref);
            const vector<Matd4x4>& elCapMatx_ref = *std::get<vector<Matd4x4>*>(T_elCapMatx_ref);

            for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                // Get the con dofs associated with the element
                for(int iElDof=0; iElDof<nElConDofs; iElDof++){

                    i1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
                    j1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);

                }

                MatSetValues(K, nElConDofs, i1, nElConDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
                MatSetValues(M, nElConDofs, i1, nElConDofs, j1, elCapMatx_ref.at(iElem).data(), ADD_VALUES);
            }

        } else if (std::holds_alternative<vector<Matd3x3>*>(T_elStiffMatx_ref)){  // Tri3 elements.
 
            const vector<Matd3x3>& elStiffMatx_ref = *std::get<vector<Matd3x3>*>(T_elStiffMatx_ref);
            const vector<Matd3x3>& elCapMatx_ref = *std::get<vector<Matd3x3>*>(T_elCapMatx_ref);

            for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                // Get the con dofs associated with the element
                for(int iElDof=0; iElDof<nElConDofs; iElDof++){

                    i1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
                    j1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
                }
                
                MatSetValues(K, nElConDofs, i1, nElConDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
                MatSetValues(M, nElConDofs, i1, nElConDofs, j1, elCapMatx_ref.at(iElem).data(), ADD_VALUES);
            }
        }

        // Free memory
        PetscFree(i1);
        PetscFree(j1);
    }

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    MatSetOption(K, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);

    // TODO: For debug!
    // // Extract a specific row
    // PetscInt ROW = 4676; // for example, extracting row 2
    // const PetscInt *COL;
    // const PetscScalar *VAL;
    // PetscInt NCOL;

    // MatGetRow(K, ROW, &NCOL, &COL, &VAL);

    // for (int w = 0; w < NCOL; w++) {
    //     cout << "Col: " << COL[w] << ", Val: " << VAL[w] << "\n";
    // }

    // // Restore the row
    // MatRestoreRow(K, ROW, &NCOL, &COL, &VAL); 

    // MatView(K, PETSC_VIEWER_STDOUT_WORLD);


    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    MatSetOption(M, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
}

void TrappingModel::InitializeDirichBC(H5IO& H5File_in){

    // Read Dirichlet BCs
    string dsetName;
    dsetName = "SimulationParameters/nPresDofs";
    nPresDofs = H5File_in.ReadScalar(dsetName);
    PetscMalloc1(nPresDofs, &presDofs);
    PetscMalloc1(nPresDofs, &presVals); 

    vector<double> dummy(3);
    for (int iPresDof=0; iPresDof<nPresDofs; iPresDof++){
        // Read values
        dsetName = "PrescribedDOFs/Prescribed_"+to_string(iPresDof);
        H5File_in.ReadFieldDoub1D(dsetName, dummy);
        // Assign values
        presDofs[iPresDof] = dummy.at(0); // nDim*iNode+dof
        presVals[iPresDof] = dummy.at(1);
        // // TODO: For debug!
        // cout << presDofs[iPresDof] << " --> " << presVals[iPresDof] << "\n";
    }

    // For the coefficient matrix
    MatZeroRows(K, nPresDofs, presDofs, 1.0, NULL, NULL);
    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
}

void TrappingModel::setDirichBC(){
    
    MatMult(M, x, F);

    VecSetValues(F, nPresDofs, presDofs, presVals, INSERT_VALUES); 
    VecAssemblyBegin(F); VecAssemblyEnd(F);

    // // TODO: For debug!
    // VecView(F, PETSC_VIEWER_STDOUT_WORLD);
}

int TrappingModel::get_nSteps() const{
    
    return nSteps;
}

Vec& TrappingModel::getF(){

    return F;
}

Vec& TrappingModel::getX(){

    return x;
}

Mat& TrappingModel::getK(){

    return K;
}

// void TrappingModel::CalcFlux(vector<BaseElemTrap*> elements, vector<void CalcFlux(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats);
// *> mats){

//     VecGetArrayRead(x, &globalBuffer);

//     for (int iSet=0; iSet<nElementSets; iSet++){
        
//         elements[iSet]->CalcFlux(mats[iSet]->getKMatx(), globalBuffer, nodFlux, nodCount);
//     }

//     // Number averaging the nodal values
//     if (nDim==2){

//         for(int iNod=0; iNod<nTotNodes; iNod++){
            
//             std::get<std::vector<ColVecd2>>(nodFlux).at(iNod) = std::get<std::vector<ColVecd2>>(nodFlux).at(iNod)/nodCount.at(iNod);
//         }

//     } else if (nDim==3){

//         for(int iNod=0; iNod<nTotNodes; iNod++){
            
//             std::get<std::vector<ColVecd3>>(nodFlux).at(iNod) = std::get<std::vector<ColVecd3>>(nodFlux).at(iNod)/nodCount.at(iNod);
//         }

//     }

//     // TODO: For debug!
//     // cout << std::get<std::vector<ColVecd3>>(nodStran).at(0) << "\n";

//     VecRestoreArrayRead(x, &globalBuffer);
// }

void TrappingModel::WriteOut(vector<BaseElemTrap*> elements, H5IO &H5File_out, const string iStep){

    // Displacements
    VecGetArrayRead(x, &globalBuffer);
    H5File_out.WriteArray_1D("Con/Step_"+iStep, nTotDofs, globalBuffer);
    VecRestoreArrayRead(x, &globalBuffer);

    // // Flux
    // if (nDim==2){
    //     H5File_out.WriteStres("Flux/Step_"+iStep, nTotNodes, 2, nodFlux);
    // } else if (nDim==3) {
    //     H5File_out.WriteStres("Flux/Step_"+iStep, nTotNodes, 3, nodFlux);
    // }

    // // set zeros
    // setZero_nodFlux();
}