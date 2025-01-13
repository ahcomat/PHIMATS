#include<iostream>
#include<variant>
#include <unordered_set>

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
    T0 = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/nExitNodes";
    nExitNodes = H5File_in.ReadScalar(dsetName);

    for (int iSet=0; iSet<nElementSets; iSet++){
        nTotGuasPts += elements[iSet]->get_nGauss()*elements[iSet]->get_nElements();
    }

    intPtFlux = vector<ColVecd3>(nTotGuasPts);

    // Read exit nodes IDs
    ExitNodeIDs.resize(nExitNodes);
    dsetName = "ExitNodes";
    H5File_in.ReadField1D(dsetName, ExitNodeIDs);

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
    MatCreateSeqAIJ(PETSC_COMM_WORLD, nTotDofs, nTotDofs, PETSC_DEFAULT, NULL, &K); // Works better for HPC

    // MatCreate(PETSC_COMM_WORLD, &K);
    // MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, nTotDofs, nTotDofs);
    // MatSetType(K, MATSEQAIJ); 

    MatDuplicate(K, MAT_DO_NOT_COPY_VALUES, &M);

    // Preallocate the coefficient matrix.
    vector<unordered_set<int>> gDofs(nTotDofs); // vector to store dofs per row.
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

    // Preallocation.
    MatSeqAIJSetPreallocation(K, PETSC_DEFAULT, nnz);
    // No new memory is allocated
    MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    MatSetOption(K,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

    MatSeqAIJSetPreallocation(M, PETSC_DEFAULT, nnz); 
    // No new memory is allocated
    MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    MatSetOption(M,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

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
        H5File_out.WriteTensor("GradPhi", nTotNodes, 2, nodGradPhi);
    } else if (nDim==3) {
        H5File_out.WriteTensor("GradPhi", nTotNodes, 3, nodGradPhi);
    }

    // laplacian
    H5File_out.WriteArray1D("LapPhi", nTotDofs, nodLapPhi);
    delete [] nodLapPhi;
}

void TrappingModel::UpdateTemp(const int iStep, double HR){

    // T = dt*(double)iStep*HR + T0;

    T += dt*HR;
}

void TrappingModel::Update_dt(vector<BaseElemTrap*> elements, double dtNew){

    dt = dtNew;

    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->set_dt(dt);
    }
}

void TrappingModel::WriteTemp(H5IO &H5File_out, const int iStep){

    H5File_out.WriteScalar("Temp/Step_"+to_string(iStep), T);
}

void TrappingModel::CalcElemStiffMatx(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats){

        for (int iSet=0; iSet<nElementSets; iSet++){
            elements[iSet]->CalcElemStiffMatx(mats[iSet], T);
        }
}

void TrappingModel::Assemble(vector<BaseElemTrap*> elements, bool updateTemp){

    for (auto* elem : elements){  // Loop through element sets

        const vector<vector<int>>& elemConDof_ptr = elem->get_elemConDof();
        nElConDofs = elem->get_nElConDofs();
        nElements = elem->get_nElements();

        //  Assemble the coefficient matrix.        
        PetscInt* i1; PetscInt* j1; // Indices for row and columns to insert.
        PetscMalloc1(nElConDofs, &i1);
        PetscMalloc1(nElConDofs, &j1);

        const T_ElStiffMatx& T_elStiffMatx_ref = elem->getElStiffMatx();

        if (!updateTemp){

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

        } else {

            if (std::holds_alternative<vector<Matd4x4>*>(T_elStiffMatx_ref)){  // Quad4 elements.
    
                const vector<Matd4x4>& elStiffMatx_ref = *std::get<vector<Matd4x4>*>(T_elStiffMatx_ref);

                for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                    // Get the con dofs associated with the element
                    for(int iElDof=0; iElDof<nElConDofs; iElDof++){

                        i1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
                        j1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);

                    }

                    MatSetValues(K, nElConDofs, i1, nElConDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
                }

            } else if (std::holds_alternative<vector<Matd3x3>*>(T_elStiffMatx_ref)){  // Tri3 elements.
    
                const vector<Matd3x3>& elStiffMatx_ref = *std::get<vector<Matd3x3>*>(T_elStiffMatx_ref);

                for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                    // Get the con dofs associated with the element
                    for(int iElDof=0; iElDof<nElConDofs; iElDof++){

                        i1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
                        j1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
                    }
                    
                    MatSetValues(K, nElConDofs, i1, nElConDofs, j1, elStiffMatx_ref.at(iElem).data(), ADD_VALUES);
                }
            } // Element type
        }

        // Free memory
        PetscFree(i1);
        PetscFree(j1);
    }

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    // For Dirichlet boundary conditions
    MatZeroRows(K, nPresDofs, presDofs, 1.0, NULL, NULL);

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

    if (!updateTemp){
        MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    }
}

void TrappingModel::ReadInitialCon(H5IO& H5File, const int iStep){

    vector<double> con_0(nTotDofs);
    H5File.ReadField1D("Con/Step_"+to_string(iStep), con_0);

    for(int iDof=0; iDof<nTotDofs; iDof++){
        VecSetValue(x, iDof, con_0.at(iDof), INSERT_VALUES);
    }
    VecAssemblyBegin(x);  VecAssemblyEnd(x);
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
        H5File_in.ReadField1D(dsetName, dummy);
        // Assign values
        presDofs[iPresDof] = dummy.at(0); // nDim*iNode+dof
        presVals[iPresDof] = dummy.at(1);
        // // TODO: For debug!
        // cout << presDofs[iPresDof] << " --> " << presVals[iPresDof] << "\n";
    }

}

void TrappingModel::setDirichBC(){

    Update_F();
    
    VecSetValues(F, nPresDofs, presDofs, presVals, INSERT_VALUES); 
    VecAssemblyBegin(F); VecAssemblyEnd(F);

    // // TODO: For debug!
    // VecView(F, PETSC_VIEWER_STDOUT_WORLD);
}

void TrappingModel::Update_F(){

    MatMult(M, x, F);
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

void TrappingModel::WriteInPtCoords(vector<BaseElemTrap*> elements, H5IO &H5File_out){

    T_nodStres glIntPtCoords; // global integration points
    glIntPtCoords = vector<ColVecd3>(nTotGuasPts);

    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->getInPtCoords(glIntPtCoords);
    }

    H5File_out.WriteTensor("IntPtCoords", nTotGuasPts, 3, glIntPtCoords);
}

void TrappingModel::CalcFlux(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats){

    // set zeros
    setZero_nodFlux();

    VecGetArrayRead(x, &globalBuffer);

    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->CalcFlux(mats[iSet], globalBuffer, nodFlux, intPtFlux, nodCount, T);
    }
    
    // Number averaging the nodal values
    if (nDim==2){

        for(int iNod=0; iNod<nTotNodes; iNod++){
            
            std::get<std::vector<ColVecd2>>(nodFlux).at(iNod) = std::get<std::vector<ColVecd2>>(nodFlux).at(iNod)/nodCount.at(iNod);
        }

    } else if (nDim==3){

        for(int iNod=0; iNod<nTotNodes; iNod++){
            
            std::get<std::vector<ColVecd3>>(nodFlux).at(iNod) = std::get<std::vector<ColVecd3>>(nodFlux).at(iNod)/nodCount.at(iNod);
        }
    }

    // TODO: For debug!
    // cout << std::get<std::vector<ColVecd3>>(nodStran).at(0) << "\n";

    VecRestoreArrayRead(x, &globalBuffer);
}

void TrappingModel::WriteFlux(H5IO &H5File_out, const string iStep){

    // Write to H5 file
    if (nDim==2){
        H5File_out.WriteTensor("Flux/Step_"+iStep, nTotNodes, 2, nodFlux);
    } else if (nDim==3) {
        H5File_out.WriteTensor("Flux/Step_"+iStep, nTotNodes, 3, nodFlux);
    }
}

void TrappingModel::WriteIntPtFlux(H5IO &H5File_out, const string iStep){

    H5File_out.WriteTensor("IntPtFlux/Step_"+iStep, nTotGuasPts, 3, intPtFlux);
}

void TrappingModel::WriteOut(H5IO &H5File_out, const string iStep){

    // Concentration
    VecGetArrayRead(x, &globalBuffer);
    H5File_out.WriteArray1D("Con/Step_"+iStep, nTotDofs, globalBuffer);
    VecRestoreArrayRead(x, &globalBuffer);
}

void TrappingModel::WriteAvCon(vector<BaseElemTrap*> elements, H5IO &H5File_out, const int iStep){

    double AvCon = 0;

    for (int iSet=0; iSet<nElementSets; iSet++){
        VecGetArrayRead(x, &globalBuffer);
        AvCon += elements[iSet]->CalcAvCon(globalBuffer);
        VecRestoreArrayRead(x, &globalBuffer);
    }

    TotTime += dt;

    H5File_out.WriteScalar("Time/Step_"+to_string(iStep), TotTime);
    H5File_out.WriteScalar("AvCon/Step_"+to_string(iStep), AvCon);
}

void TrappingModel::WriteExitFlux(H5IO &H5File_out, const int iStep){

    double sum = 0;

    for (int iExNod : ExitNodeIDs){ 
        sum += std::get<std::vector<ColVecd2>>(nodFlux).at(iExNod)[0];  // x-component
    }

    sum = sum/nExitNodes;  // Number averaging

    H5File_out.WriteScalar("AvFlux/Step_"+to_string(iStep), sum);
}