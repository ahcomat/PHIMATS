#include<iostream>
#include<variant>
#include <unordered_set>

#include "Models/TrappingModel.h"

using namespace std;

TrappingModel::TrappingModel(vector<BaseElemTrap*> elements, H5IO& H5File_in, Logger& logger)
              :logger(logger) {

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
    dsetName = "SimulationParameters/conB";
    conB= H5File_in.ReadScalar(dsetName);
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
    VecDestroy(&vecF); VecDestroy(&vecx); MatDestroy(&matK); MatDestroy(&matM);
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
    VecCreate(PETSC_COMM_WORLD, &vecF);
    VecSetSizes(vecF, PETSC_DECIDE, nTotDofs);

    VecSetType(vecF, VECSEQ);
    VecDuplicate(vecF, &vecx);      

    VecSet(vecF, 0.0); // Set all values to zero.

    VecSet(vecx, 0.0); // Set all values to zero.

    // Initialize the coefficient matrices.
    MatCreateSeqAIJ(PETSC_COMM_WORLD, nTotDofs, nTotDofs, PETSC_DEFAULT, NULL, &matK); // Works better for HPC

    // MatCreate(PETSC_COMM_WORLD, &matK);
    // MatSetSizes(matK, PETSC_DECIDE, PETSC_DECIDE, nTotDofs, nTotDofs);
    // MatSetType(matK, MATSEQAIJ); 

    MatDuplicate(matK, MAT_DO_NOT_COPY_VALUES, &matM);

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
    MatSeqAIJSetPreallocation(matK, PETSC_DEFAULT, nnz);
    // No new memory is allocated
    MatSetOption(matK, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    MatSetOption(matK,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

    MatSeqAIJSetPreallocation(matM, PETSC_DEFAULT, nnz); 
    // No new memory is allocated
    MatSetOption(matM, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    MatSetOption(matM,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

    PetscFree(nnz);
}

void TrappingModel::ReadInitialCon(H5IO& H5File, const int iStep){

    vector<double> con_0(nTotDofs);
    H5File.ReadField1D("Con/Step_"+to_string(iStep), con_0);

    for(int iDof=0; iDof<nTotDofs; iDof++){
        VecSetValue(vecx, iDof, con_0.at(iDof), INSERT_VALUES);
    }
    VecAssemblyBegin(vecx);  VecAssemblyEnd(vecx);
}

void TrappingModel::setUniformCon(double uniformCon){

    VecSet(vecx, uniformCon); // Set all values to zero.

}

void TrappingModel::InitializeBC(H5IO& H5File_in){

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

void TrappingModel::CalcEquilibriumBC(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats){

    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->CalcEquilibriumBC(mats[iSet], presVals, presDofs, nPresDofs, conB, T);
    }

}

void TrappingModel::setBC(){

    Update_F();
    
    if (nPresDofs){ // If Dirichlet BC
        VecSetValues(vecF, nPresDofs, presDofs, presVals, INSERT_VALUES); 
        VecAssemblyBegin(vecF); VecAssemblyEnd(vecF);
    }

    // // TODO: For debug!
    // VecView(vecF, PETSC_VIEWER_STDOUT_WORLD);
}

void TrappingModel::ReadNodalStress(vector<BaseElemTrap*> elements, H5IO &H5File_stress, int iStep){

    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->ReadNodalStress(H5File_stress, iStep);
    }

}

// void TrappingModel::WriteGradPhi(vector<BaseElemTrap*> elements, H5IO& H5File_out){

//     // Nodal laplacian of phi
//     double* nodLapPhi = new double[nTotNodes];
//     for (int iNod=0; iNod<nTotDofs; iNod++){
//         nodLapPhi[iNod] = 0;}
//     // Nodal gradient of phi.
//     T_nodStres nodGradPhi;     
//     // Counter for integration points surrounding nodes.
//     vector<double> nodPhiCount(nTotNodes);     
     
//     if (nDim == 2){ // Case 2D model
//         nodGradPhi = vector<ColVecd2>(nTotNodes);
//     } else if (nDim == 3){ // Case 3D
//         nodGradPhi = vector<ColVecd3>(nTotNodes);
//     }

//     // Set zeros
//     if (nDim == 2){ // Case 2D model
//         for(int iNod=0; iNod<nTotNodes; iNod++){
//             std::get<std::vector<ColVecd2>>(nodGradPhi).at(iNod).setZero();
//             nodPhiCount.at(iNod) = 0;
//         }
//     } else if (nDim == 3){ // Case 3D
//         for(int iNod=0; iNod<nTotNodes; iNod++){
//             std::get<std::vector<ColVecd3>>(nodGradPhi).at(iNod).setZero();
//             nodPhiCount.at(iNod) = 0;
//         }
//     }

//     elements[0]->CalcGrad(nodGradPhi, nodPhiCount, nodLapPhi);

//     // Number averaging the nodal values
//     if (nDim==2){
//         for(int iNod=0; iNod<nTotNodes; iNod++){
//             std::get<std::vector<ColVecd2>>(nodGradPhi).at(iNod) = std::get<std::vector<ColVecd2>>(nodGradPhi).at(iNod)/nodPhiCount.at(iNod);
//         }
//     } else if (nDim==3){
//         for(int iNod=0; iNod<nTotNodes; iNod++){
//             std::get<std::vector<ColVecd3>>(nodGradPhi).at(iNod) = std::get<std::vector<ColVecd3>>(nodGradPhi).at(iNod)/nodPhiCount.at(iNod);
//         }
//     }

//     // Flux
//     if (nDim==2){
//         H5File_out.WriteTensor("GradPhi", nTotNodes, 2, nodGradPhi);
//     } else if (nDim==3) {
//         H5File_out.WriteTensor("GradPhi", nTotNodes, 3, nodGradPhi);
//     }

//     // laplacian
//     H5File_out.WriteArray1D("LapPhi", nTotDofs, nodLapPhi);
//     delete [] nodLapPhi;
// }

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

void TrappingModel::CalcElemStiffMatx(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats, std::optional<std::vector<BaseElemPFF*>> pffElemsOpt){

    for (int iSet=0; iSet<nElementSets; iSet++){
        // Check if PFF simulation
        if (pffElemsOpt.has_value()){

            const auto& pffElems = pffElemsOpt.value();
            const std::vector<std::vector<double>>& elPhi_ptr = pffElems[iSet]->getElphi();
            elements[iSet]->CalcElemStiffMatx(mats[iSet], T, &elPhi_ptr);
            
        } else {

            elements[iSet]->CalcElemStiffMatx(mats[iSet], T);

        }
    }
}

void TrappingModel::AssembleElementMatrix(const auto* elMatx_ptr,
                                          const vector<vector<int>>& elemConDof_ptr,
                                          PetscInt nElConDofs,
                                          PetscInt nElements,
                                          Mat globalMat) {
    PetscInt* i1;
    PetscInt* j1;

    // Allocate memory for indices
    PetscMalloc1(nElConDofs, &i1);
    PetscMalloc1(nElConDofs, &j1);

    for (int iElem = 0; iElem < nElements; iElem++) {
        // Map local DOFs to global DOFs
        for (int iElDof = 0; iElDof < nElConDofs; iElDof++) {
            i1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
            j1[iElDof] = elemConDof_ptr.at(iElem).at(iElDof);
        }

        // Add matrix values
        MatSetValues(globalMat, nElConDofs, i1, nElConDofs, j1, elMatx_ptr->at(iElem).data(), ADD_VALUES);
    }

    // Free memory
    (void)PetscFree(i1);
    (void)PetscFree(j1);
}

void TrappingModel::Assemble(vector<BaseElemTrap*> elements, bool assembleM){

    // Has to be zeroed for iterative solver. 
    MatZeroEntries(matK);
    MatZeroEntries(matM);

    for (auto* elem : elements) {

        const vector<vector<int>>& elemConDof_ptr = elem->get_elemConDof();
        nElConDofs = elem->get_nElConDofs();
        nElements = elem->get_nElements();

        const T_ElStiffMatx& T_elStiffMatx_ref = elem->getElStiffMatx();
        const T_ElStiffMatx* T_elCapMatx_ref = nullptr;

        if (assembleM) {
            T_elCapMatx_ref = &elem->getElCapMatx();
        }

        auto assembleMatrix = [&](const auto* elMatx_ptr, Mat globalMat) {
            try{
                if (elMatx_ptr) {
                    AssembleElementMatrix(elMatx_ptr, elemConDof_ptr, nElConDofs, nElements, globalMat);
                } else {
                    logger.log("Unsupported element type in T_elStiffMatx_ref", "ERROR");
                    throw std::runtime_error("Unsupported element type in T_elStiffMatx_ref");
                }
            } catch (const std::runtime_error& e) {
                logger.log("\nException caught in TrappingModel::Assemble:\n", "", false);
                logger.log("    " + std::string(e.what()), "", false);
                logger.log("\nCritical error encountered. Terminating!\n", "", false);
                exit(EXIT_FAILURE);
            }
        };

        // Use std::visit to handle the variant
        std::visit(
            [&](auto&& elStiffMatx_ptr) { assembleMatrix(elStiffMatx_ptr, matK); },
            T_elStiffMatx_ref);

        // Assemble capacity matrix (if applicable)
        if (assembleM) {
            std::visit(
                [&](auto&& elCapMatx_ptr) { assembleMatrix(elCapMatx_ptr, matM); },
                *T_elCapMatx_ref);
        }
    }

    // Finalize stiffness matrix assembly
    MatAssemblyBegin(matK, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(matK, MAT_FINAL_ASSEMBLY);

    if (nPresDofs){ // If Dirichlet BC
        MatZeroRows(matK, nPresDofs, presDofs, 1.0, NULL, NULL);
    }
    
    if (assembleM) { // Only on intitializaiton 
        // Finalize capacity matrix assembly
        MatAssemblyBegin(matM, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(matM, MAT_FINAL_ASSEMBLY);
    }

    // TODO: For debug!
    // // Extract a specific row
    // PetscInt ROW = 4676; // for example, extracting row 2
    // const PetscInt *COL;
    // const PetscScalar *VAL;
    // PetscInt NCOL;

    // MatGetRow(matK, ROW, &NCOL, &COL, &VAL);

    // for (int w = 0; w < NCOL; w++) {
    //     cout << "Col: " << COL[w] << ", Val: " << VAL[w] << "\n";
    // }

    // // Restore the row
    // MatRestoreRow(matK, ROW, &NCOL, &COL, &VAL); 

    // // Save the matrix to MatrixMarket format
    // PetscViewer viewer;
    // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix_output.mtx", &viewer);
    // PetscViewerPushFormat(viewer,  PETSC_VIEWER_ASCII_MATRIXMARKET);
    // PetscViewerPushFormat(viewer,  PETSC_VIEWER_ASCII_MATLAB);

    // MatView(matK, viewer);
    // PetscViewerDestroy(&viewer);
}

void TrappingModel::Update_F(){

    MatMult(matM, vecx, vecF);
}

int TrappingModel::get_nSteps() const{
    
    return nSteps;
}

Vec& TrappingModel::getF(){

    return vecF;
}

Vec& TrappingModel::getX(){

    return vecx;
}

Mat& TrappingModel::getK(){

    return matK;
}

void TrappingModel::WriteInPtCoords(vector<BaseElemTrap*> elements, H5IO &H5File_out){

    T_nodStres glIntPtCoords; // global integration points
    glIntPtCoords = vector<ColVecd3>(nTotGuasPts);

    for (int iSet=0; iSet<nElementSets; iSet++){
        elements[iSet]->getInPtCoords(glIntPtCoords);
    }

    H5File_out.WriteTensor("IntPtCoords", nTotGuasPts, 3, glIntPtCoords);
}

void TrappingModel::CalcFlux(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats, std::optional<std::vector<BaseElemPFF*>> pffElemsOpt){

    // set zeros
    setZero_nodFlux();

    VecGetArrayRead(vecx, &globalBuffer);

    for (int iSet=0; iSet<nElementSets; iSet++){
        // Check if PFF simulation
        if (pffElemsOpt.has_value()){

            const auto& pffElems = pffElemsOpt.value();

            elements[iSet]->CalcFlux(mats[iSet], globalBuffer, nodFlux, intPtFlux, nodCount, T,  &pffElems[iSet]->getEl_gPhi_d());
            
        } else {

            elements[iSet]->CalcFlux(mats[iSet], globalBuffer, nodFlux, intPtFlux, nodCount, T);
        }
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

    VecRestoreArrayRead(vecx, &globalBuffer);
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
    VecGetArrayRead(vecx, &globalBuffer);
    H5File_out.WriteArray1D("Con/Step_"+iStep, nTotDofs, globalBuffer);
    VecRestoreArrayRead(vecx, &globalBuffer);
}

void TrappingModel::WriteAvCon(vector<BaseElemTrap*> elements, H5IO &H5File_out, const int iStep){

    double AvCon = 0;

    for (int iSet=0; iSet<nElementSets; iSet++){
        VecGetArrayRead(vecx, &globalBuffer);
        AvCon += elements[iSet]->CalcAvCon(globalBuffer);
        VecRestoreArrayRead(vecx, &globalBuffer);
    }

    TotTime += dt;

    H5File_out.WriteScalar("Time/Step_"+to_string(iStep), TotTime);
    H5File_out.WriteScalar("AvCon/Step_"+to_string(iStep), AvCon);
}

void TrappingModel::WriteExitFlux(H5IO &H5File_out, const int iStep){

    double sum = 0;

    for (int iExNod : ExitNodeIDs){ 
        sum += std::get<std::vector<ColVecd2>>(nodFlux).at(iExNod)[0];  // vecx-component
    }

    sum = sum/nExitNodes;  // Number averaging

    H5File_out.WriteScalar("AvFlux/Step_"+to_string(iStep), sum);
}