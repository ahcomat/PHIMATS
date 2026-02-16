#include<iostream>
#include<variant>
#include <unordered_set>

#include "Models/PFFModel.h"

using namespace std;

PFFModel::PFFModel(vector<BaseElemPFF*> elements, H5IO& H5File_in, Logger& logger)
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

    // Allocate memory for `Fint` and `indices`.
    PetscMalloc1(nTotDofs, &FH);
    PetscMalloc1(nTotDofs, &indices);

    // Initialize to zero.
    for (int iDof=0; iDof<nTotDofs; iDof++){
        FH[iDof] = 0;
    }

    // Set indices
    for (int iDof=0; iDof<nTotDofs; iDof++){
        indices[iDof] = iDof;
    }

    InitializePETSc(elements);
}

PFFModel::~PFFModel(){

    // Deallocate memory
    PetscFree(FH); PetscFree(indices);
    VecDestroy(&vecFp); VecDestroy(&vecx); MatDestroy(&matK);
}

void PFFModel::InitializePETSc(vector<BaseElemPFF*> elements){

    // TODO: For debug!
    // for(auto* elem : elements)
    //     for (int dispDof : elem->get_elemDispDof(0))
    //         cout << dispDof << "\n";

    // Initialize the vectors
    VecCreate(PETSC_COMM_WORLD, &vecFp);
    VecSetSizes(vecFp, PETSC_DECIDE, nTotDofs);

    VecSetType(vecFp, VECSEQ);
    VecDuplicate(vecFp, &vecx);      

    VecSet(vecFp, 0.0); // Set all values to zero.

    VecSet(vecx, 0.0); // Set all values to zero.

    // Initialize the coefficient matrices.
    MatCreateSeqAIJ(PETSC_COMM_WORLD, nTotDofs, nTotDofs, PETSC_DEFAULT, NULL, &matK); // Works better for HPC

    // MatCreate(PETSC_COMM_WORLD, &matK);
    // MatSetSizes(matK, PETSC_DECIDE, PETSC_DECIDE, nTotDofs, nTotDofs);
    // MatSetType(matK, MATSEQAIJ); 

    // Preallocate the coefficient matrix.
    vector<unordered_set<int>> gDofs(nTotDofs); // vector to store dofs per row.
    int itotv, jtotv; // for global row and colum dof.
    PetscInt *nnz; // Array for the number of zeros per row
    PetscMalloc1(nTotDofs, &nnz); // Allocates the size of nnz
    
    // Find the number of non zeros (columns) per row for preallocation.
    for (auto* elem : elements){  // Loop through element sets

        nElPhiDofs = elem->get_nElPhiDofs();
        nElements = elem->get_nElements();
        const vector<vector<int>>& elemPhiDof_ptr = elem->get_elemPhiDof();
    
        for (int iElem=0; iElem<nElements; iElem++){ // Loop through all elements per element set

            for (int idof=0; idof<nElPhiDofs; idof++){ // Row

                itotv = elemPhiDof_ptr.at(iElem).at(idof); // global row number
                
                for (int jdof=0; jdof<nElPhiDofs; jdof++){// Column.

                    jtotv = elemPhiDof_ptr.at(iElem).at(jdof); // global column number

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
    MatSetOption(matK, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    PetscFree(nnz);

}

void PFFModel::set_const_wc(vector<BaseElemPFF*> pffElem, vector<BasePFF*> pffMat){
    
    try{
        
        for (int iSet=0; iSet<nElementSets; iSet++){

            pffElem[iSet]->set_const_wc(pffMat[iSet]->get_wc());

        }

    } catch (const exception& e) {

        logger.log("Exception caught in PFFModel::set_const_wc:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("    Critical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

    }
}

void  PFFModel::set_decay_wc(vector<BaseElemPFF*> pffElem, vector<BasePFF*> pffMat, vector<BaseElemTrap*> trapElem){
    
    try{
        
        for (int iSet=0; iSet<nElementSets; iSet++){

            // Cast BasePFF material pointer to ChemoMech class
            auto* chemoMat = dynamic_cast<ChemoMech*>(pffMat[iSet]);
            
            if (chemoMat) {
                // Retrieve int-pt concentration field
                const std::vector<std::vector<double>>& el_con = trapElem[iSet]->getElCon();

                pffElem[iSet]->set_decay_wc(
                    el_con, 
                    chemoMat->get_wc(), 
                    chemoMat->get_wc_min(), 
                    chemoMat->get_beta(),
                    chemoMat->get_c_crit()
                );
            } else {
                logger.log("Material at set " + to_string(iSet) + " is not of type ChemoMech", "ERROR", true);
            }
            
        }

    } catch (const exception& e) {

        logger.log("Exception caught in PFFModel::set_decay_wc:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("    Critical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

    }
}

void PFFModel::CalcPsiSpectral(vector<BaseElemPFF*> pffElem, vector<BaseElemMech*> mechElem, vector<BaseMechanics*> mechMat){

    try{

        for (int iSet=0; iSet<nElementSets; iSet++){

            const T_elStres& T_elStrain_ref = mechElem[iSet]->getElStrain_e();

            double Gmod = mechMat[iSet]->getGmod();
            double Lambda = mechMat[iSet]->getLambda();

            pffElem[iSet]->CalcPsiSpectral(Lambda, Gmod, T_elStrain_ref);

        }  

    } catch (const exception& e) {

        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);

    }
}

void PFFModel::CalcDrivForcE(vector<BaseElemPFF*> pffElem){

    for (int iSet=0; iSet<nElementSets; iSet++){

        pffElem[iSet]->CalcDrivForcE();

    }  
}

void PFFModel::CalcDrivForcEP(vector<BaseElemPFF*> pffElem, vector<BaseElemMech*> mechElem){

    for (int iSet=0; iSet<nElementSets; iSet++){

        const std::vector<std::vector<double>>& el_wp = mechElem[iSet]->getEl_wp();

        pffElem[iSet]->CalcDrivForcEP(&el_wp);

    }
}

void PFFModel::CalcDrivForcEP_TH(vector<BaseElemPFF*> pffElem, vector<BaseElemMech*> mechElem, const double zeta){

    for (int iSet=0; iSet<nElementSets; iSet++){

        const std::vector<std::vector<double>>& el_wp = mechElem[iSet]->getEl_wp();

        pffElem[iSet]->CalcDrivForcEP_TH(&el_wp, zeta);

    } 
}

void PFFModel::CalcElemStiffMatx(vector<BaseElemPFF*> elements){

        for (int iSet=0; iSet<nElementSets; iSet++){
            elements[iSet]->CalcElemStiffMatx();
        }
}

void PFFModel::AssembleElementMatrix(const auto* elMatx_ptr,
                                          const vector<vector<int>>& elemPhiDof_ptr,
                                          PetscInt nElPhiDofs,
                                          PetscInt nElements,
                                          Mat globalMat) {
    PetscInt* i1;
    PetscInt* j1;

    // Allocate memory for indices
    PetscMalloc1(nElPhiDofs, &i1);
    PetscMalloc1(nElPhiDofs, &j1);

    for (int iElem = 0; iElem < nElements; iElem++) {
        // Map local DOFs to global DOFs
        for (int iElDof = 0; iElDof < nElPhiDofs; iElDof++) {
            i1[iElDof] = elemPhiDof_ptr.at(iElem).at(iElDof);
            j1[iElDof] = elemPhiDof_ptr.at(iElem).at(iElDof);
        }

        // Add matrix values
        MatSetValues(globalMat, nElPhiDofs, i1, nElPhiDofs, j1, elMatx_ptr->at(iElem).data(), ADD_VALUES);
    }

    // Free memory
    (void)PetscFree(i1);
    (void)PetscFree(j1);
}

void PFFModel::Assemble(vector<BaseElemPFF*> pffElem){

    // Has to be zeroed for iterative solver. 
    MatZeroEntries(matK);

    for (auto* elem : pffElem) {

        const vector<vector<int>>& elemPhiDof_ptr = elem->get_elemPhiDof();
        nElPhiDofs = elem->get_nElPhiDofs();
        nElements = elem->get_nElements();

        const T_ElStiffMatx& T_elStiffMatx_ref = elem->getElStiffMatx();

        auto assembleMatrix = [&](const auto* elMatx_ptr, Mat globalMat) {
            try{
                if (elMatx_ptr) {
                    AssembleElementMatrix(elMatx_ptr, elemPhiDof_ptr, nElPhiDofs, nElements, globalMat);
                } else {
                    logger.log("Unsupported element type in T_elStiffMatx_ref", "ERROR");
                    throw std::runtime_error("Unsupported element type in T_elStiffMatx_ref");
                }
            } catch (const std::runtime_error& e) {
                logger.log("\nException caught in PFFModel::Assemble:\n", "", false);
                logger.log("    " + std::string(e.what()), "", false);
                logger.log("\nCritical error encountered. Terminating!\n", "", false);
                exit(EXIT_FAILURE);
            }
        };

        // Use std::visit to handle the variant
        std::visit(
            [&](auto&& elStiffMatx_ptr) { assembleMatrix(elStiffMatx_ptr, matK); },
            T_elStiffMatx_ref);
    }

    // Finalize stiffness matrix assembly
    MatAssemblyBegin(matK, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(matK, MAT_FINAL_ASSEMBLY);

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

void PFFModel::CalcFH(vector<BaseElemPFF*> pffElem){


    // Has to be set to zero before any calculation, otherwise it accumulates.
    for (int iDof=0; iDof<nTotDofs; iDof++){
        FH[iDof] = 0;
    }

    for (int iSet=0; iSet<nElementSets; iSet++){
    
        // Calculate internal force
        pffElem[iSet]->CalcFH(FH);
        VecSetValues(vecFp, nTotDofs, indices, FH, INSERT_VALUES); 
        VecAssemblyBegin(vecFp); VecAssemblyEnd(vecFp);
    }

}

void PFFModel::Calc_gPhi_d(vector<BaseElemPFF*> pffElem){

    VecGetArrayRead(vecx, &globalBuffer);
    for (int iSet=0; iSet<nElementSets; iSet++){
        pffElem[iSet]->Calc_gPhi_d(globalBuffer);
    }
    VecRestoreArrayRead(vecx, &globalBuffer);
}


Vec& PFFModel::getF(){

    return vecFp;
}

Vec& PFFModel::getX(){

    return vecx;
}

Mat& PFFModel::getK(){

    return matK;
}

void PFFModel::WriteOut(H5IO &H5File_out, const string iStep){

    // Concentration
    VecGetArrayRead(vecx, &globalBuffer);
    H5File_out.WriteArray1D("Phi/Step_"+iStep, nTotDofs, globalBuffer);
    VecRestoreArrayRead(vecx, &globalBuffer);
}