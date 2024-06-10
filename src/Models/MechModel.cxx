#include<iostream>
#include<variant>

#include"Models/MechModel.h"

using namespace std;

MechModel::MechModel(vector<BaseElemMech*> elements, H5IO& H5File_in){

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

    // Initialize to zeros
    setZero_nodStres();

    // Allocate memory for `Fint`.
    PetscMalloc1(nTotDofs, &Fint);
    // Initialize to zeros, otherwise will get garbage memory values.
    for(int iDof=0; iDof<nTotDofs; iDof++){
        Fint[iDof] = 0;
    }

    InitializePETSc(elements);
}

MechModel::~MechModel(){

    // Deallocate memory.
    PetscFree(presDofs); PetscFree(presVals); PetscFree(Fint);
    VecDestroy(&b); VecDestroy(&x); MatDestroy(&A);
    // Exit message
    cout << "MechModel elements exited correctly" << "\n";
}

void MechModel::setZero_nodStres(){

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
    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b, PETSC_DECIDE, nTotDofs);

    // Since we are interested in sequential implementation for now.
    VecSetType(b, VECSEQ);
    // VecSetFromOptions(b); =// for the general case.
    VecDuplicate(b, &x);      // Initialize the solution vector

    VecSet(b, 0.0); // Set all values to zero.
    VecAssemblyBegin(b); VecAssemblyEnd(b);

    // Initialize the coefficient matrix.
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nTotDofs, nTotDofs);
    
    // Since we are interested in only sequential in this implementation. 
    MatSetType(A, MATSEQAIJ);
    // MatSetFromOptions(A);  // for command line options, but we dont do it here.

    // MAT_SYMMETRIC: symmetric in terms of both structure and value
    MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);

    // Preallocate the coefficient matrix.
    vector<vector<int>> gDofs(nTotDofs); // vector to store dofs per row.
    int itotv, jtotv; // for global row and colum dof.
    PetscInt *nnz; // Array for the number of zeros per row
    PetscMalloc1(nTotDofs, &nnz); // Allocates the size of nnz
    
    // Find the number of non zeros (columns) per row for preallocation.
    // TODO: Make a check for elements.size()==nElementSets?
    for (auto* elem : elements){  // Loop through element sets

        nElDispDofs = elem->get_nElDispDofs();
        nElements = elem->get_nElements();
        const vector<vector<int>>& elemDispDof_ptr = elem->get_elemDispDof();
    
        for (int iElem=0; iElem<elem->get_nElements(); iElem++){ // Loop through all elements per element set

            for (int idof=0; idof<nElDispDofs; idof++){ // Row

                itotv = elemDispDof_ptr.at(iElem).at(idof); // global row number
                
                for (int jdof=0; jdof<nElDispDofs; jdof++){// Column.

                    jtotv = elemDispDof_ptr.at(iElem).at(jdof); // global column number

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
    MatSeqAIJSetPreallocation(A, PETSC_DEFAULT, nnz); 
    PetscFree(nnz); 

    // Throw error if unallocated entry is accessed.
    MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE); 
}

void MechModel::CalcElemStiffMatx(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats){

    for (int iSet=0; iSet<nElementSets; iSet++)
      elements[iSet]->CalcElemStiffMatx(mats[iSet]->getDMatx());
}

void MechModel::Assemble(vector<BaseElemMech*> elements){

    MatZeroEntries(A);  // Set all entries to zero.

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

        PetscScalar* vals;   // values.
        PetscMalloc1(nElDispDofs*nElDispDofs, &vals);

        // `T_ElStiffMatx` variant that holds a pointer to the vector.
        const T_ElStiffMatx& elStiffMatx_ptr = elem->getElStiffMatx();

        if (std::holds_alternative<vector<Matd8x8>*>(elStiffMatx_ptr)){  // Quad4 elements.

            /* The "*" operator dereferences the pointer to get the actual variable that the 
            pointer is pointing to. const auto& binds the dereferenced vector to a constant reference. 
            This way avoids copying the vector and can safely read from it without modifying it.
            */ 
            const vector<Matd8x8> elStiffMatx_ref = *std::get<vector<Matd8x8>*>(elStiffMatx_ptr);

            Matd8x8 dummyVals;

            for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                // Get the disp dofs associated with the element
                for(int iElDof=0; iElDof<nElDispDofs; iElDof++){

                    i1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                    j1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                }

                // Get the element stiffness matrix
                dummyVals = elStiffMatx_ref.at(iElem);

                Matd8x8::Map(vals, dummyVals.rows(), dummyVals.cols()) = dummyVals;
                MatSetValues(A, nElDispDofs, i1, nElDispDofs, j1, vals, ADD_VALUES);
            }
        } else if (std::holds_alternative<vector<Matd6x6>*>(elStiffMatx_ptr)){  // Tri3 elements.
 
            const vector<Matd6x6> elStiffMatx_ref = *std::get<vector<Matd6x6>*>(elStiffMatx_ptr);

            Matd6x6 dummyVals;

            for (int iElem =0; iElem<nElements; iElem++){ // Loop through elements

                // Get the disp dofs associated with the element
                for(int iElDof=0; iElDof<nElDispDofs; iElDof++){

                    i1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                    j1[iElDof] = elemDispDof_ptr.at(iElem).at(iElDof);
                }

                // Get the element stiffness matrix
                dummyVals = elStiffMatx_ref.at(iElem);

                Matd6x6::Map(vals, dummyVals.rows(), dummyVals.cols()) = dummyVals;
                MatSetValues(A, nElDispDofs, i1, nElDispDofs, j1, vals, ADD_VALUES);
            }
        }

        // Free memory
        PetscFree(i1);
        PetscFree(j1);
        PetscFree(vals);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

void MechModel::InitializeDirichBC(H5IO& H5File_in){

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
        presDofs[iPresDof] = nDim*dummy.at(0)+dummy.at(1); // nDim*iNode+dof
        presVals[iPresDof] = dummy.at(2)/nSteps;
        // TODO: For debug!
        // cout << presDofs[iPresDof] << "-->" << presVals[iPresDof] << "\n";
    }
}

void MechModel::setDirichBC(){

    // MatZeroRowsColumns(A, nPresDofs, presDofs, 1.0, NULL, NULL);
    MatZeroRows(A, nPresDofs, presDofs, 1.0, NULL, NULL);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // TODO: Check for incremental loads "ADD_VALUES"
    VecSetValues(b, nPresDofs, presDofs, presVals, ADD_VALUES); 
    VecAssemblyBegin(b); VecAssemblyEnd(b);
}

int MechModel::get_nSteps() const{
    
    return nSteps;
}

Vec& MechModel::getB(){

    return b;
}

Vec& MechModel::getX(){

    return x;
}

Mat& MechModel::getA(){

    return A;
}

void MechModel::CalcStres(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats){

    VecGetArrayRead(x, &globalBuffer);

    for (int iSet=0; iSet<nElementSets; iSet++){
        
        elements[iSet]->CalcStres(mats[iSet]->getDMatx(), globalBuffer, Fint, nodStres, nodStran, nodCount);
    }

    // Number averaging the nodal values
    for(int iNod=0; iNod<nTotNodes; iNod++){
        
        std::get<std::vector<ColVecd3>>(nodStran).at(iNod) = std::get<std::vector<ColVecd3>>(nodStran).at(iNod)/nodCount.at(iNod);
        std::get<std::vector<ColVecd3>>(nodStres).at(iNod) = std::get<std::vector<ColVecd3>>(nodStres).at(iNod)/nodCount.at(iNod);
    }

    // TODO: For debug!
    // cout << std::get<std::vector<ColVecd3>>(nodStran).at(0) << "\n";

    VecRestoreArrayRead(x, &globalBuffer);
}

void MechModel::WriteOut(vector<BaseElemMech*> elements, H5IO &H5File_out, const string iStep){

    // Displacements
    VecGetArrayRead(x, &globalBuffer);
    H5File_out.WriteArray_1D("Disp/Step_"+iStep, nTotDofs, globalBuffer);
    VecRestoreArrayRead(x, &globalBuffer);
    H5File_out.WriteArray_1D("Force/Step_"+iStep, nTotDofs, Fint);

    // Stresses and strains
    if (nDim==2){
        H5File_out.WriteStres("Strain/Step_"+iStep, nTotNodes, 3, nodStran);
        H5File_out.WriteStres("Stress/Step_"+iStep, nTotNodes, 3, nodStres);
    } else if (nDim==3) {
        H5File_out.WriteStres("Strain/Step_"+iStep, nTotNodes, 6, nodStran);
        H5File_out.WriteStres("Stress/Step_"+iStep, nTotNodes, 6, nodStres);
    }

    // set zeros
    setZero_nodStres();
}

