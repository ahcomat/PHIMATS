#include<iostream>
#include"Models/MechModel.h"

using namespace std;

MechModel::MechModel(BaseElemMech* elements){

    nTotDof = elements->get_nTotDof();
    nElements = elements->get_nElements();
    nElDispDofs = elements->get_nElDispDofs();
    nDim = elements->get_nDim();

    InitializePETSc(elements);
}

MechModel::~MechModel(){

    // Deallocate memory.
    PetscFree(presDofs); PetscFree(presVals);
    VecDestroy(&b); VecDestroy(&x); MatDestroy(&A);
    // Exit message
    cout << "MechModel elements exited correctly" << "\n";
}

void MechModel::CalcElemStiffMatx(BaseElemMech* elements, T_DMatx DMatx){

    elements->CalcElemStiffMatx(DMatx);
}

void MechModel::InitializePETSc(BaseElemMech* elements){

    const vector<vector<int>>& elemDispDof = elements->get_elemDispDof();

    // for (int dispDof : elemDispDof[0])
    //     cout << dispDof << "\n";

    // Initialize the vectors
    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b, PETSC_DECIDE, nTotDof);

    // Since we are interested in sequential implementation for now.
    VecSetType(b, VECSEQ);
    // VecSetFromOptions(b); =// for the general case.
    VecDuplicate(b, &x);      // Initialize the solution vector

    VecSet(b, 0.0); // Set all values to zero.
    VecAssemblyBegin(b); VecAssemblyEnd(b);

    // Initialize the coefficient matrix.
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nTotDof, nTotDof);
    
    // Since we are interested in only sequential in this implementation. 
    MatSetType(A, MATSEQAIJ);
    // MatSetFromOptions(A);  // for command line options, but we dont do it here.

    // MAT_SYMMETRIC: symmetric in terms of both structure and value
    MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);

    // Preallocate the coefficient matrix.
    vector<vector<int>> gDofs(nTotDof); // vector to store dofs per row.
    int itotv, jtotv; // for global row and colum dof.
    PetscInt *nnz; // Array for the number of zeros per row
    PetscMalloc1(nTotDof, &nnz); // Allocates the size of nnz
    
    // Find the number of non zeros (columns) per row for preallocation. 
    for (int iElem=0; iElem<nElements; iElem++){
        // Row.
        for (int idof=0; idof<nElDispDofs; idof++){

            itotv = elemDispDof.at(iElem).at(idof); // global row number
            // Column.
            for (int jdof=0; jdof<nElDispDofs; jdof++){

                jtotv = elemDispDof.at(iElem).at(jdof); // global column number

                // If jtotv is not in gDofs.at(itotv)
                if (find(gDofs.at(itotv).begin(), gDofs.at(itotv).end(), jtotv) == gDofs.at(itotv).end()){
                    gDofs.at(itotv).push_back(jtotv);
                }
            }
        }
    }

    // Add the number of non zero columns to nnz
    for (int iDof=0; iDof<nTotDof; iDof++){
        nnz[iDof] = gDofs.at(iDof).size();
    }

    // Preallocate the stiffness matrix.
    MatSeqAIJSetPreallocation(A, PETSC_DEFAULT, nnz); 
    PetscFree(nnz); 

    // Throw error if unallocated entry is accessed.
    MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE); 
}

void MechModel::Assemble(BaseElemMech* elements){

    const vector<vector<int>>& elemDispDof = elements->get_elemDispDof();

    MatZeroEntries(A);  // Set all entries to zero.

    //  Assemble the coefficient matrix.
    PetscInt   i1[nElDispDofs], j1[nElDispDofs]; // Indices for row and columns to insert.
    PetscScalar vals[nElDispDofs*nElDispDofs];   // values.

    // `T_ElStiffMatx` variant that holds a pointer to the vector.
    const T_ElStiffMatx& elStiffMatx_ptr = elements->getElStiffMatx();

    if (std::holds_alternative<vector<Matd8x8>*>(elStiffMatx_ptr)){

        /* The "*" operator dereferences the pointer to get the actual variable that the 
           pointer is pointing to. const auto& binds the dereferenced vector to a constant reference. 
           This way avoids copying the vector and can safely read from it without modifying it.
        */ 
        const vector<Matd8x8> elStiffMatx_ref = *std::get<vector<Matd8x8>*>(elStiffMatx_ptr);

        Matd8x8 dummyVals;

        for (int iElem =0; iElem<nElements; iElem++){
            // Get the disp dofs associated with the element
            for(int iElDof=0; iElDof<nElDispDofs; iElDof++){
                i1[iElDof] = elemDispDof.at(iElem).at(iElDof);
                j1[iElDof] = elemDispDof.at(iElem).at(iElDof);
            }

            // Get the element stiffness matrix
            dummyVals = elStiffMatx_ref.at(iElem);

            Matd8x8::Map(vals, dummyVals.rows(), dummyVals.cols()) = dummyVals;
            MatSetValues(A, nElDispDofs, i1, nElDispDofs, j1, vals, ADD_VALUES);
        }
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
    dsetName = "SimulationParameters/nSteps";
    nSteps = H5File_in.ReadScalar(dsetName);

    vector<double> dummy(3);
    for (int iPresDof=0; iPresDof<nPresDofs; iPresDof++){
        // Read values
        dsetName = "PrescribedDOFs/Prescribed_"+to_string(iPresDof);
        H5File_in.ReadFieldDoub1D(dsetName, dummy);
        // Assign values
        presDofs[iPresDof] = nDim*dummy.at(0)+dummy.at(1); // nDim*iNode+dof
        presVals[iPresDof] = dummy.at(2);
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

void MechModel::CalcStres(BaseElemMech* elements, T_DMatx DMatx){

    VecGetArrayRead(x, &globalBuffer);
    elements->CalcStres(DMatx, globalBuffer);
    VecRestoreArrayRead(x, &globalBuffer);
}

void MechModel::WriteOut(BaseElemMech* elements, H5IO &H5File_out){

    // Displacements
    VecGetArrayRead(x, &globalBuffer);
    H5File_out.WriteArray_1D("Disp", nTotDof, globalBuffer);
    VecRestoreArrayRead(x, &globalBuffer);

    elements->WriteOut(H5File_out);
}

