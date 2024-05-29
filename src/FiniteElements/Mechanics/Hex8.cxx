/**
 * @file Hex8.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing hex8 elements. 
 * @date 2024-05-23
 * 
 * @copyright Copyright (c) 2024
 * 
 * Updates (when, what and who)
 * 
 */

#include<iostream>

#include"FiniteElements/Mechanics/Hex8.h"

#ifndef DEBUG
#define at(x) operator[](x)
#endif

/*
 Unifying indices:
 i -> number of nodes.
 j -> spatial dimensions.
 k -> number of stresses.
 l -> total displacement dofs.
*/

Hex8::Hex8(H5IO &H5File_in, Nodes &Nodes)
    : BaseElemMech(2, 4, 2, 3, 8, 4){ // nDim, nElNodes, dispDofs, nStres, nElDispDofs, nGauss

    // InitShapeFunc();
    // ReadElementsData(H5File_in);
    // InitializeElements(Nodes);
    // InitPETSC();
}

Hex8::~Hex8(){

    // PetscFree(presDofs); PetscFree(presVals); PetscFree(Fint);
    // VecDestroy(&b); MatDestroy(&A);
    // Exit message
    cout << "Elements exited correctly" << "\n";
}