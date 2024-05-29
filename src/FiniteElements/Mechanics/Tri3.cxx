/**
 * @file Tri3.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing tri3 elements.
 * @date 2024-05-23
 * 
 * @copyright Copyright (c) 2024
 * 
 * Updates (when, what and who)
 * 
 */

#include<iostream>
#include<algorithm>

#include"FiniteElements/Mechanics/Tri3.h"

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

Tri3::Tri3(H5IO &H5File_in, Nodes &Nodes)
    : BaseElemMech(2, 3, 2, 3, 6, 1){ // nDim, nElNodes, dispDofs, nStres, nElDispDofs, nGauss

    InitShapeFunc();
    ReadElementsData(H5File_in);
    // InitializeElements(Nodes);
    // InitPETSC();
}

Tri3::~Tri3(){

    // PetscFree(presDofs); PetscFree(presVals); PetscFree(Fint);
    // VecDestroy(&b); MatDestroy(&A);
    // Exit message
    cout << "Tri3 elements exited correctly" << "\n";
}

void Tri3::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    vector<double> ip = {1.0/3.0, 1.0/3.0};
    gaussPts.push_back(ip);

    // for (auto& s : gaussPts[0])
    //     cout << s << "\n";

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nElNodes);
    shapeFuncDeriv.resize(nElNodes);

    for(int i=0; i<nGauss; i++){
        shapeFunc.at(i) = getShapeFunc(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
        shapeFuncDeriv.at(i) = getShapeFuncDeriv(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
    }

}

RowVecd3 Tri3::getShapeFunc(double xi, double eta){

    // N_i
    RowVecd3 shape;

    shape(0) = xi;
    shape(1) = eta;
    shape(2) = 1 - xi - eta;

    return shape;
}

Matd2x3 Tri3::getShapeFuncDeriv(double xi, double eta){

    // dN_ji
    Matd2x3 shapeDeriv;

    shapeDeriv(0,0) = 1;
    shapeDeriv(0,1) = 0;
    shapeDeriv(0,2) = -1;

    shapeDeriv(1,0) = 0;
    shapeDeriv(1,1) = 1;
    shapeDeriv(1,2) = -1;
             
    return shapeDeriv;
}

void Tri3::ReadElementsData(H5IO &H5File_in){

        string dsetName;
    dsetName = "SimulationParameters/nElements";
    nElements = H5File_in.ReadScalar(dsetName);

    // Initialize the size.
    elemNodeConn.resize(nElements);  
    elemDispDof.resize(nElements);

    dsetName = "SimulationParameters/nElementSets";
    nElementSets = H5File_in.ReadScalar(dsetName);

    // Read node connectivity.
    vector<int> dummy(nElNodes);
    for (int iElem=0; iElem<nElements; iElem++){
        dsetName = "NodeConnectivity/Element_"+to_string(iElem);
        H5File_in.ReadFieldInt1D(dsetName, dummy);

        elemNodeConn.at(iElem) = dummy;
        elemDispDof.at(iElem) = getElemDispDof(iElem);
    }

    for (auto& s : elemNodeConn[0])
        cout << s << "\n";
}

vector<int> Tri3::getElemDispDof(int iElem){

    vector<int> dispDof(nElDispDofs);
    for(int iNod=0; iNod<nElNodes; iNod++){

        dispDof.at(2*iNod) = 2*elemNodeConn.at(iElem).at(iNod);
        dispDof.at(2*iNod+1) = 2*elemNodeConn.at(iElem).at(iNod)+1;
    }

    return dispDof;
}
