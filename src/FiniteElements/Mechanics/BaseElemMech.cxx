#include "FiniteElements/Mechanics/BaseElemMech.h"

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

void BaseElemMech::ReadElementsData(H5IO &H5File_in, int iSet){

    string dsetName;
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nElements";
    nElements = H5File_in.ReadScalar(dsetName);
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nNodes";
    nNodes = H5File_in.ReadScalar(dsetName);

    // Initialize the size.
    elemNodeConn.resize(nElements);  
    elemDispDof.resize(nElements);
    elemIDs.resize(nElements);
    // Read global element IDs.
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/ElementSetIDs";
    H5File_in.ReadFieldInt1D(dsetName, elemIDs);

    // Read node connectivity.
    vector<int> dummy(nElNodes);

    for (int iElem=0; iElem<nElements; iElem++){

        dsetName = "NodeConnectivity/Element_"+to_string(elemIDs.at(iElem));
        H5File_in.ReadFieldInt1D(dsetName, dummy);

        elemNodeConn.at(iElem) = dummy;
        elemDispDof.at(iElem) = CalcElemDispDof(iElem);
    }

    // TODO: For debug!
    // for (auto& s : elemNodeConn[0])
    //     cout << s << "\n"; 
}

vector<int> BaseElemMech::CalcElemDispDof(int iElem){

    vector<int> dispDof(nElDispDofs);

    if (nElDim==2){

        for(int iNod=0; iNod<nElNodes; iNod++){

            dispDof.at(nElDim*iNod) = nElDim*elemNodeConn.at(iElem).at(iNod);
            dispDof.at(nElDim*iNod+1) = nElDim*elemNodeConn.at(iElem).at(iNod)+1;
        }

    } else if (nElDim==3){

        for(int iNod=0; iNod<nElNodes; iNod++){

            dispDof.at(nElDim*iNod) = nElDim*elemNodeConn.at(iElem).at(iNod);
            dispDof.at(nElDim*iNod+1) = nElDim*elemNodeConn.at(iElem).at(iNod)+1;
            dispDof.at(nElDim*iNod+2) = nElDim*elemNodeConn.at(iElem).at(iNod)+2;

        }
    }

    return dispDof;
}