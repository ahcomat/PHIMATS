#include "FiniteElements/Trapping/BaseElemTrap.h"

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

void BaseElemTrap::ReadElementsData(H5IO &H5File_in, int iSet){

    string dsetName;
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nElements";
    nElements = H5File_in.ReadScalar(dsetName);
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nNodes";
    nNodes = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/dt";
    dt = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/Trapping";
    Trapping = H5File_in.ReadScalar(dsetName);
    
    // Initialize the size.
    elemNodeConn.resize(nElements);  
    elemConDof.resize(nElements);
    elemIDs.resize(nElements);

    // Read global element IDs
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/ElementSetIDs";
    H5File_in.ReadField1D(dsetName, elemIDs);

    // Read node connectivity
    dsetName = "SimulationParameters/nTotElements";
    int totElements = H5File_in.ReadScalar(dsetName);  // Total number of elements
    vector<vector<int>> totElNodConnectivity(totElements);  // Node connectivity for all elements
    dsetName = "NodeConnectivity";
    H5File_in.ReadField2D(dsetName, totElements, nElNodes, totElNodConnectivity);

    // Node connectivity for the element set
    for (int iElem=0; iElem<nElements; iElem++){

        elemNodeConn.at(iElem) = totElNodConnectivity.at(elemIDs.at(iElem));
        elemConDof.at(iElem) = totElNodConnectivity.at(elemIDs.at(iElem));
    }

    // TODO: For debug!
    // for (auto& s : elemNodeConn[0])
    //     cout << s << "\n"; 
}