#include "FiniteElements/PFF/BaseElemPFF.h"

/*
 Unifying indices:
 i -> number of nodes.
 j -> spatial dimensions.
 k -> number of stresses.
 l -> total displacement dofs.
*/

void BaseElemPFF::ReadElementsData(H5IO &H5File_in, H5IO &H5File_mesh, int iSet){

    string dsetName;
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nElements";
    nElements = H5File_mesh.ReadScalar(dsetName);
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nNodes";
    nNodes = H5File_mesh.ReadScalar(dsetName);
    const_ell = H5File_in.ReadScalar("Materials/Material_" + to_string(iSet)+"/PFF/const_ell");
    
    // Initialize the size.
    elemNodeConn.resize(nElements);  
    elemPhiDof.resize(nElements);
    elemIDs.resize(nElements);

    // Read global element IDs
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/ElementSetIDs";
    H5File_mesh.ReadField1D(dsetName, elemIDs);

    // Read node connectivity
    dsetName = "SimulationParameters/nTotElements";
    int totElements = H5File_in.ReadScalar(dsetName);  // Total number of elements
    vector<vector<int>> totElNodConnectivity(totElements);  // Node connectivity for all elements
    dsetName = "NodeConnectivity";
    H5File_mesh.ReadField2D(dsetName, totElements, nElNodes, totElNodConnectivity);

    // Node connectivity for the element set
    for (int iElem=0; iElem<nElements; iElem++){

        elemNodeConn.at(iElem) = totElNodConnectivity.at(elemIDs.at(iElem));
        elemPhiDof.at(iElem) = totElNodConnectivity.at(elemIDs.at(iElem));
    }

    // TODO: For debug!
    // for (auto& s : elemNodeConn[0])
    //     cout << s << "\n"; 
}

const std::vector<std::vector<double>>& BaseElemPFF::getEl_gPhi_d() const {
    assert(el_gPhi_d_ptr != nullptr && "el_gPhi_d_ptr is null!");
    return *el_gPhi_d_ptr;
}

const std::vector<std::vector<double>>& BaseElemPFF::getElphi() const {
    assert(elPhi_ptr != nullptr && "elPhi_ptr is null!");
    return *elPhi_ptr;
}
