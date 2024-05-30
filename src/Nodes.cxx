#include "Nodes.h"
#include <iostream>

#ifndef DEBUG
#define at(x) operator[](x)
#endif

Nodes::~Nodes(){

    // Exit message
    cout << "Nodes exited correctly" << "\n";
}

void Nodes::ReadNodes(H5IO &H5File){

    string dsetName;
    dsetName = "SimulationParameters/nNodes";
    nNodes = H5File.ReadScalar(dsetName);

    dsetName = "SimulationParameters/nDim";
    nDim = H5File.ReadScalar(dsetName);

    vector<double> dummy1(nDim);

    for (int iNod=0; iNod<nNodes; iNod++){
        dsetName = "NodeCoordinates/Node_"+to_string(iNod);
        H5File.ReadFieldDoub1D(dsetName, dummy1);
        nodeCoordinates.push_back(dummy1);
    }
}

vector<double> Nodes::getNodCoord(int nod){

    return nodeCoordinates.at(nod); 
}

int Nodes::getNNodes(){

    return nNodes;
}