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
    dsetName = "SimulationParameters/nTotNodes";
    nTotNodes = H5File.ReadScalar(dsetName);

    dsetName = "SimulationParameters/nDim";
    nDim = H5File.ReadScalar(dsetName);

    nodeCoordinates.resize(nTotNodes);

    dsetName = "NodeCoordinates";
    H5File.ReadFieldFloat2D(dsetName, nTotNodes, nDim, nodeCoordinates);
}

vector<double> Nodes::getNodCoord(int nod){

    return nodeCoordinates.at(nod); 
}

int Nodes::get_nNodes(){

    return nTotNodes;
}