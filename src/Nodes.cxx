#include "Nodes.h"
#include <iostream>

#ifndef DEBUG
#define at(x) operator[](x)
#endif

Nodes::~Nodes(){}

void Nodes::ReadNodes(H5IO &H5File_in, H5IO &H5File_mesh){

    string dsetName;
    dsetName = "SimulationParameters/nTotNodes";
    nTotNodes = H5File_in.ReadScalar(dsetName);

    dsetName = "SimulationParameters/nDim";
    nDim = H5File_in.ReadScalar(dsetName);

    nodeCoordinates.resize(nTotNodes);

    dsetName = "NodeCoordinates";
    H5File_mesh.ReadField2D(dsetName, nTotNodes, nDim, nodeCoordinates);
}

vector<double> Nodes::getNodCoord(int nod){

    return nodeCoordinates.at(nod); 
}

int Nodes::get_nNodes(){

    return nTotNodes;
}