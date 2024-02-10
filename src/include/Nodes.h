#ifndef NODES_H
#define NODES_H

#include <vector>
#include "H5IO.h"

using namespace std;

class Nodes{

public:

Nodes(){};

~Nodes();

/**
 * @brief Read node coordinates from hdf5 file. 
 * 
 * @param H5File 
 */
void ReadNodes(H5IO &H5File);

/**
 * @brief get the coordinates of node "nod".
 * 
 * @param nod : Node number.
 * @return vector<double>.
 */
vector<double> getNodCoord(int nod);

/**
 * @brief returns the total number of nodes.
 * 
 */
int getNNodes();

private:

int nNodes;      /// Total number of nodes.
int nDim;        /// Simulation dimensions.

vector<vector<double>> nodeCoordinates; 

};
#endif