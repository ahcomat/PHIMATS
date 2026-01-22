/**
 * @file Nodes.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A class for reading and storing nodal coordinates. Works for 2D and 3D.
 * @date 2024-05-23
 * 
 * @copyright Copyright (C) 2025 Abdelrahman Hussein
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *   
 */

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
void ReadNodes(H5IO &H5File_in, H5IO &H5File_mesh);

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
int get_nNodes();

private:

int nTotNodes;      /// Total number of nodes.
int nDim;        /// Simulation dimensions.

vector<vector<double>> nodeCoordinates; 

};
#endif