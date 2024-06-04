/**
 * @file Hex8.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing hex8 elements. 
 * @date 2024-05-23
 * 
 * @copyright Copyright (c) 2024
 *  
 */

#ifndef HEX8_H
#define HEX8_H

#include "BaseElemMech.h"
#include "Nodes.h"

using namespace std;

class Hex8: public BaseElemMech{

public:

Hex8(H5IO &H5File_in, Nodes &Nodes);   
~Hex8() override ;

/**
 * @brief Initializes the data for the shape functions.
 * 
 */
void InitShapeFunc();

/**
 * @brief Returns the integration pointe values of shape functions in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return RowVecd4
 */
RowVecd4  CalcShapeFunc(double xi, double eta);

/**
 * @brief Returns the integration pointe values of of shape function derivatives in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd2x4 
 */
Matd2x4 CalcShapeFuncDeriv(double xi, double eta);

private:

const vector<double> wts{1.0, 1.0, 1.0, 1.0};  /// @brief Weights of the gauss points.

vector<RowVecd4> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates.
vector<Matd2x4> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates. 
vector<Matd4x2> elemNodCoord;     /// @brief Node Coordinates. 
vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements. 

};
#endif