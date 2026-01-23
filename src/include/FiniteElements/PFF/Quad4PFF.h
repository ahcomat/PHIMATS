/**
 * @file Quad4PFF.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Quad4PFF elements phase-field fracture. 
 * @date 2025-09-21
 * 
 * @copyright Copyright (C) 2025 Abdelrahman Hussein
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *s
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  
 */

#ifndef QUAD4PFF_H
#define QUAD4PFF_H

#include "BaseElemPFF.h"
#include "Nodes.h"

class Quad4PFF: public BaseElemPFF{

public:

Quad4PFF(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, Logger& logger);   

~Quad4PFF() override ;

/**
 * @brief Initializes the data for the shape functions.
 * 
 */
void InitShapeFunc();

/**
 * @brief Returns the shape functions of int-pts in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return RowVecd4
 */
RowVecd4  CalcShapeFunc(double xi, double eta);

/**
 * @brief Returns the shape function derivatives of int-pts in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd2x4 
 */
Matd2x4 CalcShapeFuncDeriv(double xi, double eta);

/**
 * @brief Initializes the data for all elements: `elemNodCoord`, `gaussPtCart`, `BMat`
 *  and `intPtVol`.
 * 
 * @param Nodes 
 */
void InitializeElements(Nodes& Nodes, H5IO &H5File_in);

/**
 * @brief Get the cartesian coordinates of gauss points `N_i x_ij`.
 * 
 * @param elCoord 
 * @param sFunc 
 * @return RowVecd2 
 */
RowVecd2 getGaussCart(RowVecd4& sFunc, Matd4x2& elCoord);

/**
 * @brief Calculates `intPtVol` and `BMat`.
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param wt 
 * @param intVol 
 * @param cartDeriv 
 */
void CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv);

void CalcPsiSpectral(double lam, double Gmod, const T_elStres& elStrain_e) override;

void CalcElemStiffMatx() override;

void CalcFH(double* FH) override;

void Calc_gPhi_d(const double* globalBuffer) override;

private:

/// @brief Weights of the gauss points [nElGauss].
const vector<double> wts{1.0, 1.0, 1.0, 1.0};  

/// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<RowVecd4> shapeFunc; 

/// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd2x4> shapeFuncDeriv; 

/// @brief Node Coordinates [nElNodes, nElDim]. 
vector<Matd4x2> elemNodCoord;   

/// @brief Cartesian coordinates of Gauss points for all elements [nElDim].
vector<vector<RowVecd2>> gaussPtCart;   

/// @brief Derivatives (scalar) matrix [nElDim, nElNodes]
vector<vector<Matd2x4>> BMat;              

/// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].
vector<Matd4x4> elStiffMatx;           

};
#endif