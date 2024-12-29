/**
 * @file Quad4.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing quad4 elements. 
 * @date 2024-05-22
 * 
 * @copyright Copyright (C) 2024 Abdelrahman Hussein
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

#ifndef QUAD4_H
#define QUAD4_H

#include"BaseElemMech.h"
#include"Nodes.h"

class Quad4: public BaseElemMech{

public:

Quad4(H5IO &H5File_in, Nodes &Nodes, int iSet, string matModel);   

~Quad4() override ;

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
 * @brief Initializes the data for all elements: `elemNodCoord`, `gaussPtCart`, `BMat`, `BuMat`
 *  and `intPtVol`.
 * 
 * @param Nodes 
 */
void InitializeElements(Nodes& Nodes);

/**
 * @brief Get the cartesian coordinates of gauss points `N_i x_ij`.
 * 
 * @param elCoord 
 * @param sFunc 
 * @return RowVecd2 
 */
RowVecd2 getGaussCart(RowVecd4& sFunc, Matd4x2& elCoord);

/**
 * @brief Calculates `intPtVol`, `BMat` and `BuMat`.
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param intVol 
 * @param cartDeriv 
 * @param strainMat 
 */
void CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv, Matd3x8& strainMat);

/**
 * @brief Calculates the element stiffness matrix for all elements.
 * 
 * @param DMatx 
 */
void CalcElemStiffMatx(T_DMatx DMatx) override ;

void CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount) override;

void CalcNodVals(T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<int>& nodCount) override;

void CalcElStran(const double* globalBuffer) override;

void CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep) override;

void CalcFint(double* Fint) override;

private:

/// @brief Weights of the gauss points [nElGauss].
const vector<double> wts{1.0, 1.0, 1.0, 1.0};  

/// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<RowVecd4> shapeFunc; 

/// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd2x4> shapeFuncDeriv;  

/// @brief Node Coordinates [nElDim, nElNodes]. 
vector<Matd4x2> elemNodCoord;  

/// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 
vector<vector<RowVecd2>> gaussPtCart;  

/// @brief Int-pt strains [nElStres].
vector<vector<ColVecd3>> elStran;   

/// @brief Int-pt stresses [nElStres].
vector<vector<ColVecd3>> elStres; 

/// @brief Int-pt elastic strain [nElStres].
vector<vector<ColVecd3>> elStran_e;   

/// @brief Int-pt plastic strain [nElStres].
vector<vector<ColVecd3>> elStran_p;  

/// @brief Int-pt equivalent plastic strain [nElStres].
vector<vector<double>> elStran_eq;
    
/// @brief Int-pt equivalent stress (von Mises) [nElStres].
vector<vector<double>> elStres_eq;

/// @brief Int-pt hydrostatic stress [nElStres].
vector<vector<double>> elStres_h;

/// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<Matd2x4>> BMat;

/// @brief Strain matrix [nElStres, nElDispDofs].
vector<vector<Matd3x8>> BuMat; 

/// @brief Int-pt volume.
vector<vector<double>> intPtVol;  

/// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].
vector<Matd8x8> elStiffMatx;        

};
#endif