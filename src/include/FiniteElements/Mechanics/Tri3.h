/**
 * @file Tri3.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Tri3 elements.
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

#ifndef TRI3_H
#define TRI3_H

#include"BaseElemMech.h"
#include"Nodes.h"

using namespace std;

class Tri3: public BaseElemMech{

public:

Tri3(H5IO &H5File_in, Nodes &Nodes, int iSet, string matModel, Logger& logger);   

~Tri3() override ;

/**
 * @brief Initializes the data for the shape functions.
 * 
 */
void InitShapeFunc();

/**
 * @brief Returns the int-pt values of shape functions in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return RowVecd3
 */
RowVecd3  CalcShapeFunc(double xi, double eta);

/**
 * @brief Returns the int-pt values of of shape function derivatives in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd2x3 
 */
Matd2x3 CalcShapeFuncDeriv(double xi, double eta);

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
RowVecd2 getGaussCart(RowVecd3& sFunc, Matd3x2& elCoord);

/**
 * @brief Calculates `intPtVol`, `BMat` and `BuMat`.
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param intVol 
 * @param cartDeriv 
 * @param strainMat 
 */
void CalcCartDeriv(Matd3x2& elNodCoord, Matd2x3& sFuncDeriv, const double& wt, double& intVol, Matd2x3& cartDeriv, Matd3x6& strainMat);

void CalcElemStiffMatx(T_DMatx DMatx) override ;

void CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount) override;

void CalcNodVals(T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<double>& nodStres_h, vector<double>& nodRho, vector<int>& nodCount) override;


void CalcElStran(const double* globalBuffer) override;

void CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep) override;

private:

/// @brief Weights of the gauss points [nElGauss].
const vector<double> wts{0.5};  

/// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<RowVecd3> shapeFunc; 

/// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd2x3> shapeFuncDeriv; 

/// @brief Node Coordinates [nElDim, nElNodes]. 
vector<Matd3x2> elemNodCoord;     

/// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 
vector<vector<RowVecd2>> gaussPtCart;  

/// @brief Int-pt strains [nElStres].
vector<vector<ColVecd3>> elStran;   

/// @brief Int-pt stresses [nElStres].
vector<vector<ColVecd3>> elStres;   

/// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<Matd2x3>> BMat;       

/// @brief Strain matrix [nElStres, nElDispDofs].
vector<vector<Matd3x6>> BuMat;      

/// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].
vector<Matd6x6> elStiffMatx;        

};
#endif