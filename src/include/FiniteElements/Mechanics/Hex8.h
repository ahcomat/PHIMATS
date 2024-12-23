/**
 * @file Hex8.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing hex8 elements. 
 * @date 2024-05-23
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

#ifndef HEX8_H
#define HEX8_H

#include "BaseElemMech.h"
#include "Nodes.h"

using namespace std;

class Hex8: public BaseElemMech{

public:

Hex8(H5IO &H5File_in, Nodes &Nodes, int iSet, string matModel = "Elastic");   
~Hex8() override;

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
 * @return RowVecd8
 */
RowVecd8  CalcShapeFunc(double xi, double eta, double zeta);

/**
 * @brief Returns the integration pointe values of of shape function derivatives in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd3x8
 */
Matd3x8 CalcShapeFuncDeriv(double xi, double eta, double zeta);

/**
 * @brief Initializes the data for all elements: `elemNodCoord`, `gaussPtCart`, `BMat`, `BuMat`
 *  and `intPtVol`.
 * 
 * @param Nodes 
 */
void InitializeElements(Nodes &Nodes);

/**
 * @brief Get the cartesian coordinates of gauss points `N_i x_ij`.
 * 
 * @param elCoord 
 * @param sFunc 
 * @return RowVecd3 
 */
RowVecd3 getGaussCart(RowVecd8& sFunc, Matd8x3& elCoord);

/**
 * @brief Calculates `intPtVol`, `BMat` and `BuMat`.
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param intVol 
 * @param cartDeriv 
 * @param strainMat 
 */
void CalcCartDeriv(Matd8x3& elNodCoord, Matd3x8& sFuncDeriv, const double& wt, double& intVol, Matd3x8& cartDeriv, Matd6x24& strainMat);

void CalcElemStiffMatx(T_DMatx DMatx) override;

void CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount) override;

void CalcElStran(const double* globalBuffer) override;

void CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep) override;

void CalcFint(double* Fint) override;

private:

const vector<double> wts{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};  /// @brief Weights of the gauss points [nElGauss].

vector<RowVecd8> shapeFunc;      /// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<Matd3x8> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd8x3> elemNodCoord;          /// @brief Node Coordinates [nElDim, nElNodes]. 
vector<vector<RowVecd3>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 

vector<vector<ColVecd6>> elStran;   /// @brief Int-pt total strains [nElStres].
vector<vector<ColVecd6>> elStres;   /// @brief Int-pt stresses [nElStres].

vector<vector<ColVecd6>> elStran_e;   /// @brief Int-pt elastic strain [nElStres].
vector<vector<ColVecd6>> elStran_p;   /// @brief Int-pt plastic strain [nElStres].
vector<vector<double>> elStran_eq;    /// @brief Int-pt equivalent plastic strain [nElStres].
vector<vector<double>> elStres_eq;    /// @brief Int-pt equivalent stress (vin Mises) [nElStres].

vector<vector<Matd3x8>> BMat;       /// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<Matd6x24>> BuMat;     /// @brief Strain matrix [nElStres, nElDispDofs].
vector<vector<double>> intPtVol;    /// @brief Int-pt volume.
vector<Matd24x24> elStiffMatx;      /// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].
};
#endif