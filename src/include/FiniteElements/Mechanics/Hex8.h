/**
 * @file Hex8.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Hex8 elements. 
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

#ifndef HEX8_H
#define HEX8_H

#include "BaseElemMech.h"
#include "Nodes.h"

using namespace std;

class Hex8: public BaseElemMech{

public:

/**
 * @brief Constructor. Reads data for element set `iSet` in the elements vector and initializes the element set attributes.
 * 
 * @param H5File_in Input hdf5 file
 * @param Nodes Node object. Holds element node coordinates. 
 * @param iSet Element set number.
 * @param matModel material model for the element set. Options [`Elastic`, `ElasoPlastic`]. Default `Elastic`.
 */
Hex8(H5IO &H5File_in, Nodes &Nodes, int iSet, string matModel, Logger& logger);   
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

void CalcNodVals(T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<double>& nodStres_h, vector<double>& nodRho, vector<int>& nodCount) override;

void CalcElDStran(const double* globalBuffer) override;

void CalcElStran(const double* globalBuffer) override;

void CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep) override;

void CalcFint(double* Fint) override;

void getNew() override;

private:

/// @brief Weights of the gauss points [nElGauss].
const vector<double> wts{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};  

/// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<RowVecd8> shapeFunc; 

/// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd3x8> shapeFuncDeriv;  

/// @brief Node Coordinates [nElDim, nElNodes]. 
vector<Matd8x3> elemNodCoord;  

/// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 
vector<vector<RowVecd3>> gaussPtCart;  

/// @brief Int-pt total strains [nElStres].
vector<vector<ColVecd6>> elStran;   

/// @brief Int-pt strain increment [nElStres].
vector<vector<ColVecd6>> elDStran; 

/// @brief Int-pt stresses [nElStres].
vector<vector<ColVecd6>> elStres;   

/// @brief Int-pt elastic strain [nElStres].
vector<vector<ColVecd6>> elStran_e;   

/// @brief Int-pt plastic strain [nElStres].
vector<vector<ColVecd6>> elStran_p;  

/// @brief Int-pt elastic strain (last convergged increment) [nElStres].
vector<vector<ColVecd6>> elStran_e_old;   

/// @brief Int-pt plastic strain (last convergged increment) [nElStres].
vector<vector<ColVecd6>> elStran_p_old;  

/// @brief Int-pt equivalent plastic strain (last convergged increment) [nElStres].
vector<vector<double>> elStran_eq_old;

/// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<Matd3x8>> BMat;   

/// @brief Strain matrix [nElStres, nElDispDofs].
vector<vector<Matd6x24>> BuMat;  

/// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].
vector<Matd24x24> elStiffMatx;  

};
#endif