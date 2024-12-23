/**
 * @file Quad4TH.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Quad4TH elements for hydrogen diffusion and trapping. 
 * @date 2024-06-26
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

#ifndef QUAD4TH_H
#define QUAD4TH_H

#include "BaseElemTrap.h"
#include "Nodes.h"

class Quad4TH: public BaseElemTrap{

public:

Quad4TH(H5IO &H5File_in, Nodes &Nodes, int iSet);   

~Quad4TH() override ;

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
 * @brief Calculates `intPtVol`, `BMat` and `BuMat`.
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param intVol 
 * @param cartDeriv 
 * @param strainMat 
 */
void CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv);

/**
 * @brief Evaluates the gradients of phi at the int-points and maps them to the nodes. 
 * 
 * @param nodGrad
 * @param nodCount 
 * 
 * @todo Move to BaseElemTrap.
 */
void CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi) override;

/**
 * @brief Get the int-pt coordinates
 * 
 * @param glIntPtCoords 
 */
virtual void getInPtCoords(T_nodStres& glIntPtCoords) override;

/**
 * @brief Calculates the element stiffness matrix for all elements.
 * 
 * @param DMatx 
 */
void CalcElemStiffMatx(BaseTrapping* mat, const double T) override ;

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
void CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T) override;

private:

const vector<double> wts{1.0, 1.0, 1.0, 1.0};  /// @brief Weights of the gauss points [nElGauss].


vector<RowVecd4> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<Matd2x4> shapeFuncDeriv; /// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 

vector<Matd4x2> elemNodCoord;   /// @brief Node Coordinates [nElDim, nElNodes]. 

vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 

vector<vector<ColVecd2>> elFlux; /// @brief Int-pt flux [nElStres]
vector<vector<double>> elPhi;    /// @brief Int-pt phi [nElStres]. 

vector<double> nodPhi;  /// @brief nodal values of phi [nTotNodes]

vector<vector<Matd2x4>> BMat;   /// @brief Derivatives (scalar) matrix [nElDim, nElNodes].

vector<vector<double>> intPtVol;    /// @brief Int-pt volume.       

vector<Matd4x4> elStiffMatx;    /// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].    
vector<Matd4x4> elCapMatx;      /// @brief Element capacitance matrix [nElDispDofs, nElDispDofs].   
vector<Matd4x4> elMKTMatx;       /// @brief Element trapping matrix [nElDispDofs, nElDispDofs].   

};
#endif