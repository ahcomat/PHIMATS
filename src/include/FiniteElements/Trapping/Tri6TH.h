/**
 * @file Tri6TH.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Tri6TH elements for full-fileld hydrogen trapping. 
 * @date 2025-01-16
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

#ifndef TRI6TH_H
#define TRI6TH_H

#include "BaseElemTrap.h"
#include "Nodes.h"
#include "Version.h"

class Tri6TH: public BaseElemTrap{

public:

Tri6TH(H5IO &H5File_in, Nodes &Nodes, int iSet, Logger& logger);   

~Tri6TH() override ;

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
 * @return RowVecd3
 */
RowVecd6 CalcShapeFunc(double xi, double eta);

/**
 * @brief Returns the shape function derivatives of int-pts in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd2x3 
 */
Matd2x6 CalcShapeFuncDeriv(double xi, double eta);

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
RowVecd2 getGaussCart(RowVecd6& sFunc, Matd6x2& elCoord);

/**
 * @brief Calculates `intPtVol`, `BMat` and `BuMat`.
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param intVol 
 * @param cartDeriv 
 * @param strainMat 
 */
void CalcCartDeriv(Matd6x2& elNodCoord, Matd2x6& sFuncDeriv, const double& wt, double& intVol, Matd2x6& cartDeriv);

/**
 * @brief Returns the int-pt cartesian coordinates
 * 
 * @param glIntPtCoords 
 */
void getInPtCoords(T_nodStres& glIntPtCoords) override;

// /**
//  * @brief Evaluates the gradients of phi at the int-points and maps them to the nodes. 
//  * 
//  * @param nodGrad
//  * @param nodCount 
//  * 
//  * @todo Move to BaseElemTrap.
//  */
// void CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi) override;

/**
 * @brief Calculates the element stiffness and capacitance matrix for all elements.
 * 
 * @param mat Material
 * @param T Current temperature
 */
void CalcElemStiffMatx(BaseTrapping* mat, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr = nullptr) override;

/**
 * @brief Calculates the volume averaged concentration
 * 
 * @param globalBuffer Solution vector
 * @return Volume averaged concentration
 */
double CalcAvCon(const double* globalBuffer) override;

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
void CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr = nullptr) override;

void CalcFsrc(const double conB, BaseTrapping* mat, double* FsrcBuffer, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr) override;


private:

/// @brief Weights of the gauss points [nElGauss].
const vector<double> wts{1.0/6.0, 1.0/6.0, 1.0/6};  

/// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<RowVecd6> shapeFunc; 

/// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd2x6> shapeFuncDeriv; 

/// @brief Node Coordinates [nElNodes, nElDim]. 
vector<Matd6x2> elemNodCoord;   

/// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 
vector<vector<RowVecd2>> gaussPtCart;  

/// @brief Int-pt flux [nElDim]
vector<vector<ColVecd2>> elFlux;         
   
/// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<Matd2x6>> BMat;   

/// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].
vector<Matd6x6> elStiffMatx;  

/// @brief Element capacitance matrix [nElDispDofs, nElDispDofs].
vector<Matd6x6> elCapMatx;         

};
#endif