/**
 * @file Quad4T.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Quad4T elements for heat and mass transfer. 
 * @date 2024-06-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef QUAD4T_H
#define QUAD4T_H

#include"BaseElemTrans.h"
#include"Nodes.h"

class Quad4T: public BaseElemTrans{

public:

Quad4T(H5IO &H5File_in, Nodes &Nodes, int iSet);   

~Quad4T() override ;

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
void CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv);

/**
 * @brief Calculates the element stiffness matrix for all elements.
 * 
 * @param DMatx 
 */
void CalcElemStiffMatx(T_DMatx DMatx, double s) override ;

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
void CalcFlux(T_DMatx KMatx, const double* globalBuffer, T_nodFlux& nodFlux, vector<int>& nodCount) override;

private:

const vector<double> wts{1.0, 1.0, 1.0, 1.0};  /// @brief Weights of the gauss points [nElGauss].

vector<RowVecd4> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<Matd2x4> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd4x2> elemNodCoord;     /// @brief Node Coordinates [nElDim, nElNodes]. 
vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 

vector<vector<ColVecd2>> elFlux;   /// @brief Int-pt strains [nElStres].

vector<vector<Matd2x4>> BMat;       /// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<double>> intPtVol;    /// @brief Int-pt volume.
vector<Matd4x4> elStiffMatx;        /// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].
vector<Matd4x4> elCapMatx;          /// @brief Element capacitance matrix [nElDispDofs, nElDispDofs].

};
#endif