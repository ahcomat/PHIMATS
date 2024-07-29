/**
 * @file Tri3T.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Tri3T elements for heat and mass transfer with 3 integration points. 
 * @date 2024-07-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef TRI3T_H
#define TRI3T_H

#include"BaseElemTransport.h"
#include"Nodes.h"

class Tri3T: public BaseElemTransport{

public:

Tri3T(H5IO &H5File_in, Nodes &Nodes, int iSet);   

~Tri3T() override ;

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
RowVecd3  CalcShapeFunc(double xi, double eta);

/**
 * @brief Returns the shape function derivatives of int-pts in natural coordinates.
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
void CalcCartDeriv(Matd3x2& elNodCoord, Matd2x3& sFuncDeriv, const double& wt, double& intVol, Matd2x3& cartDeriv);

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
void CalcFlux(T_DMatx KMatx, const double* globalBuffer, T_nodStres& nodFlux, vector<double>& nodCount) override;

/**
 * @brief Calculates the volume averaged concentration
 * 
 * @param globalBuffer Solution vector
 * @return double Volume averaged concentration
 */
double CalcAvCon(const double* globalBuffer) override;

private:

const vector<double> wts{1.0/6.0, 1.0/6.0, 1.0/6};  /// @brief Weights of the gauss points [nElGauss].

vector<RowVecd3> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<Matd2x3> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd3x2> elemNodCoord;     /// @brief Node Coordinates [nElNodes, nElDim]. 
vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 

vector<vector<ColVecd2>> elFlux;   /// @brief Int-pt flux [nElStres].

vector<vector<Matd2x3>> BMat;       /// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<double>> intPtVol;    /// @brief Int-pt volume.
vector<Matd3x3> elStiffMatx;        /// @brief Element stiffness matrix [nElConDofs, nElConDofs].
vector<Matd3x3> elCapMatx;          /// @brief Element capacitance matrix [nElConDofs, nElConDofs].

};
#endif