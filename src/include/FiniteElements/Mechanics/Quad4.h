/**
 * @file Quad4.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing quad4 elements. 
 * @date 2024-05-22
 * 
 * @copyright Copyright (c) 2024
 *
 */

#ifndef QUAD4_H
#define QUAD4_H

#include"BaseElemMech.h"
#include"Nodes.h"

class Quad4: public BaseElemMech{

public:

Quad4(H5IO &H5File_in, Nodes &Nodes, int iSet);   

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

/**
 * @brief Calculates the Fint, strains and stresses. Also Calculates the stress nodal values 
 *        if `nodStresFlag=true`.
 * 
 */
void CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount) override;

void CalcElStran(const double* globalBuffer) override;

private:

const vector<double> wts{1.0, 1.0, 1.0, 1.0};  /// @brief Weights of the gauss points [nElGauss].

vector<RowVecd4> shapeFunc;      /// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<Matd2x4> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd4x2> elemNodCoord;     /// @brief Node Coordinates [nElDim, nElNodes]. 
vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 

vector<vector<ColVecd3>> elStran;   /// @brief Int-pt strains [nElStres].
vector<vector<ColVecd3>> elStres;   /// @brief Int-pt stresses [nElStres].

vector<vector<Matd2x4>> BMat;       /// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<Matd3x8>> BuMat;      /// @brief Strain matrix [nElStres, nElDispDofs].
vector<vector<double>> intPtVol;    /// @brief Int-pt volume.
vector<Matd8x8> elStiffMatx;        /// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].

};
#endif