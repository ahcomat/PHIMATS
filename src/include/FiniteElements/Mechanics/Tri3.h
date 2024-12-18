/**
 * @file Tri3.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing tri3 elements.
 * @date 2024-05-23
 * 
 * @copyright Copyright (c) 2024
 *  
 */

#ifndef TRI3_H
#define TRI3_H

#include"BaseElemMech.h"
#include"Nodes.h"

using namespace std;

class Tri3: public BaseElemMech{

public:

Tri3(H5IO &H5File_in, Nodes &Nodes, int iSet);   

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

void CalcElStran(const double* globalBuffer) override;

void CalcRetrunMapping(vector<BaseMechanics*> mats) override;

private:

const vector<double> wts{0.5};  /// @brief Weights of the gauss points [nElGauss].

vector<RowVecd3> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<Matd2x3> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 
vector<Matd3x2> elemNodCoord;     /// @brief Node Coordinates [nElDim, nElNodes]. 
vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 

vector<vector<ColVecd3>> elStran;   /// @brief Int-pt strains [nElStres].
vector<vector<ColVecd3>> elStres;   /// @brief Int-pt stresses [nElStres].

vector<vector<Matd2x3>> BMat;       /// @brief Derivatives (scalar) matrix [nElDim, nElNodes].
vector<vector<Matd3x6>> BuMat;      /// @brief Strain matrix [nElStres, nElDispDofs].
vector<vector<double>> intPtVol;    /// @brief Int-pt volume.
vector<Matd6x6> elStiffMatx;        /// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].

};
#endif