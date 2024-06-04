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

using namespace std;

class Quad4: public BaseElemMech{

public:

Quad4(H5IO &H5File_in, Nodes &Nodes);   

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
 * @brief Reads the data `nElements`, `nElementSets` and `elemNodeConn` from hdf5 file.
 * 
 * @param H5File_in  
 */
void ReadElementsData(H5IO &H5File_in);

/**
 * @brief Returns the displacement dofs associated with element `iElem`.
 * 
 * @param iElem 
 * @return vector<int> 
 */
vector<int> CalcElemDispDof(int iElem);

// /**
//  * @brief Returns the node connectivity of element iElem.
//  * 
//  * @param iElem 
//  * @return vector<int> 
//  * @todo Possibly move to base `Elements`.
//  */
// vector<int> getNodeConnect(int iElem);

/**
 * @brief Initializes the data for all elements: `elemNodCoord`, `gaussPtCart`, `BMat`, `BuMat`
 *  and `intPtVol`.
 * 
 * @param Nodes 
 */
void InitializeElements(Nodes& Nodes);

// /**
//  * @brief Get the node coordinates of a given element.
//  * 
//  * @param iElem 
//  * @return Matd4x2 
//  */
// Matd4x2 getElemNodCoord(int iElem);

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
void CalcStres(T_DMatx DMatx, const double* globalBuffer) override;

/**
 * @brief Write nodal stresses and strains.
 * 
 * @param H5File_out 
 */
void WriteOut(H5IO &H5File_out) override;

private:

const vector<double> wts{1.0, 1.0, 1.0, 1.0};  /// @brief Weights of the gauss points.

vector<RowVecd4> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates.
vector<Matd2x4> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates. 
vector<Matd4x2> elemNodCoord;     /// @brief Node Coordinates. 
vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements. 

// vector<vector<Eigen::Vector<double, 4>>> ShapeFuncCart;    /// Cartesian shape functions at integration points.

vector<vector<ColVecd3>> elStran;   /// @brief Int-pt strains.
vector<vector<ColVecd3>> elStres;   /// @brief Int-pt stresses.
vector<ColVecd3> nodStran;          /// @brief Nod-pt strains.
vector<ColVecd3> nodStres;          /// @brief Nod-pt stresses.
vector<int> nodCount;     /// @brief Counter for integration points surrounding nodes.

vector<vector<Matd2x4>> BMat;       /// @brief Derivatives (scalar) matrix.
vector<vector<Matd3x8>> BuMat;      /// @brief Strain matrix.
vector<vector<double>> intPtVol;    /// @brief Int-pt volume.
vector<Matd8x8> elStiffMatx;        /// @brief Element stiffness matrix.

// Petsc ----------------------------

PetscScalar* Fint;     // Array to hold the internal force vector.

};
#endif