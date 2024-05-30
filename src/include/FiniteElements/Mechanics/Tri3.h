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

Tri3(H5IO &H5File_in, Nodes &Nodes);   

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
RowVecd3 getShapeFunc(double xi, double eta);

/**
 * @brief Returns the int-pt values of of shape function derivatives in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd2x3 
 */
Matd2x3 getShapeFuncDeriv(double xi, double eta);

/**
 * @brief Reads the data `nElements`, `nElementSets` and `elemNodeConn` from hdf5 file.
 * 
 * @param H5File_in 
 * @todo Consider moving to base `Elements`. `getElemDispDof` will probably have to be pure `virtual` function.
 */
void ReadElementsData(H5IO &H5File_in);

/**
 * @brief Returns the displacement dofs associated with element `iElem`.
 * 
 * @param iElem 
 * @return vector<int> 
 */
vector<int> getElemDispDof(int iElem);

/**
 * @brief Initializes the data for all elements: `elemNodCoord`, `gaussPtCart`, `BMat`, `BuMat`
 *  and `intPtVol`.
 * 
 * @param Nodes 
 */
void InitializeElements(Nodes& Nodes);

void CalcElemStiffMatx(T_DMatx DMatx) override;

void CalcStres(T_DMatx DMatx, const double* globalBuffer, bool nodStresFlag=false) override;

void WriteOut(H5IO &H5File_out) override;

private:

const vector<double> wts{0.5};  /// @brief Weights of the gauss points.

vector<RowVecd3> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates.
vector<Matd2x3> shapeFuncDeriv;  /// @brief Values of the shape function derivatives at integration points in natural coordinates. 
vector<Matd3x2> elemNodCoord;     /// @brief Node Coordinates. 
vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements. 

vector<vector<ColVecd3>> elStran;   /// @brief Int-pt strains.
vector<vector<ColVecd3>> elStres;   /// @brief Int-pt stresses.
vector<ColVecd3> nodStran;          /// @brief Nod-pt strains.
vector<ColVecd3> nodStres;          /// @brief Nod-pt stresses.
vector<int> nodCount;     /// @brief Counter for integration points surrounding nodes.

vector<vector<Matd2x3>> BMat;       /// @brief Derivatives (scalar) matrix.
vector<vector<Matd3x6>> BuMat;      /// @brief Strain matrix.
vector<vector<double>> intPtVol;    /// @brief Int-pt volume.
vector<Matd6x6> elStiffMatx;        /// @brief Element stiffness matrix.

// Petsc ----------------------------

PetscScalar* Fint;     // Array to hold the internal force vector.

};
#endif