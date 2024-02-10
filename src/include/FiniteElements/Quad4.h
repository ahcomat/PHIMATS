#ifndef QUAD4_H
#define QUAD4_H

#include "Eigen/Dense"
#include "petsc.h"
#include "Nodes.h"
#include "H5IO.h"

// /**
//  * @brief Storage for element data. 
//  * 
// //  */
// struct elemData {
//   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//   Eigen::Vector<double, 4> shapFunc;
//   Eigen::Matrix<double, 2, 4> shapeFuncDeriv;
//   vector<int> nConnect;
//   vector<vector<double>> nCoord;
// };

// Element Specific data.
const int nDim = 2;           /// Spatial dimensions of the element.
const int nElNodes = 4;       /// Number of nodes per element.
const int dispDofs = 2;       /// Number of displacement dofs. 
const int nStres = 3;         /// Stress/strain components.
const int nElDispDofs = 8;    /// Number of element displacement dofs.
const int nGauss = 4;         /// Number of gauss points.

/**
 * Eigen Vector/Matrix typedefs
 */
typedef Eigen::RowVector<double, nDim> RowVecd2;                // nDim 
typedef Eigen::Vector<double, nStres> ColVecd3;                 // nStres 
typedef Eigen::RowVector<double, nElNodes> RowVecd4;            // nElNodes 
typedef Eigen::Vector<double, nElDispDofs> ColVecd8;            // nElDispDofs 


typedef Eigen::Matrix<double, nDim, nDim> Matd2x2;                // nDim x nDim 
typedef Eigen::Matrix<double, nStres, nStres> Matd3x3;            // nDim x nDim 
typedef Eigen::Matrix<double, nElNodes, nDim> Matd4x2;            // nElNodes x nDim 
typedef Eigen::Matrix<double, nDim, nElNodes> Matd2x4;            // nDim x nElNodes
typedef Eigen::Matrix<double, nStres, nElDispDofs> Matd3x8;       // nStres x nElDispDofs
typedef Eigen::Matrix<double, nElDispDofs, nElDispDofs,           // nStres x nElDispDofs
Eigen::RowMajor> Matd8x8;  


using namespace std;


class Quad4
{

public:

Quad4(H5IO &H5File_in, Nodes &Nodes);   

~Quad4();

/**
 * @brief Initializes the data for the shape functions.
 * 
 */
void InitShapeFunc();

/**
 * @brief Returns the values of shape functions in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return RowVecd4
 */
RowVecd4 getShapeFunc(double xi, double eta);

/**
 * @brief Returns the derivatives of the shape functions in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd2x4 
 */
Matd2x4 getShapeFuncDeriv(double xi, double eta);

/**
 * @brief Reads the data `nElements`, `nElementSets` and `elemNodeConn` from hdf5 file.
 * 
 * @param H5File_in 
 */
void ReadElementsData(H5IO &H5File_in);

/**
 * @brief Get the #DOfs
 * 
 * @return int 
 */
int getTotDof();

/**
 * @brief Returns the node connectivity of element iElem.
 * 
 * @param iElem 
 * @return vector<int> 
 */
vector<int> getNodeConnect(int iElem);

/**
 * @brief Return the displacement dofs associated with element `iElem`.
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
void InitializeElements(Nodes &Nodes);

/**
 * @brief Get the node coordinates of a given element.
 * 
 * @param iElem 
 * @return Matd4x2 
 */
Matd4x2 getElemNodCoord(int iElem);

/**
 * @brief Get the cartesian coordinates of gauss points.
 * 
 * @param elCoord 
 * @param sFunc 
 * @return RowVecd2 
 */
RowVecd2 getGaussCart(Matd4x2 elCoord, RowVecd4 sFunc);

/**
 * @brief 
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param intVol 
 * @param cartDeriv 
 * @param strainMat 
 */
void CalcCartDeriv(Matd4x2 elNodCoord, Matd2x4 sFuncDeriv, double &intVol, Matd2x4 &cartDeriv, Matd3x8 &strainMat);

/**
 * @brief Evaluates the element stiffness matrix for all elements.
 * 
 * @param DMatx 
 */
void CalcElemStiffMatx(Matd3x3 DMatx);

/**
 * @brief Return the vector of element stiffness matrix.
 * 
 * @return vector<Matd8x8> 
 */
vector<Matd8x8>getElemStiffMatx();

/**
 * @brief Initializes Petsc data structures.
 * 
 */
void InitPETSC();

/**
 * @brief Assemble the global stiffness matrix
 * 
 */
void Assemble();

/**
 * @brief Set Dirichlet boundary conditions.
 * 
 */
void setDirichBC();

Vec& getB();

Vec& getX();

Mat& getA();

/**
 * @brief Calculates the strains and stresses.
 * 
 */
void CalcStres(Matd3x3 DMatx, bool nodStrFlag=false);


/**
 * @brief Write simulation data to hdf5 file.
 * 
 * @param H5File_out 
 */
void WriteOut(H5IO &H5File_out);

private:

// Input data size.
int nElementSets;   /// Number of element sets
int nElements;      /// Total number of elements.

int nNodes;   /// Total number of nodes.
int nTotDof;  /// Total number of DOFs.

int nPresDofs; /// Number of prescribed displacement dofs.

vector<vector<int>> elemNodeConn; /// Node connectivity.
vector<vector<int>> elemDispDof;  /// Element displacement dofs
vector<Matd4x2> elemNodCoord;  /// Node Coordinates. 

vector<vector<double>> gaussPts;    /// Gauss points in natural coordinates. 
vector<double> wts{1.0, 1.0, 1.0, 1.0};  /// Weights of the gauss points.
vector<vector<RowVecd2>> gaussPtCart;  /// Gauss points in cartesian coordinates of all elements. 

vector<RowVecd4> shapeFunc; /// Shape functions at integration points.
vector<Matd2x4> shapeFuncDeriv;  // Shape functions derivatives at integration points. 

// vector<vector<Eigen::Vector<double, 4>>> ShapeFuncCart;    /// Cartesian shape functions at integration points.

vector<vector<Matd2x4>> BMat;       /// Derivatives matrix.
vector<vector<Matd3x8>> BuMat;      /// Strain matrix.
vector<vector<double>> intPtVol;    /// Integration point volume.

vector<Matd8x8> elStiffMatx;    /// Element stiffness matrix.

vector<vector<ColVecd3>> elStres, elStran;  /// Integration point stresses and strains.

vector<ColVecd3> nodStres, nodStran;  /// Nodal point stresses and strains.
vector<int> nodCount;                 /// Counter for nodes.


// Petsc ----------------------------

const PetscScalar* globalBuffer;

// Dirichlet boundary conditions
PetscInt  *presDofs;    // Array to hold the prescribed dofs.
PetscScalar  *presVals; // Array to hold the prescribed values.
PetscScalar  *Fint;     // Array to hold the internal force vector.

Vec x, b; // The RHS and solution.

Mat A; // The global coefficient (stiffness) matrix.

};
#endif