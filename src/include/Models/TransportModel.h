/**
 * @file TransportModel.h
 * @brief A class based on heat and mass transport models to interface with PETSc based global
 * capacitance matrix `C`, conductivity matrix `K`, solution vector `x` and RHS `F`.
 * 
 * @details Main functions:
 * - Calls `BaseElemTrans::CalcElemStiffMatx` to build the local diffusivity and capacitance matrix.
 * - Assembles the global capacitance and conductivity matrix.
 * - Initializes and applies boundary conditions.
 * - Manages the output by writing to H5file_out.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2024-06-15
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include "H5IO.h"
#include "FiniteElements/Transport/BaseElemTrans.h"
#include "Materials/Transport/BaseTransport.h"

class TransportModel{

public:

TransportModel(vector<BaseElemTrans*> elements, H5IO& H5File_in);
~TransportModel();

/**
 * @brief Set nodal quantities `nodStres`, `nodStran` and `nodCount` to zeros.
 * 
 */
void setZero_nodFlux();

/**
 * @brief Initializes and preallocates the RHS `F`, solution `x` and 
 *        conductivity matrix `K` and capacitance matrix `C`.
 * 
 * @param elements 
 */
void InitializePETSc(vector<BaseElemTrans*> elements);

/**
 * @brief Calculates the element stiffness matrix. 
 * 
 * @param elements 
 * @param mats 
 */
void CalcElemStiffMatx(vector<BaseElemTrans*> elements, vector<BaseTransport*> mats);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 */
void Assemble(vector<BaseElemTrans*> elements);

/**
 * @brief Reads and initializes Dirichlet BCs.
 * 
 * @param H5File_in 
 */
void InitializeDirichBC(H5IO& H5File_in);

/**
 * @brief Set Dirichlet boundary conditions in the RHS `b` and global stiffness matrix `A`.
 * 
 */
void setDirichBC();

/**
 * @brief Get the number of steps to apply load.
 * 
 * @return int 
 */
int get_nSteps() const;

/**
 * @brief Pass reference of RHS (to solver).
 * 
 * @return Vec& 
 */
Vec& getF();

/**
 * @brief Pass a reference of the solution vector (to solver).
 * 
 * @return Vec& 
 */
Vec& getX();

/**
 * @brief Pass a reference fo the stiffness matrix (to solver).
 * 
 * @return Mat& 
 */
Mat& getK();

/**
 * @brief Calculates the flux. Also Calculates the maps the flux to nodal values.
 * 
 * @param elements 
 */
void CalcFlux(vector<BaseElemTrans*> elements, vector<BaseTransport*> mats);

/**
 * @brief Write nodal values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(vector<BaseElemTrans*> elements, H5IO &H5File_out, const string iStep);

private:

int nElementSets;   /// @brief Number of element sets
int nTotNodes;      /// @brief Total number of nodes.
int nTotDofs;       /// @brief Total number of DOFs.
int nTotElements;   /// @brief Total number of elements.
int nDim;           /// @brief Spatial dimensions of the model.
int nElConDofs;     /// @brief Number of element concentration (temp) dofs.
int nElements;      /// @brief Number of elements per element set.
int nSteps;         /// @brief Number of steps to apply the load.
double dt;          /// @brief Time increment.

T_nodStres nodFlux;        /// @brief Nodal flux.
vector<double> nodCount;     /// @brief Counter for integration points surrounding nodes.

// PETSc ------------------------

const double* globalBuffer;  /// @brief buffer array for PETSc data

// Boundary conditions
PetscInt nPresDofs;            /// @brief number of prescribed dofs.
PetscInt* presDofs = NULL;     /// @brief Array to hold the prescribed dofs.
PetscScalar* presVals = NULL;  /// @brief Array to hold the prescribed values.

Vec F;   /// @brief RHS vector.
Vec x;   /// @brief solution vector.
Mat K;   /// @brief The global conductivity matrix.
Mat M;   /// @brief The global capacitance matrix.

};
#endif