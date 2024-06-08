/**
 * @file MechModel.h
 * @brief A class based on mechanical models to interface with PETSc based global
 * stiffness matrix `A`, solution vector `x` and RHS `a`.
 * 
 * @details Main functions:
 * - Calls `BaseElemMech::CalcElemStiffMatx` to build the local stiffness matrix.
 * - Assembles the global stiffness matrix.
 * - Initializes and applies boundary conditions.
 * - Manages the output by writing to H5file_out.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2024-05-24
 * 
 * @copyright Copyright (c) 2024
 * 
 * @todo - Use a vector of `BaseElemMech` for multi-material models.
 * - Consider inheritance.
 *  
 */

#ifndef MECHMODEL_H
#define MECHMODEL_H

#include "petsc.h"
#include "H5IO.h"
#include "FiniteElements/Mechanics/BaseElemMech.h"
#include "Materials/Mechanics/BaseMechanics.h"

class MechModel{

public:

MechModel(vector<BaseElemMech*> elements, H5IO& H5File_in);
~MechModel();

/**
 * @brief Initializes and preallocates the RHS `b`, solution `x` and 
 *        stiffness matrix `A.
 * 
 * @param elements 
 * @param H5File_in 
 */
void InitializePETSc(vector<BaseElemMech*> elements);

/**
 * @brief Calculates the element stiffness matrix. 
 * 
 * @param elements 
 * @param mats 
 */
void CalcElemStiffMatx(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 */
void Assemble(vector<BaseElemMech*> elements);

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
Vec& getB();

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
Mat& getA();

/**
 * @brief Calculates the Fint, strains and stresses. Also Calculates the stress nodal values 
 *        if `nodStresFlag=true`.
 * 
 * @param elements 
 * @param nodStresFlag 
 */
void CalcStres(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats);

/**
 * @brief Write nodal values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(vector<BaseElemMech*> elements, H5IO &H5File_out);

private:

int nElementSets;   /// @brief Number of element sets
int nTotNodes;       /// @brief Total number of DOFs.
int nTotDofs;       /// @brief Total number of DOFs.
int nTotElements;   /// @brief Total number of elements.
int nDim;           /// @brief Spatial dimensions of the model.
int nElDispDofs;    /// @brief Number of element displacement dofs.
int nElements;      /// @brief Number of elements per element set. 
// int nSteps;      /// @brief Number of steps to apply the load.

T_nodStres nodStres;      /// @brief Nodal stress.
T_nodStres nodStran;      /// @brief Nodal strain.
vector<int> nodCount;     /// @brief Counter for integration points surrounding nodes.

// PETSc ------------------------

const double* globalBuffer;  /// @brief buffer array for PETSc data
double* Fint = NULL;         /// @brief For calculating the internal force vector.

// Boundary conditions
PetscInt nPresDofs;            /// @brief number of prescribed dofs.
PetscInt* presDofs = NULL;     /// @brief Array to hold the prescribed dofs.
PetscScalar* presVals = NULL;  /// @brief Array to hold the prescribed values.

Vec b;   /// @brief RHS vector.
Vec x;   /// @brief solution vector.
Mat A;   /// @brief The global coefficient (stiffness) matrix.

};
#endif