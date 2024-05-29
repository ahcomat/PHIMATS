/**
 * @file MechModel.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A class based on mechanical models to interface with PETSc based global
 *        stiffness matrix, solution vector and boundary conditions. It calls
 *        `BaseElemMech::CalcElemStiffMatx` to build the local stiffness matrix. 
 *        It also manages the output by writing to H5file_out. 
 *          
 * @date 2024-05-24
 * 
 * @copyright Copyright (c) 2024
 * 
 * @todo - Use a vector of `BaseElemMech` for multi-material models.
 * - Consider inheritance.
 * 
 * Updates (when, what and who)
 * 
 */

#ifndef MECHMODEL_H
#define MECHMODEL_H

#include"petsc.h"
#include"H5IO.h"
#include"FiniteElements/Mechanics/BaseElemMech.h"

class MechModel
{

public:

MechModel(BaseElemMech* elements);
~MechModel();

/**
 * @brief Initializes and preallocates the RHS `b`, solution `x` and 
 *        stiffness matrix `A.
 * 
 * @param elements 
 * @param H5File_in 
 */
void InitializePETSc(BaseElemMech* elements);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 */
void Assemble(BaseElemMech* elements);

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
 * @brief Calculates the Fint, strains and stresses. Also evaluates the stress nodal values 
 *        if `nodStresFlag=true`.
 * 
 * @param elements 
 * @param nodStresFlag 
 */
void CalcStres(BaseElemMech* elements, T_DMatx DMatx, bool nodStresFlag=false);

/**
 * @brief Write nodal values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(BaseElemMech* elements, H5IO &H5File_out);

private:

int nTotDof;        /// @brief Total number of DOFs.
int nElDispDofs;    /// @brief Number of element displacement dofs.
int nElements;      /// @brief Total number of elements.
int nDim;           /// @brief Spatial dimensions of the model.

// PETSc ------------------------

const double* globalBuffer;

// Boundary conditions
PetscInt nPresDofs;     /// @brief number of prescribed dofs.
PetscInt  *presDofs;    /// @brief Array to hold the prescribed dofs.
PetscScalar  *presVals; /// @brief Array to hold the prescribed values.

Vec b;   /// @brief RHS vector.
Vec x;   /// @brief solution vector.
Mat A;   /// @brief The global coefficient (stiffness) matrix.

};
#endif