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
 * @todo - Consider inheritance.
 *  
 */

#ifndef MECHMODEL_H
#define MECHMODEL_H

#include <petscvec.h>
#include <petscmat.h>
#include "BaseModel.h"
#include "H5IO.h"
#include "FiniteElements/Mechanics/BaseElemMech.h"
#include "Materials/Mechanics/BaseMechanics.h"

class MechModel: public BaseModel{

public:

MechModel(vector<BaseElemMech*> elements, H5IO& H5File_in);
~MechModel();

/**
 * @brief Set nodal quantities `nodStres`, `nodStran`, `nodCount` and `Fint` to zeros.
 * 
 */
void setZeroNodVals();

/**
 * @brief Initializes and preallocates the RHS `b`, solution `x` and 
 *        stiffness matrix `A.
 * 
 * @param elements 
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
 * @brief Calculates the Fint, strains and stresses. Also Calculates the stress nodal values.
 * 
 * @param elements 
 */
void CalcStres(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats);

/**
 * @brief Write nodal values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(vector<BaseElemMech*> elements, H5IO &H5File_out, const string iStep);

private:

const double tol = 1e-6;    /// @brief Tolerance for Newton-Raphson iterations.
const int max_iter = 10;    /// @brief Maximum number of iterations.
const int NR_freq = 3;      /// @brief Frequency of updating the stiffness matrix.

int nElDispDofs;    /// @brief Number of element displacement dofs.

T_nodStres nodStres;      /// @brief Nodal stress.
T_nodStres nodStran;      /// @brief Nodal strain.
vector<int> nodCount;     /// @brief Counter for integration points surrounding nodes.

// PETSc ------------------------

double* Fint = NULL;         /// @brief For calculating the internal force vector.

/// @brief RHS vector.
Vec b; 

/// @brief solution vector.         
Vec vecDisp;

/// @brief solution vector.         
Vec vecFint;

/// @brief The global coefficient (stiffness) matrix.
Mat matA;       

};
#endif