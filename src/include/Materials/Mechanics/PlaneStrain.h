/**
 * @file PlaneStrain.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Elastic plane-strain model.
 * 
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PLANESTRAIN_H
#define PLANESTRAIN_H

#include "Matrix.h"
#include "BaseMechanics.h"
#include "H5IO.h"

class PlaneStrain: public BaseMechanics{

public:

/**
 * @brief Constructor, reads elastic plane-strain parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
PlaneStrain(H5IO &H5File, int iSet, string isoType="Isotropic");

/**
 * @brief Destructor.
 */
~PlaneStrain() override {}

/**
 * @brief Returns the 2D stiffness matrix in Voigt notation.
 */
T_DMatx getDMatx() override;

private:

/**
 * @brief The 2D elastic stiffness matrix in Voigt notation.
 * 
 */
Matd3x3 DMatx;

};
#endif