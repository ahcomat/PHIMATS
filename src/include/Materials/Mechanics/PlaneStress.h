/**
 * @file PlaneStress.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Elastic plane-stress model.
 * 
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 *
 */

#ifndef PLANESTRESS_H
#define PLANESTRESS_H

#include "BaseMechanics.h"
#include "H5IO.h"

class PlaneStress: public BaseMechanics{

public:

/**
 * @brief Constructor, reads elastic plane-strain parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
PlaneStress(H5IO &H5File, int iSet, string isoType="Isotropic");

/**
 * @brief Destructor.
 */
~PlaneStress() override {}

/**
 * @brief Returns the 2D stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
T_DMatx getDMatx() override;

private:

Matd3x3 DMatx;      /// @brief The 2D elastic stiffness matrix for plane-stress in Voigt notation.

};
#endif