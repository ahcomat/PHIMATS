/**
 * @file Elastic3D.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief 3D Elastic tensor in Voigt notation.
 * 
 * The class supports the following isotropies:
 * `Isotropic` 
 * `Cubic`
 * 
 * @date 2024-05-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ELASTIC3D_H
#define ELASTIC3D_H

#include "BaseMechanics.h"
#include "H5IO.h"

class Elastic3D: public BaseMechanics{

public:

/**
 * @brief Constructor, reads elastic parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
Elastic3D(H5IO &H5File, int iSet, string isoType="Isotropic");

/**
 * @brief Returns the 3D stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
T_DMatx getDMatx() override;

private:

/**
 * @brief The 3D elastic stiffness matrix in Voigt notation.
 */
Matd6x6 DMatx;

};
#endif