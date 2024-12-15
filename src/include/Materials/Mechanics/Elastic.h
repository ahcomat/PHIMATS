/**
 * @file Elastic.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Elastic tensor in Voigt notation.
 * 
 * The class supports the following isotropies:
 * `Isotropic` 
 * `Cubic`
 * 
 * And Elasticities:
 * `3D`
 * `PlaneStrain`
 * 'PlaneStress'
 * 
 * @date 2024-05-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ELASTIC_H
#define ELASTIC_H

#include "BaseMechanics.h"
#include "H5IO.h"

class Elastic: public BaseMechanics{

public:

/**
 * @brief Constructor, reads elasticity parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
Elastic(string dimensions, H5IO& H5File, int iSet);

void InitializeIsoElasticityMatrix(const string& elasticity, double Emod, double nu, double ho, double uo);

void InitializeCubicElasticityMatrix(const string& elasticity, double C11, double C12, double C44);

/**
 * @brief Returns the 3D stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
T_DMatx getDMatx() const override;

private:

T_DMatx DMatx;      /// @brief The 3D elastic stiffness matrix in Voigt notation.

};
#endif