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

/**
 * @brief Initialize the isotropic elastic stiffness matrix in Voigt notation
 * 
 * @param elasticity 
 * @param Emod 
 * @param nu 
 * @param ho 
 * @param uo 
 */
void InitializeIsoElasticityMatrix(const string& elasticity, double Emod, double nu, double ho, double uo);

/**
 * @brief Initialize the cubic elastic stiffness matrix in Voigt notation
 * 
 * @param elasticity 
 * @param C11 
 * @param C12 
 * @param C44 
 */
void InitializeCubicElasticityMatrix(const string& elasticity, double C11, double C12, double C44);

/**
 * @brief Returns the stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
T_DMatx getDMatx() const override;

private:

protected:

/// @brief Analysis type: `3D`, `PlaneStrain` or `PlaneStress`
string analysisType;    

/// @brief Isotropic elasticity parameters.
double uo = 0.0;
double ho = 0.0;

/// @brief Cubic elasticity parameters.
double C11 = 0.0;
double C12 = 0.0;
double C44 = 0.0;

/// @brief Elastic stiffness matrix in Voigt notation.
T_DMatx DMatx_e;

};
#endif