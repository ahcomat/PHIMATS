/**
 * @file LinearElastic.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for linear elasticity. Supports `3D`, `PlaneStress` and `PlaneStrain`. Also works for `Isotropic` and `Cubic` anisotropy.
 * 
 * The class supports the following isotropies:
 * `Isotropic` 
 * `Cubic`
 * 
 * And elsticity types:
 * `3D`
 * `PlaneStrain`
 * 'PlaneStress'
 * 
 * @date 2024-05-20
 * 
 * @copyright Copyright (C) 2025 Abdelrahman Hussein
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  
 */

#ifndef LINEARELASTIC_H
#define LINEARELASTIC_H

#include "BaseMechanics.h"
#include "H5IO.h"

class LinearElastic: public BaseMechanics{

public:

/**
 * @brief Constructor, reads elasticity parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
LinearElastic(string dimensions, H5IO& H5File, int iSet, Logger& logger);

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
double Emod = 0.0;
double nu = 0.0;

/// @brief Cubic elasticity parameters.
double C11 = 0.0;
double C12 = 0.0;
double C44 = 0.0;

/// @brief Elastic stiffness matrix in Voigt notation.
T_DMatx DMatx_e;

};
#endif