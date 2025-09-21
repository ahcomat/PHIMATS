/**
 * @file BaseMechanics.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types for mechanics.
 * @date 2024-05-18
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

#ifndef BASEMECHANICS_H
#define BASEMECHANICS_H

#include "Matrix.h"
#include "Materials/BaseMaterial.h"

class BaseMechanics: public BaseMaterial{

public:

BaseMechanics(string dimensions, Logger& logger): BaseMaterial(dimensions, logger) {};

/**
 * @brief Get the Lambda constant.
 * 
 * @return double 
 */
double getLambda() const { return ho; }

/**
 * @brief Get the Gmod constant.
 * 
 * @return double 
 */
double getGmod() const { return uo; }

/**
 * @brief Returns a stiffness matrix variant.
 * 
 * @return T_DMatx 
 */
virtual T_DMatx getCMatx() const = 0;

protected:

/// @brief Isotropic elasticity parameters.
double uo = 0.0;
double ho = 0.0;
double Emod = 0.0;
double nu = 0.0;

/// @brief Cubic elasticity parameters.
double C11 = 0.0;
double C12 = 0.0;
double C44 = 0.0;

};
#endif