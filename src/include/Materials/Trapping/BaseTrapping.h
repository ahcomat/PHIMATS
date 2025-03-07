/**
 * @file BaseTrapping.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types for heat and mass transport.
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

#ifndef BaseTrapping_H
#define BaseTrapping_H

#include "Matrix.h"
#include "Materials/BaseMaterial.h"

class BaseTrapping: public BaseMaterial{

public:

BaseTrapping(string dimensions, Logger& logger): BaseMaterial(dimensions, logger) {};
virtual ~BaseTrapping() override {};

/**
 * @brief Calculates the spatially varying diffusivity matrix.
 * 
 * @param phi Spatial field parameter.
 * @param T Temperature [K].
 * @return T_DMatx Diffusivity variant.
 */
virtual T_DMatx CalcDMatx(const double phi, const double T) = 0;

protected:

const double R = 8.31446261815324; /// @brief Universal gas constant [J/mol.K]

};
#endif