/**
 * @file BasePFFs.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types for phase-field fracture.
 * @date 2024-09-20
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

#ifndef BASEPFF_H
#define BASEPFF_H

#include "Matrix.h"
#include "Materials/BaseMaterial.h"

class BasePFF: public BaseMaterial{

public:

BasePFF(string dimensions, Logger& logger): BaseMaterial(dimensions, logger) {};

/**
 * @brief Get critical work density.
 * 
 * @return double 
 */
double get_wc() const {return wc;}

/**
 * @brief Get const_ell.
 * 
 * @return double 
 */
double get_const_ell() const {return const_ell;}

protected:

/// @brief Critial work density.
double wc = 0.0;

/// @brief PFF internal length scale.
double const_ell = 0.0;

};
#endif