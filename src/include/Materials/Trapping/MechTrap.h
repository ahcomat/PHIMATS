/**
 * @file MechTrap.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for hydrostatic stress trapping. 
 * @date 2024-06-23
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

#ifndef MECHTRAP_H
#define MECHTRAP_H

#include "BaseTrapping.h"
#include "H5IO.h"

class MechTrap: public BaseTrapping{

public:

/**
 * @brief Constructor, reads mechanical trapping parameters from HDF5 file.
 * 
 * @param dimensions Dimensions of the model. Options ["2D" or "3D"].
 * @param H5File Input HDF5 file.
 * @param iSet Set number. 
 * @param logger Logger object. 
 */
MechTrap(string dimensions, H5IO& H5File, int iSet, Logger& logger);

T_DMatx CalcDMatx(const double phi, const double T) override;

/**
 * @brief Get the trapping capacity `s`.
 * 
 * @return double
 */
double get_s() const { return s; };

/**
 * @brief Get the partial molar volume `Vh`.
 * 
 * @return double
 */
    double get_Vh() const { return Vh; };

/**
 * @brief Get trapping parameter `zeta_rho`.
 * 
 * @return double
 */
double get_zeta_rho() const { return zeta_rho; };

/**
 * @brief Crack source/sink term.
 * 
 * @return double
 */
double get_Zd() const { return Zd; };

private:

/// @brief Anisotropic diffusivity parameters
double D0x, D0y, D0z, DQx, DQy, DQz; 

/// @brief Partial molar volume for hydrogen [m]
double Vh = 0;       

/// @brief Dislocation to lattice diffusivity ratio  
double m = 0;

/// @brief Trapping capacity  
double s = 0;

/// @brief Dislocation trapping parameter 
double zeta_rho = 0;

/// @brief Source/sink term at the damaged region  
double Zd = 0;

};
#endif