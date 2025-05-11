/**
 * @file TrapHLGB.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for grain boundary interacion taking into account the grain boundary character. 
 * @date 2025-05-11
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

#ifndef TRAPHLGB_H
#define TRAPHLGB_H

#include "BaseTrapping.h"
#include "H5IO.h"

class TrapHLGB: public BaseTrapping{

public:

/**
 * @brief Constructor, reads GB diffusion and interaction parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
TrapHLGB(string dimensions, H5IO& H5File, int iSet, Logger& logger);

/**
 * @brief Function for calculating the diffusivity
 * 
 * @param phi 
 * @param T 
 * @return T_DMatx 
 */
T_DMatx CalcDMatx(const double phi, const double T) override;

double get_zeta_HAGB() const {return zeta_HAGB;};

double get_zeta_LAGB() const {return zeta_LAGB;};

private:

/// @ brief diffusivity parameters 
double D0x1, D0y1, DQx1, DQy1;  

/// @brief HAGB interaction parameter  
double zeta_HAGB;  

/// @brief LAGB interaction parameter  
double zeta_LAGB; 

};
#endif