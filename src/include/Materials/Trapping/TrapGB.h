/**
 * @file TrapGB.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for grain boundary trapping. 
 * @date 2024-06-21
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

#ifndef TRAPGB_H
#define TRAPGB_H

#include "BaseTrapping.h"
#include "H5IO.h"

class TrapGB: public BaseTrapping{

public:

/**
 * @brief Constructor, reads GB diffusion and trapping parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
TrapGB(string dimensions, H5IO& H5File, int iSet, Logger& logger);

T_DMatx CalcDMatx(const double phi, const double T) override;

/**
 * @brief Calculates the phase-field dependent trapping matrix.
 * 
 */
T_DMatx CalcTMatx(const double gPhi, const double T);

private:

/// @ brief diffusivity parameters @todo have to modify
double D0x1, D0y1, DQx1, DQy1, D0x2, D0y2, DQx2, DQy2;  

/// @brief GB trapping parameter  
double zeta_GB;      

};
#endif