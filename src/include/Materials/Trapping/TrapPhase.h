/**
 * @file TrapPhase.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for phase trapping. 
 * @date 2024-07-17
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

#ifndef TRAPPHASE_H
#define TRAPPHASE_H

#include "BaseTrapping.h"
#include "H5IO.h"

class TrapPhase: public BaseTrapping{

public:

/**
 * @brief Constructor, reads GB diffusion and trapping parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
TrapPhase(string dimensions, H5IO &H5File, int iSet, Logger& logger);

T_DMatx CalcDMatx(const double phi, const double T) override;

double get_zeta_j() const {return zeta_j;};

double get_zeta_jj() const {return zeta_jj;};

double get_zeta_ij() const {return zeta_ij;};

double get_zeta_ii() const {return zeta_ii;};


private:

/// @ brief diffusivity parameters @todo have to modify
double D0x1, D0y1, DQx1, DQy1, D0x2, D0y2, DQx2, DQy2;  

double zeta_j;       /// @brief Phase phi_j trapping parameter  
double zeta_jj;      /// @brief gPhi_jj interface trapping parameter  
double zeta_ij;      /// @brief gPhi_ij trapping parameter  
double zeta_ii;      /// @brief gPhi_ii trapping parameter  

};
#endif