/**
 * @file PFF.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for phase-field fracture material.
 * 
 * 
 * @date 2025-09-20
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

#ifndef PFF_H
#define PFF_H

#include "BasePFF.h"
#include "H5IO.h"

class PFF: public BasePFF{

public:

/**
 * @brief Constructor, reads elasticity parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
PFF(string dimensions, H5IO& H5File, int iSet, Logger& logger);

private:      


};
#endif