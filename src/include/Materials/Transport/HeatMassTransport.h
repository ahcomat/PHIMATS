/**
 * @file HeatMassTransport.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for heat and mass transport. 
 * @date 2024-06-13
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

#ifndef HEATMASSTRANSPORT_H
#define HEATMASSTRANSPORT_H

#include "BaseTransport.h"
#include "H5IO.h"

class HeatMassTransport: public BaseTransport{

public:

/**
 * @brief Constructor, reads heat/mass transport parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
HeatMassTransport(string dimensions, H5IO& H5File, int iSet);

/**
 * @brief Returns the diffusivity (conductivity) matrix.
 * 
 * @return T_DMatx 
 */
T_DMatx getKMatx() const override;

/**
 * @brief Returns the diffusive (heat) capacity.
 * 
 * @return double 
 */
double getCapacity() const;

private:

double s;           /// @brief (diffusive) capacity
T_DMatx KMatx;      /// @brief Variant for storing the diffusivity (conductivity) matrix

};
#endif