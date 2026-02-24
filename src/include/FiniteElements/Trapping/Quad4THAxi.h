/**
 * @file Quad4THAxi.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Quad4THAxi elements for axisymmetric full-fileld hydrogen trapping. 
 * @date 2026-02-24
 * 
 * @copyright Copyright (C) 2026 Abdelrahman Hussein
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

#include "FiniteElements/Trapping/Quad4TH.h"

#ifndef QUAD4THAXI_H
#define QUAD4THAXI_H

class Quad4THAxi: public Quad4TH{

public:

Quad4THAxi(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, Logger& logger, H5IO* H5File_rve = nullptr);   

void CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const RowVecd4& shFunc, 
                       const double& wt, double& intVol, Matd2x4& cartDeriv) override;

};
#endif