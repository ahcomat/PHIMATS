/**
 * @file LinearTransport.cxx
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A wrapper for PETSc linear solver `ksp` mainly for linear heat and mass transport.
 * @date 2024-06-18
 * 
 * @copyright Copyright (C) 2024 Abdelrahman Hussein
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

#ifndef LINEARTRANSPORT_H
#define LINEARTRANSPORT_H

#include "BaseSolver.h"
#include "Logger.h"

#include <string>

using namespace std;

class LinearTransport: public BaseSolver{

public:

LinearTransport(Mat &A, Logger& logger, string solverType="DIRECT");

~LinearTransport() override;

void UpdateKSP(Mat &A);

/**
 * @brief Solve the linear system `Ax=b`.
 * 
 * @param A 
 * @param x 
 * @param b 
 */
void SolveTransport(Vec &x, Vec &F);

};
#endif