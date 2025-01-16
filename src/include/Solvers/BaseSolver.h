/**
 * @file BaseSolver.cxx
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A wrapper for PETSc linear solver `ksp` mainly for linear elasticity.
 * @date 2024-05-28
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

#ifndef BASESOLVER_H
#define BASESOLVER_H

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
class BaseSolver{

public:

virtual ~BaseSolver() = default;

protected:

KSP ksp;        /// @brief `KSP` object.
PC pc;          /// @brief Pre-conditioner.

};
#endif