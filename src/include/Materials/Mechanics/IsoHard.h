/**
 * @file IsoHard.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief J2 plasticity with isotropic hardening.
 * 
 * @date 2024-12-12
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

#ifndef ISOHARD_H
#define ISOHARD_H

#include "LinearElastic.h"
#include "H5IO.h"

struct PlaneStrain {};
struct PlaneStress {};

class IsoHard: public LinearElastic{

public:

enum class AnalysisType{
    PlaneStrain,
    PlaneStress
};

/**
 * @brief Constructor, reads plasticity parameters from hdf5 file.
 * 
 * @param dimensions Dimensions of the model. 
 * @param H5File Input hdf5 file. 
 * @param iSet Element set number. 
 */
IsoHard(string dimensions, H5IO& H5File, int iSet);

/**
 * @brief Power-law hardening
 * 
 * @param eps_eq 
 * @return double 
 */
double R_pow(const double& eps_eq);

/**
 * @brief Derivative of power-law hardening
 * 
 * @param eps_eq 
 * @return double 
 */
double dR_pow(const double& eps_eq);

/**
 * @brief Von Mises stress for 3D models. 
 * 
 * @param sig3D Stress vector.
 * @return double 
 */
double Mises3D(const ColVecd6& sig3D);

/**
 * @brief Return mapping algorithm for isotropic hardening plasticity. 
 * 
 * @param sig Stress tensor.
 * @param eps Total strain tensor.
 * @param eps_e Elastic strain.
 * @param eps_p Plastic strain.
 * @param eqp Equivalent plastic strain.
 * @param sig_eq Equivalent stress.
 * @param iStep Step number. 
 */
void ReturnMapping3D(ColVecd6& sig, ColVecd6& eps, ColVecd6& eps_e, ColVecd6& eps_p, double& eqp, double& sig_eq, const int iStep);

void ReturnMapping2D(ColVecd3& sig, ColVecd3& eps, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, const int iStep);

/**
 * @brief Returns the stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
T_DMatx getDMatx() const override;

private:

/// @brief Stores the analysis type
AnalysisType analysis2D; 

/// @brief Tolerance for return mapping algorithm.
const double tol = 1e-6;    

/// @brief Maximum number of iterations.
const int max_iter = 10;

/// @brief 3D Identity tensor in Voigt notation.
const ColVecd6 I6 = (ColVecd6() << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0).finished();

/// @brief 2D Identity tensor in Voigt notation.
const ColVecd3 I3 = (ColVecd3() << 1.0, 1.0, 0.0).finished();

/// @brief Plasticity type.
string Platicity;         

/// @brief Initial yield stress. 
double sig_y0;              

 /// @brief Hardening law. 
string HardLaw;     

/// @brief Strength coefficient.       
double K_hard = 0.0;      

/// @brief Strain hardening exponent.
double n_pow = 0.0;  

/// @brief Elastioplastic stiffness matrix in Voigt notation.
T_DMatx DMatx_ep;        

/**
 * @brief Template for 2D von Mises stress. Works for plane stress and plane strain.  
 * 
 * @tparam AnalysisType `enmum` class to store analysis type. 
 * @param sig2D Stress vector
 * @return double 
 */
template <typename AnalysisType>
double Mises2D(const ColVecd3& sig2D);

template <typename AnalysisType>
void ReturnMapping2D(ColVecd3& sig, ColVecd3& eps, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, int iStep);

};

#endif