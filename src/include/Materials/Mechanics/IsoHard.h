/**
 * @file IsoHard.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief J2 plasticity with isotropic hardening. Inherits from `LinearElastic`.
 * 
 * @date 2024-12-12
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

#ifndef ISOHARD_H
#define ISOHARD_H

#include "LinearElastic.h"
#include "H5IO.h"

/// @cond
struct PlaneStrain {};
struct PlaneStress {};

struct PowerLaw {};
struct Voce {};
struct KME {};
/// @endcond

class IsoHard: public LinearElastic{

public:

enum class AnalysisType{
    PlaneStrain,
    PlaneStress
};

enum class HardeningLaw{
    PowerLaw,
    Voce,
    KME
};

/**
 * @brief Constructor, reads plasticity parameters from hdf5 file.
 * 
 * @param dimensions Dimensions of the model. 
 * @param H5File Input hdf5 file. 
 * @param iSet Element set number. 
 */
IsoHard(string dimensions, H5IO& H5File, int iSet, Logger& logger);

/**
 * @brief Von Mises stress for 3D models. 
 * 
 * @param sig3D Stress vector.
 * @return double 
 */
double Mises3D(const ColVecd6& sig3D);

/**
 * @brief Return mapping algorithm for 3D isotropic hardening plasticity. 
 * 
 * @param deps Strain increment
 * @param sig Stress tensor.
 * @param eps_e Elastic strain.
 * @param eps_p Plastic strain.
 * @param eqp Equivalent plastic strain
 * @param sig_eq Equivalent stress.
 * @param sig_h Hydrostatic stress. 
 * @param rho Normalized dislocation density. 
 * @param sig_old 
 * @param eps_e_old 
 * @param eps_p_old 
 * @param eqp_old 
 * @param iStep Step number.
 */
void ReturnMapping3D(ColVecd6& deps, ColVecd6& sig, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd6& eps_e_old, const ColVecd6& eps_p_old, const double& eps_eq_old, const int iStep);

/**
 * @brief Return mapping algorithm for 2D isotropic hardening plasticity. 
 * 
 * @param deps Strain increment
 * @param sig Stress tensor.
 * @param eps_e Elastic strain.
 * @param eps_p Plastic strain.
 * @param eps_eq Equivalent plastic strain
 * @param sig_eq Equivalent stress.
 * @param sig_h Hydrostatic stress.
 * @param rho Normalized dislocation density.
 * @param eps_e_old 
 * @param eps_p_old 
 * @param eps_eq_old 
 * @param iStep 
 */
void ReturnMapping2D(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep);

void ReturnMapping2D_PFF(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep, const double gPhi_d, const double wp_old, double wp);

/**
 * @brief Returns the stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
T_DMatx getCMatx() const override;

private:

/// @brief Stores the analysis type.
AnalysisType analysis2D; 

/// @brief Stores hardening law.
HardeningLaw hardening;

/// @brief Tolerance for return mapping algorithm.
const double tol = 1e-6;    

/// @brief Maximum number of iterations.
const int max_iter = 20;

/// @brief 3D Identity tensor in Voigt notation.
const ColVecd6 I6 = (ColVecd6() << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0).finished();

/// @brief 2D Identity tensor in Voigt notation.
const ColVecd3 I3 = (ColVecd3() << 1.0, 1.0, 0.0).finished();

/// @brief Plasticity type.
string Platicity;         

/// @brief Initial yield stress. 
double sig_y0;              

 /// @brief Hardening law. 
string hardLaw;     

/// @brief Strength coefficient.       
double K_hard = 0.0;      

/// @brief Strain hardening exponent.
double n_pow = 0.0;  

/// @brief Initial dislocation density
double rho_0 = 0.0;

/// @brief Taylor factor
double M = 0.0;

/// @brief Dislocation interaction constant
double alpha = 0.0;

/// @brief Burgers vector 
double b = 2.5e-10;   

/// @brief Dislocation multiplication coefficient
double k1 = 0.0;

/// @brief Dislocation recovery coefficient
double k2 = 0.0;

/// @brief Saturation dislocation density
double rho_s = 0.0;

/// @brief Constant for KME model
double C_prime = 0.0;

/// @brief Elastoplastic stiffness matrix in Voigt notation.
T_DMatx CMatx_ep;

/**
 * @brief Function to calculate update in plastic stress and hardening modulus
 * 
 * @tparam HardeningLaw type of hardening law
 * @param eqpl Equivalent plastic strain
 * @param syield Current yield stress
 * @param rho Current dislocaiton density
 * @param hard Hardening modulus
 */
template <typename HardeningLaw>
void UHard(const double& eqpl, double& syield, double& rho, double& hard);

/**
 * @brief 3D Return-mapping to handle different hardening laws
 * 
 * @tparam HardeningLaw 
 * @param deps 
 * @param sig 
 * @param eps_e 
 * @param eps_p 
 * @param eps_eq 
 * @param sig_eq 
 * @param eps_e_old 
 * @param eps_p_old 
 * @param eps_eq_old 
 * @param iStep 
 */
template <typename HardeningLaw>
void RM3D(ColVecd6& deps, ColVecd6& sig, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd6& eps_e_old, const ColVecd6& eps_p_old, const double& eps_eq_old, const int iStep);

using RM3DFn = void (IsoHard::*)(ColVecd6&, ColVecd6&, ColVecd6&, ColVecd6&, double&, double&, double&, double&, const ColVecd6&, const ColVecd6&, const double&, const int);

// Function pointer for the selected ReturnMapping variant
RM3DFn selectedRM3D;

// Function to map HardeningLaw to the appropriate ReturnMapping
static RM3DFn selectRM3D(HardeningLaw hardening) {
    switch (hardening) {
        case HardeningLaw::PowerLaw:
            return &IsoHard::RM3D<PowerLaw>;
        case HardeningLaw::Voce:
            return &IsoHard::RM3D<Voce>;
        case HardeningLaw::KME:
            return &IsoHard::RM3D<KME>;
        default:
            throw std::runtime_error("Unsupported hardening law.");
    }
}

/**
 * @brief Template for 2D von Mises stress. Works for plane stress and plane strain.  
 * 
 * @tparam AnalysisType `enmum` class to store analysis type. 
 * @param sig2D Stress vector
 * @return double 
 */
template <typename AnalysisType>
double Mises2D(const ColVecd3& sig2D);

/**
 * @brief Template for 2D hydrostatic stress. Works for plane stress and plane strain.  
 * 
 * @tparam AnalysisType `enmum` class to store analysis type. 
 * @param sig2D Stress vector
 * @return double 
 */
template <typename AnalysisType>
double Shydro2D(const ColVecd3& sig2D);

// 2D Return-mapping to handle different hardening laws and stress state.
template <typename AnalysisType, typename HardeningLaw>
void RM2D(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep);

using RM2DFn = void (IsoHard::*)(ColVecd3&, ColVecd3&, ColVecd3&, ColVecd3&, double&, double&, double&, double&, const ColVecd3&, const ColVecd3&, const double&, const int);

// Function pointer for the selected ReturnMapping variant
RM2DFn selectedRM2D;

// Function to map HardeningLaw to the appropriate ReturnMapping
static RM2DFn selectRM2D(AnalysisType analysis2D, HardeningLaw hardening) {
    
    switch (analysis2D){

        case AnalysisType::PlaneStrain:
            switch (hardening) {
                case HardeningLaw::PowerLaw:
                    return &IsoHard::RM2D<PlaneStrain, PowerLaw>;
                case HardeningLaw::Voce:
                    return &IsoHard::RM2D<PlaneStrain, Voce>;
                case HardeningLaw::KME:
                    return &IsoHard::RM2D<PlaneStrain, Voce>;
                default:
                    throw std::runtime_error("Unsupported hardening law.");
            }
        
        case AnalysisType::PlaneStress:
            switch (hardening) {
                case HardeningLaw::PowerLaw:
                    return &IsoHard::RM2D<PlaneStress, PowerLaw>;
                case HardeningLaw::Voce:
                    return &IsoHard::RM2D<PlaneStress, Voce>;
                case HardeningLaw::KME:
                    return &IsoHard::RM2D<PlaneStress, KME>;
                default:
                    throw std::runtime_error("Unsupported hardening law.");
            }

        default:
            throw std::runtime_error("Unsupported analysis type.");

    }
}

// 2D Return-mapping to handle different hardening laws and stress state for PFF.
template <typename AnalysisType, typename HardeningLaw>
void RM2DPFF(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep, const double gPhi_d, const double wp_old, double wp);

using RM2DFnPFF = void (IsoHard::*)(ColVecd3&, ColVecd3&, ColVecd3&, ColVecd3&, double&, double&, double&, double&, const ColVecd3&, const ColVecd3&, const double&, const int, const double gPhi_d, const double wp_old, double wp);

// Function pointer for the selected ReturnMapping variant
RM2DFnPFF selectedRM2DPFF;

// Function to map HardeningLaw to the appropriate ReturnMapping
static RM2DFnPFF selectRM2DPFF(AnalysisType analysis2D, HardeningLaw hardening) {
    
    switch (analysis2D){

        case AnalysisType::PlaneStrain:
            switch (hardening) {
                case HardeningLaw::PowerLaw:
                    return &IsoHard::RM2DPFF<PlaneStrain, PowerLaw>;
                case HardeningLaw::Voce:
                    return &IsoHard::RM2DPFF<PlaneStrain, Voce>;
                case HardeningLaw::KME:
                    return &IsoHard::RM2DPFF<PlaneStrain, Voce>;
                default:
                    throw std::runtime_error("Unsupported hardening law.");
            }
        
        case AnalysisType::PlaneStress:
            switch (hardening) {
                case HardeningLaw::PowerLaw:
                    return &IsoHard::RM2DPFF<PlaneStress, PowerLaw>;
                case HardeningLaw::Voce:
                    return &IsoHard::RM2DPFF<PlaneStress, Voce>;
                case HardeningLaw::KME:
                    return &IsoHard::RM2DPFF<PlaneStress, KME>;
                default:
                    throw std::runtime_error("Unsupported hardening law.");
            }

        default:
            throw std::runtime_error("Unsupported analysis type.");

    }
}

};

#endif