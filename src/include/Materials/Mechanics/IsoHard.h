/**
 * @file IsoHard.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief J2 plasticity with isotropic hardening.
 * 

 * 
 * @date 2024-12-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ISOHARD_H
#define ISOHARD_H

#include "Elastic.h"
#include "H5IO.h"

class IsoHard: public Elastic{

public:

/**
 * @brief Constructor, reads plasticity parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
IsoHard(string dimensions, H5IO& H5File, int iSet);

/**
 * @brief Power-law hardening
 * 
 * @param eps_eq 
 * @return double 
 */
double R_pow(double eps_eq);

/**
 * @brief Derivative of power-law hardening
 * 
 * @param eps_eq 
 * @return double 
 */
double dR_pow(double eps_eq);

/**
 * @brief Von Mises stress for 3D models. 
 * 
 * @param sig3D 
 * @return double 
 */
double Mises3D(ColVecd6& sig3D);

// /**
//  * @brief Equivalent plastic strain for 3D models. 
//  * 
//  * @param eps3D 
//  * @return double 
//  */
// double PEEQ3D(ColVecd6& eps3D);

void ReturnMapping2D();

/**
 * @brief Returns the stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
T_DMatx getDMatx() const override;

private:

const double tol = 1e-6;    /// @brief Tolerance for return mapping algorithm.
const int max_iter = 10;    /// @brief Maximum number of iterations. 

string Platicity;           /// @brief Plasticity type.

string HardLaw;             /// @brief Hardening law. 
double K_hard = 0.0;        /// @brief Strength coefficient.
double n_pow = 0.0;         /// @brief Strain hardening exponent.

T_DMatx DMatx_e;           /// @brief Elastic stiffness matrix in Voigt notation.
T_DMatx DMatx_ep;           /// @brief The elastioplastic stiffness matrix in Voigt notation.

};
#endif