/**
 * @file MechTrap.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for hydrostatic stress trapping. 
 * @date 2024-06-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef MECHTRAP_H
#define MECHTRAP_H

#include "BaseTrapping.h"
#include "H5IO.h"

class MechTrap: public BaseTrapping{

public:

/**
 * @brief Constructor, reads mechanical trapping parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
MechTrap(string dimensions, H5IO &H5File, int iSet, string isoType="Isotropic");

/**
 * @brief Calculates the phase-field dependent diffusivity matrix.
 * 
 */
T_DMatx CalcKMatx(const double T);

/**
 * @brief Calculates the phase-field dependent trapping matrix.
 * 
 */
T_DMatx CalcTMatx(const double T);

private:

/// @ brief diffusivity parameters @todo have to modify
double D0x1, D0y1, DQx1, DQy1; 

/// @brief Phase trapping parameter  
double Vh; 

// /// @brief Variant for storing the diffusivity (conductivity) matrix
// T_DMatx KMatx; 

// /// @brief Variant for storing the trapping matrix D*zeta/(RT)
// T_DMatx TMatx;      

};
#endif