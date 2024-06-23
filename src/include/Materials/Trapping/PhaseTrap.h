/**
 * @file PhaseTrap.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for phase trapping. 
 * @date 2024-06-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PHASETRAP_H
#define PHASETRAP_H

#include "BaseTrapping.h"
#include "H5IO.h"

class PhaseTrap: public BaseTrapping{

public:

/**
 * @brief Constructor, reads heat/mass transport parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
PhaseTrap(string dimensions, H5IO &H5File, int iSet, string isoType="Isotropic");

/**
 * @brief Calculates the phase-field dependent diffusivity matrix.
 * 
 */
T_DMatx CalcKMatx(const double phiL, const double phiT, const double T);

/**
 * @brief Calculates the phase-field dependent trapping matrix.
 * 
 */
T_DMatx CalcTMatx(const double phiL, const double phiT, const double T);

private:

/// @ brief diffusivity parameters @todo have to modify
double D0x1, D0y1, DQx1, DQy1, D0x2, D0y2, DQx2, DQy2;  

/// @brief Phase trapping parameter  
double zeta;      

};
#endif