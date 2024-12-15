/**
 * @file TrapGB.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for grain boundary trapping. 
 * @date 2024-06-21
 * 
 * @copyright Copyright (c) 2024
 * s
 */

#ifndef TRAPGB_H
#define TRAPGB_H

#include "BaseTrapping.h"
#include "H5IO.h"

class TrapGB: public BaseTrapping{

public:

/**
 * @brief Constructor, reads GB diffusion and trapping parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
TrapGB(string dimensions, H5IO& H5File, int iSet);

/**
 * @brief Calculates the phase-field dependent diffusivity matrix.
 * 
 */
T_DMatx CalcDMatx(const double gPhi, const double T);

/**
 * @brief Calculates the phase-field dependent trapping matrix.
 * 
 */
T_DMatx CalcTMatx(const double gPhi, const double T);

private:

/// @ brief diffusivity parameters @todo have to modify
double D0x1, D0y1, DQx1, DQy1, D0x2, D0y2, DQx2, DQy2;  

/// @brief GB trapping parameter  
double kappa_GB;      

};
#endif