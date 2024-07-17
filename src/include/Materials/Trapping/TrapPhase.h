/**
 * @file TrapPhase.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for phase trapping. 
 * @date 2024-07-17
 * 
 * @copyright Copyright (c) 2024
 * s
 */

#ifndef TRAPPHASE_H
#define TRAPPHASE_H

#include "BaseTrapping.h"
#include "H5IO.h"

class TrapPhase: public BaseTrapping{

public:

/**
 * @brief Constructor, reads GB diffusion and trapping parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
TrapPhase(string dimensions, H5IO &H5File, int iSet, string isoType="Isotropic");

/**
 * @brief Calculates the phase-field dependent diffusivity matrix.
 * 
 */
T_DMatx CalcDMatx(const double phi, const double T);

double get_zeta_M() const {return zeta_M;};

double get_zeta_MM() const {return zeta_MM;};

double get_zeta_fM() const {return zeta_fM;};

double get_zeta_ff() const {return zeta_ff;};


private:

/// @ brief diffusivity parameters @todo have to modify
double D0x1, D0y1, DQx1, DQy1, D0x2, D0y2, DQx2, DQy2;  

double zeta_M;       /// @brief martensite phase trapping parameter  
double zeta_MM;      /// @brief martensite-martensite interface trapping parameter  
double zeta_fM;      /// @brief ferrite-martensite trapping parameter  
double zeta_ff;      /// @brief ferrite-ferrite trapping parameter  

};
#endif