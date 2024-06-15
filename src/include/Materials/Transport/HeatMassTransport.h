/**
 * @file HeatMassTransport.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A model for heat and mass transport. 
 * @date 2024-06-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef HEATMASSTRANSPORT_H
#define HEATMASSTRANSPORT_H

#include "BaseTransport.h"
#include "H5IO.h"

class HeatMassTransport: public BaseTransport{

public:

/**
 * @brief Constructor, reads heat/mass transport parameters from hdf5 file.
 * 
 * @param H5File Input file.
 * @param matType Material isotropy.
 */
HeatMassTransport(string dimensions, H5IO &H5File, int iSet, string isoType="Isotropic");

/**
 * @brief Returns the diffusivity (conductivity) matrix.
 * 
 * @return T_DMatx 
 */
T_DMatx getKMatx() const override;

/**
 * @brief Returns the diffusive (heat) capacity matrix.
 * 
 * @return T_DMatx 
 */
T_DMatx getKMatx() const;

private:

double rho;         /// @brief Mass density 
double c;           /// @brief Specific heat (diffusive) capacity
T_DMatx KMatx;      /// @brief Variant for storing the diffusivity (conductivity) matrix

};
#endif