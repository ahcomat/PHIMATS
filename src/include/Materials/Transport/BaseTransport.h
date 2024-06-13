/**
 * @file BaseTransport.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types for heat and mass transport.
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASETRANSPORT_H
#define BASETRANSPORT_H

#include "Matrix.h"
#include "Materials/BaseMaterial.h"

class BaseTransport: public BaseMaterial{

public:

BaseTransport(string isoType, string dimensions): BaseMaterial(isoType, dimensions) {};
virtual ~BaseTransport() override {};

/**
 * @brief Returns the 3D stiffness matrix in Voigt notation.
 * 
 * @return T_DMatx 
 */
virtual T_DMatx getKMatx() const = 0;

};
#endif