/**
 * @file BaseMechanics.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types for mechanics.
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASEMECHANICS_H
#define BASEMECHANICS_H

#include "Matrix.h"
#include "Materials/BaseMaterial.h"
class BaseMechanics: public BaseMaterial{

public:

BaseMechanics(string isoType, string dimensions): BaseMaterial(isoType, dimensions) {};

/**
 * @brief Returns a stiffness matrix variant.
 * 
 * @return T_DMatx 
 */
virtual T_DMatx getDMatx() const = 0;

};
#endif