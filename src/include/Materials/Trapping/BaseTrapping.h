/**
 * @file BaseTrapping.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types for heat and mass transport.
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BaseTrapping_H
#define BaseTrapping_H

#include "Matrix.h"
#include "Materials/BaseMaterial.h"

class BaseTrapping: public BaseMaterial{

public:

BaseTrapping(string isoType, string dimensions): BaseMaterial(isoType, dimensions) {};
virtual ~BaseTrapping() override {};

protected:

const double R = 8.31446261815324; /// @brief Universal gas constant [J/mol.K]

};
#endif