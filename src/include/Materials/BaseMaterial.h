/**
 * @file BaseMaterial.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types.
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 *
 */

#ifndef BASEMATERIAL_H
#define BASEMATERIAL_H

#include <string>

using namespace std;

class BaseMaterial{

public:

BaseMaterial(string dimensions): dims(dimensions) {};
virtual ~BaseMaterial() = default;

string getDims(){ return dims; }

protected:

const string dims;      /// @brief Dimensions of the material model.

};
#endif