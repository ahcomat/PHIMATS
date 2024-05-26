/**
 * @file BaseMaterial.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types.
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 * 
 * Updates (when, what and who)
 * 
 */

#ifndef BASEMATERIAL_H
#define BASEMATERIAL_H

#include <string>

using namespace std;

class BaseMaterial{

public:

BaseMaterial(string isoType, string dimensions): isotropy(isoType), dims(dimensions) {};
virtual ~BaseMaterial() = default;

string getIsotropy(){ return isotropy; }
string getDims(){ return dims; }

protected:
/**
 * @brief Isotropy type of the material. 
 * 
 * Available types are: `Isotropic` and `Cubic`.
 * 
 */
const string isotropy; 

/**
 * @brief Dimensions of the material tensor.
 * 
 */
const string dims;

};
#endif