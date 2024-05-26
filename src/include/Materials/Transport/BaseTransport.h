/**
 * @file BaseTransport.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for material types for heat and mass transport.
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 * 
 * Updates (when, what and who)
 * 
 */

#ifndef BASETRANSPORT_H
#define BASETRANSPORT_H

#include "Matrix.h"
#include "Materials/BaseMaterial.h"

using namespace std;

class BaseTransport: public BaseMaterial{

public:

BaseTransport(string isoType, string dimensions): BaseMaterial(isoType, dimensions) {};
virtual ~BaseTransport() override {};

/*
    NOTE: Don't use pure virtual functions or some derived classed will 
    be abstract classes and will not be called. Consider using the 3D 
    version and only loop over the `dims`.
*/ 

virtual Matd2x2 getDMat2(){};
virtual Matd3x3 getDMat3(){};

};
#endif