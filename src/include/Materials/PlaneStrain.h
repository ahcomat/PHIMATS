#ifndef PLANESTRAIN_H
#define PLANESTRAIN_H

#include "Eigen/Dense"

// #include "NodeSets.h"
#include "H5IO.h"

class PlaneStrain{

public:

void ReadInput(H5IO &H5File);

void CalcDMatx();

Eigen::Matrix<double, 3,3> getDMatx();

private:

double nu;
double Emod;

Eigen::Matrix<double, 3,3> DMatx;

};
#endif