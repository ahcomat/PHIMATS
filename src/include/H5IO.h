#ifndef H5IO_H
#define H5IO_H

#include <string>
#include <vector>
#include "Eigen/Dense"
#include "hdf5.h"

typedef Eigen::Vector<double, 3> ColVecd3;       // nStres 
typedef Eigen::Vector<double, 6> ColVecd6;       // nStres 

using namespace std;


class H5IO{

public:

H5IO(std::string H5FileName);
~H5IO();

/**
 * @brief Reads a scalar input from HDF5 file. 
 * 
 * @param dsetName  
 */
double ReadScalar(std::string dsetName);

/**
 * @brief Reads a 1D field input from HDF5 file. 
 * 
 * @param dsetName 
 * @return double/int 
 */
// void ReadFieldInt1D(std::string dsetName, Array1D<int>& Field);

void ReadFieldInt1D(string dsetName, vector<int>& Field);

void ReadFieldDoub1D(string dsetName, vector<double>& Field);

void WriteArray_1D(string dsetName, int xSize, const double *Array);

void WriteStres3(string dsetName, int xSize, int ySize, vector<ColVecd3> &Array);
void WriteStres6(string dsetName, int xSize, int ySize, vector<ColVecd6> &Array);


// void WriteArray_3D(std::string dsetName, int xSize, int ySize, int zSize, const double *Array);

private:

const string H5FileName;         /// Name of HDF file. 


};

#endif