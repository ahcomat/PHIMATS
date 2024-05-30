/**
 * @file H5IO.cxx
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Wrapper class for HDF5 files.
 * @date 2024-05-23
 * 
 * @copyright Copyright (c) 2024
 * 
 * @todo Update for C++ API
 *  
 */

#ifndef H5IO_H
#define H5IO_H

#include <string>
#include <vector>
#include "Matrix.h"
#include "hdf5.h"

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