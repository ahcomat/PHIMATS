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

void ReadFieldFloat2D(string dsetName, const int row, const int col, vector<vector<double>>& Field);

void ReadFieldInt2D(string dsetName, const int row, const int col, vector<vector<int>>& Field);

void ReadFieldInt1D(string dsetName, vector<int>& Field);

void ReadFieldDoub1D(string dsetName, vector<double>& Field);

void WriteArray_1D(string dsetName, int xSize, const double *Array);

void WriteStres(string dsetName, int xSize, int ySize, const T_nodStres& Array);

private:

const string H5FileName;         /// @brief Name of HDF file. 

};

#endif