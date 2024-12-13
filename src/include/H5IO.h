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

/**
 * @brief Constructor. Checks if the file exists. 
 * 
 * @param H5FileName 
 */
H5IO(string H5FileName);
~H5IO();

/**
 * @brief Reads a scalar input from HDF5 file. 
 * 
 * @param dsetName  
 */
double ReadScalar(const string& dsetName);

/**
 * @brief Write scalar `val`.
 * 
 * @param dsetName 
 * @param val 
 */
void WriteScalar(const string& dsetName, double val);

/**
 * @brief Reads a string variable. 
 * 
 * @param dsetName 
 * @return string 
 */
string ReadString(const string& dsetName);

void ReadFieldFloat2D(const string &dsetName, const int row, const int col, vector<vector<double>>& Field);

void ReadFieldInt2D(const string& dsetName, const int row, const int col, vector<vector<int>>& Field);

void ReadFieldInt1D(const string& dsetName, vector<int>& Field);

void ReadFieldDoub1D(const string& dsetName, vector<double>& Field);

void WriteArray_1D(const string& dsetName, const int xSize, const double *Array);

void WriteStres(const string& dsetName, const int xSize, const int ySize, const T_nodStres& Array);

private:

const string H5FileName;         /// @brief Name of HDF file. 

};

#endif