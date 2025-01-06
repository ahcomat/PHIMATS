/**
 * @file H5IO.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Wrapper class for HDF5 files.
 * @date 2024-05-23
 * 
 * @todo Update for C++ API
 *  
 * @copyright Copyright (C) 2024 Abdelrahman Hussein
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  
 */

#ifndef H5IO_H
#define H5IO_H

#include <iostream>
#include <exception>
#include <string>
#include <vector>
#include "Matrix.h"
#include "hdf5.h"


#ifndef DEBUG
#define at(x) operator[](x)
#endif

using namespace std;

class H5IO{

public:

/**
 * @brief Constructor. Checks if the file exists. 
 * 
 * @param H5FileName Name of hdf5 file. 
 */
H5IO(string H5FileName);
~H5IO();

/**
 * @brief Reads a scalar input from hdf5 file. 
 * 
 * @param dsetName Name of dataset.
 */
double ReadScalar(const string& dsetName);

/**
 * @brief Write scalar.
 * 
 * @param dsetName Name of dataset.
 * @param val Value.
 */
void WriteScalar(const string& dsetName, double val);

/**
 * @brief Reads a string variable. 
 * 
 * @param dsetName 
 * @return string 
 */
string ReadString(const string& dsetName);

/**
 * @brief Reads a 1D field data.
 * 
 * @tparam T Template type (`int` or `double`).
 * @param dsetName Name of dataset.
 * @param Field Initialized vector.
 */
template <typename T>
void ReadField1D(const std::string& dsetName, std::vector<T>& Field) {
    try {

        int Mx = Field.size();

        // Allocate buffer for the dataset
        std::vector<T> buffer(Mx);

        // Open HDF5 file
        hid_t file_id = H5Fopen(this->H5FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            throw std::runtime_error("Failed to open HDF5 file: " + this->H5FileName);
        }

        // Open dataset
        hid_t dataset_id = H5Dopen2(file_id, dsetName.c_str(), H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: " + dsetName);
        }

        // Determine HDF5 data type for T
        hid_t h5_type;
        if constexpr (std::is_same<T, double>::value) {
            h5_type = H5T_NATIVE_DOUBLE;
        } else if constexpr (std::is_same<T, int>::value) {
            h5_type = H5T_NATIVE_INT;
        } else {
            H5Dclose(dataset_id);
            H5Fclose(file_id);
            throw std::invalid_argument("Unsupported data type for ReadField1D.");
        }

        // Read the dataset into the buffer
        herr_t status = H5Dread(dataset_id, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
        if (status < 0) {
            H5Dclose(dataset_id);
            H5Fclose(file_id);
            throw std::runtime_error("Failed to read dataset: " + dsetName);
        }

        // Copy buffer to Field
        for (int i = 0; i < Mx; i++) {
            Field.at(i) = buffer[i];
        }

        // Close dataset and file
        H5Dclose(dataset_id);
        H5Fclose(file_id);

    } catch (const std::exception& e) {
        std::cerr << "ERROR in ReadField1D: " << e.what() << std::endl;
        throw; // Rethrow for higher-level handling
    }
}

/**
 * @brief Reads a 2D field data. 
 * 
 * @tparam T Template type (`int` or `double`)
 * @param dsetName Name of dataset. 
 * @param dim1 Number of rows.
 * @param dim2 Number of columns. 
 * @param Field Initialized vector. 
 */
template <typename T>
void ReadField2D(const string& dsetName, const int dim1, const int dim2, vector<vector<T>>& Field) {

    try {
        // Open HDF5 file
        hid_t file_id = H5Fopen(H5FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            throw std::runtime_error("Failed to open HDF5 file: " + H5FileName);
        }

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id, dsetName.c_str(), H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Fclose(file_id);
            throw std::runtime_error("Failed to open dataset: " + dsetName);
        }

        // Determine HDF5 data type for T
        hid_t h5_type;
        if constexpr (std::is_same<T, double>::value) {
            h5_type = H5T_NATIVE_DOUBLE;
        } else if constexpr (std::is_same<T, int>::value) {
            h5_type = H5T_NATIVE_INT;
        } else {
            H5Dclose(dataset_id);
            H5Fclose(file_id);
            throw std::invalid_argument("Unsupported data type for ReadField2D.");
        }

        // Allocate buffer for reading the dataset
        std::vector<T> buffer(dim1 * dim2);

        // Read the dataset into the 1D buffer
        herr_t status = H5Dread(dataset_id, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
        if (status < 0) {
            H5Dclose(dataset_id);
            H5Fclose(file_id);
            throw std::runtime_error("Failed to read dataset: " + dsetName);
        }

        vector<T> dummy1(dim2);
        // Resize and populate the 2D vector
        // Field.resize(dim1, std::vector<T>(dim2));
        for (int i = 0; i < dim1; ++i) {
            for (int j = 0; j < dim2; ++j) {
                dummy1.at(j) = buffer.at(i * dim2 + j);
            }
            Field.at(i) = dummy1;
        }

        // Close the dataset and file
        H5Dclose(dataset_id);
        H5Fclose(file_id);

    } catch (const std::exception& e) {
        std::cerr << "ERROR in ReadField2D: " << e.what() << std::endl;
        throw; // Rethrow for higher-level handling if needed
    }
}

/**
 * @brief Writes array.
 * 
 * @param dsetName Name of dataset. 
 * @param xSize Size of array.
 * @param Array The array to write.
 */
void WriteArray1D(const string& dsetName, const int xSize, const double *Array);

/**
 * @brief Writes a tensor field in Voigt notation.
 * 
 * @param dsetName Name of data set.
 * @param nNodes Total number of nodes. 
 * @param nStres Sizer of tensor in Voigt notation (`2`, `3` or `6`)
 * @param Array The tensor field container.
 */
void WriteTensor(const string& dsetName,const int nNodes, const int nStres, const T_nodStres& Array);

private:

/// @brief Name of HDF file. 
const string H5FileName;         

};

#endif