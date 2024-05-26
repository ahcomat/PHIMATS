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
 * Updates (when, what and who)
 * 
 */

#include <iostream>
#include <exception>
#include "H5IO.h"

#ifndef DEBUG
#define at(x) operator[](x)
#endif

H5IO::H5IO(std::string H5FName)
    : H5FileName(H5FName) {

}

H5IO::~H5IO(){

    // Exit message
    cout << "H5IO exited correctly" << "\n";
}

double H5IO::ReadScalar(string dsetName){

    hid_t  file_id, dataset_id;
    herr_t status;

    // Buffer for getting data.
    double Var[1];

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();
    
    // Open existing file
    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
    // Open existing data set
    dataset_id = H5Dopen2(file_id, dataSetName, H5P_DEFAULT);

    // Read dataset buffer
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Var);
    // Close dataset
    status = H5Dclose(dataset_id);
    // Close file
    status = H5Fclose(file_id);

    return Var[0];
}

void H5IO::ReadFieldInt1D(string dsetName,  vector<int> &Field){

    hid_t  file_id, dataset_id;
    herr_t status;

    int Mx = Field.size();
    double BufferField[Mx];

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    // Open existing file
    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
    // Open existing data set
    dataset_id = H5Dopen2(file_id, dataSetName, H5P_DEFAULT);
    // Read dataset buffer
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, BufferField);
    // Close dataset
    status = H5Dclose(dataset_id);
    // Close file
    status = H5Fclose(file_id);

    for (int i=0; i<Mx; i++){
        Field.at(i) = BufferField[i];
    }
}

void H5IO::ReadFieldDoub1D(string dsetName, vector<double>& Field){

    hid_t  file_id, dataset_id;
    herr_t status;

    int Mx = Field.size();
    double BufferField[Mx];

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    // Open existing file
    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
    // Open existing data set
    dataset_id = H5Dopen2(file_id, dataSetName, H5P_DEFAULT);
    // Read dataset buffer
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, BufferField);
    // Close dataset
    status = H5Dclose(dataset_id);
    // Close file
    status = H5Fclose(file_id);

    for (int i=0; i<Mx; i++){
        Field.at(i) = BufferField[i];
    }
}


void H5IO::WriteArray_1D(std::string dsetName, int xSize, const double *Array){

    hid_t  file_id, dataset_id, dataspace_id;
    herr_t status;
    hsize_t dims[1];

    dims[0] = xSize;

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(1, dims, NULL);

    dataset_id = H5Dcreate2(file_id, dataSetName, H5T_NATIVE_DOUBLE, dataspace_id, 
                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, Array);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);
}

void H5IO::WriteStres3(string dsetName, int nNodes, int nStres, vector<ColVecd3> &Array){

    hid_t  file_id, dataset_id, dataspace_id;
    herr_t status;
    hsize_t dims[2];

    dims[0] = nNodes;
    dims[1] = nStres;

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    double BufferField [nNodes][nStres];

    for (int i=0; i<nNodes; i++){
        for (int j=0; j<nStres; j++){
            BufferField[i][j] = Array.at(i)(j);
        }
    }

    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(2, dims, NULL);

    dataset_id = H5Dcreate2(file_id, dataSetName, H5T_NATIVE_DOUBLE, dataspace_id, 
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, BufferField);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);
}

// void H5IO::WriteArray_3D(std::string dsetName, int xSize, int ySize, int zSize, const double *Array){

//     hid_t  file_id, dataset_id, dataspace_id;
//     herr_t status;
//     hsize_t dims[3];

//     dims[0] = xSize;
//     dims[1] = ySize;
//     dims[2] = zSize;

//     float *outBuffer = new float[xSize*ySize*zSize];

//     for (int i=0; i<xSize*ySize*zSize; i++){

//         outBuffer[i] = (float)Array[i];
//     }

//     const char* fileName = this->H5FileName.c_str();
//     const char* dataSetName = dsetName.c_str();

//     file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

//     dataspace_id = H5Screate_simple(3, dims, NULL);

//     dataset_id = H5Dcreate2(file_id, dataSetName, H5T_NATIVE_FLOAT, dataspace_id, 
//                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
             
//     status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
//              H5P_DEFAULT, outBuffer);

//     status = H5Dclose(dataset_id);
//     status = H5Sclose(dataspace_id);
//     status = H5Fclose(file_id);
// }
