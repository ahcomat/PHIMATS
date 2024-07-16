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

void H5IO::WriteScalar(string dsetName, double val){

    hid_t  file_id, dataset_id, dataspace_id;
    herr_t status;
    hsize_t dims[1];

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    dims[0] = 1;
    double Var[1];
    Var[0] = val;

    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(1, dims, NULL);

    dataset_id = H5Dcreate2(file_id, dataSetName, H5T_NATIVE_DOUBLE, dataspace_id, 
                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, Var);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);
}

void H5IO::ReadFieldFloat2D(string dsetName, const int row, const int col, vector<vector<double>>& Field){

    hid_t  file_id, dataset_id;
    herr_t status;

    double BufferField[row][col];

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

    vector<double> dummy1(col);

    for (int i=0; i<row; i++){
        for (int j=0; j<col; j++){
            dummy1.at(j) = BufferField[i][j];
        }
        Field.at(i) = dummy1;
    }
}

void H5IO::ReadFieldInt2D(string dsetName, const int row, const int col, vector<vector<int>>& Field){

    hid_t  file_id, dataset_id;
    herr_t status;

    int BufferField[row][col];

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    // Open existing file
    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
    // Open existing data set
    dataset_id = H5Dopen2(file_id, dataSetName, H5P_DEFAULT);
    // Read dataset buffer
    status = H5Dread(dataset_id, H5T_NATIVE_INT , H5S_ALL, H5S_ALL, H5P_DEFAULT, BufferField);
    // Close dataset
    status = H5Dclose(dataset_id);
    // Close file
    status = H5Fclose(file_id);

    vector<int> dummy1(col);

    for (int i=0; i<row; i++){
        for (int j=0; j<col; j++){
            dummy1.at(j) = BufferField[i][j];
        }
        Field.at(i) = dummy1;
    }
}

void H5IO::ReadFieldInt1D(string dsetName,  vector<int> &Field){

    hid_t  file_id, dataset_id;
    herr_t status;

    int Mx = Field.size();
    int BufferField[Mx];

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    // Open existing file
    file_id = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
    // Open existing data set
    dataset_id = H5Dopen2(file_id, dataSetName, H5P_DEFAULT);
    // Read dataset buffer
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, BufferField);
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

void H5IO::WriteArray_1D(std::string dsetName, const int xSize, const double *Array){

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

void H5IO::WriteStres(string dsetName, const int nNodes, const int nStres, const T_nodStres& Array){

    hid_t  file_id, dataset_id, dataspace_id;
    herr_t status;
    hsize_t dims[2];

    dims[0] = nNodes;
    dims[1] = nStres;

    const char* fileName = this->H5FileName.c_str();
    const char* dataSetName = dsetName.c_str();

    double BufferField [nNodes][nStres];

    if(nStres==2){

        for (int i=0; i<nNodes; i++){
            for (int j=0; j<nStres; j++){
                BufferField[i][j] = std::get<std::vector<ColVecd2>>(Array).at(i)(j);
            }
        }
    } else if(nStres==3){

        for (int i=0; i<nNodes; i++){
            for (int j=0; j<nStres; j++){
                BufferField[i][j] = std::get<std::vector<ColVecd3>>(Array).at(i)(j);
            }
        }
    } else if(nStres==6) {

        for (int i=0; i<nNodes; i++){
            for (int j=0; j<nStres; j++){
                BufferField[i][j] = std::get<std::vector<ColVecd6>>(Array).at(i)(j);
            }
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