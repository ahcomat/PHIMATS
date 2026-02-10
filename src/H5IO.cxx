#include "H5IO.h"

#ifndef DEBUG
#define at(x) operator[](x)
#endif

H5IO::H5IO(string H5FName, Logger& logger)
    : H5FileName(H5FName), logger(logger) {

    // Attempt to open the file
    hid_t file_id = H5Fopen(H5FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0) {
        throw std::runtime_error("Error: Unable to open file " + H5FileName);
    } 
    // Close the file after the check
    H5Fclose(file_id);

}

H5IO::~H5IO(){}

double H5IO::ReadScalar(const string& dsetName){

    hid_t  file_id = -1, dataset_id = -1;

    // Buffer for data.
    double value = 0.0;

    try {
        // Open the file
        file_id = H5Fopen(H5FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) throw std::runtime_error("Failed to open file "+ H5FileName);

        // Open the dataset
        dataset_id = H5Dopen2(file_id, dsetName.c_str(), H5P_DEFAULT);
        if (dataset_id < 0) throw std::runtime_error("Failed to open dataset " + dsetName);

        // Read the data
        if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) < 0)
            throw std::runtime_error("Failed to read dataset " + dsetName);

    } catch (const std::runtime_error& e) {
        logger.log("Exception caught in H5IO::ReadScalar:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("    Critical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
    }

    // Ensure resources are released
    if (dataset_id >= 0) H5Dclose(dataset_id);
    if (file_id >= 0) H5Fclose(file_id);

    return value;
}

void H5IO::WriteScalar(const std::string& dsetName, double val) {

    hid_t file_id = -1, dataset_id = -1, dataspace_id = -1;

    try {
        // Open the HDF5 file
        file_id = H5Fopen(H5FileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (file_id < 0) throw std::runtime_error("Failed to open file " + H5FileName);

        // Create a simple dataspace
        hsize_t dims[1] = {1};
        dataspace_id = H5Screate_simple(1, dims, NULL);
        if (dataspace_id < 0) throw std::runtime_error("Failed to create dataspace");

        // Create the dataset
        dataset_id = H5Dcreate2(file_id, dsetName.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id < 0) throw std::runtime_error("Failed to create dataset " + dsetName);

        // Write data 
        if (H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val) < 0)
            throw std::runtime_error("Failed to write data " + dsetName);

    } catch (const std::runtime_error& e) {
        logger.log("Exception caught in H5IO::WriteScalar:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
    }

    if (dataset_id >= 0) H5Dclose(dataset_id);
    if (dataspace_id >= 0) H5Sclose(dataspace_id);
    if (file_id >= 0) H5Fclose(file_id);
}

string H5IO::ReadString(const string& dsetName) {

    hid_t file_id = -1, dataset_id = -1, datatype = -1, dataspace = -1;
    herr_t status;

    char* buffer = nullptr;
    string result;

    try {
        file_id = H5Fopen(H5FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) throw std::runtime_error("Failed to open HDF5 file " + H5FileName);

        dataset_id = H5Dopen2(file_id, dsetName.c_str(), H5P_DEFAULT);
        if (dataset_id < 0) throw std::runtime_error("Failed to open dataset " + dsetName);

        // Get datatype and dataspace
        datatype = H5Dget_type(dataset_id);
        dataspace = H5Dget_space(dataset_id);

        // Ensure the datatype is a string
        if (H5Tget_class(datatype) != H5T_STRING) {
            throw std::runtime_error("Dataset is not a string type!");
        }

        // Get string size and allocate memory
        size_t str_size = H5Tget_size(datatype);
        buffer = new char[str_size + 1];

        // Read the dataset
        status = H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
        if (status < 0) throw std::runtime_error("Failed to read dataset " + dsetName);

        buffer[str_size] = '\0';  // Null-terminate the string
        result = string(buffer);

    } catch (const std::exception& e) {
        logger.log("Exception caught in H5IO::ReadString:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
    }

    // Cleanup resources
    if (buffer) delete[] buffer;
    if (datatype >= 0) H5Tclose(datatype);
    if (dataspace >= 0) H5Sclose(dataspace);
    if (dataset_id >= 0) H5Dclose(dataset_id);
    if (file_id >= 0) H5Fclose(file_id);
    return result;
}

void H5IO::WriteArray1D(const string& dsetName, const int xSize, const double *Array, int dscale){

    try{

        hid_t  file_id, dataset_id, dataspace_id, dcpl;
        hsize_t dims[1]{ (hsize_t)xSize };
        hsize_t chunk_dims[1] = { (hsize_t)xSize }; // Chunk size equals total size for 1D

        file_id = H5Fopen(H5FileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (file_id < 0) throw std::runtime_error("Failed to open HDF5 file: " + H5FileName);

        dataspace_id = H5Screate_simple(1, dims, NULL);

        // --- Compression Setup ---
        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(dcpl, H5D_CHUNKED);
        H5Pset_chunk(dcpl, 1, chunk_dims);
        H5Pset_scaleoffset(dcpl, H5Z_SO_FLOAT_DSCALE, dscale);
        H5Pset_shuffle(dcpl);
        H5Pset_deflate(dcpl, 7); // Level 6 is a good balance of speed/compression

        dataset_id = H5Dcreate2(file_id, dsetName.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
                                H5P_DEFAULT, dcpl, H5P_DEFAULT);
        
        if (dataset_id < 0) {
            H5Pclose(dcpl);
            H5Sclose(dataspace_id);
            H5Fclose(file_id);
            throw std::runtime_error("Failed to create dataset: " + dsetName);
        }

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Array);

        H5Dclose(dataset_id);
        H5Pclose(dcpl);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);

    } catch (const std::exception& e) {
        logger.log("Exception caught in H5IO::WriteArray1D:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
    }
}

void H5IO::WriteTensor(const string& dsetName, const int nNodes, const int nStres, const T_nodStres& Array, int dscale){

    try {

        hid_t file_id, dataset_id, dataspace_id, dcpl;
        hsize_t dims[2] = { (hsize_t)nNodes, (hsize_t)nStres };
        // Chunking by rows (e.g., 1000 nodes at a time) or the whole field
        hsize_t chunk_dims[2] = { (hsize_t)nNodes, (hsize_t)nStres };

        std::vector<double> BufferField(nNodes * nStres);

        if(nStres==2){

            for (int i=0; i<nNodes; i++){
                for (int j=0; j<nStres; j++){
                    BufferField[i*nStres + j] = std::get<std::vector<ColVecd2>>(Array).at(i)(j);
                }
            }
        } else if(nStres==3){

            for (int i=0; i<nNodes; i++){
                for (int j=0; j<nStres; j++){
                    BufferField[i*nStres + j] = std::get<std::vector<ColVecd3>>(Array).at(i)(j);
                }
            }

        } else if(nStres==4){

            for (int i=0; i<nNodes; i++){
                for (int j=0; j<nStres; j++){
                    BufferField[i*nStres + j] = std::get<std::vector<ColVecd4>>(Array).at(i)(j);
                }
            }

        } else if(nStres==6) {

            for (int i=0; i<nNodes; i++){
                for (int j=0; j<nStres; j++){
                    BufferField[i*nStres + j] = std::get<std::vector<ColVecd6>>(Array).at(i)(j);
                }
            }        
        }

        // Open the HDF5 file
        file_id = H5Fopen(H5FileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        dataspace_id = H5Screate_simple(2, dims, NULL);

        // --- Compression Setup ---
        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(dcpl, H5D_CHUNKED);   
        H5Pset_chunk(dcpl, 2, chunk_dims);
        // Scale-Offset (Lossy compression: keeps 7 decimal places)
        // This reduces entropy by discarding unnecessary precision
        H5Pset_scaleoffset(dcpl, H5Z_SO_FLOAT_DSCALE, dscale);
        // Shuffle (Reorders bytes for better pattern matching)
        H5Pset_shuffle(dcpl);
        H5Pset_deflate(dcpl, 7); 
        // --------------------------

        dataset_id = H5Dcreate2(file_id, dsetName.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
                                H5P_DEFAULT, dcpl, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, BufferField.data());

        H5Dclose(dataset_id);
        H5Pclose(dcpl);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);

    } catch (const std::runtime_error& e) {
        logger.log("\nException caught in H5IO::WriteTensor:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
    }

}