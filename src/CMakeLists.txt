cmake_minimum_required(VERSION 3.16)
project(PHIMATS VERSION 1.3.0 LANGUAGES CXX)

# Force C++20 and Position Independent Code (-fPIC) for the 90MB library
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Map Environment Variables
set(EIGEN_DIR $ENV{EIGEN})
set(HDF5_INC  $ENV{H5ID})
set(HDF5_LIB  $ENV{H5LD})
set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})

# Include Directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN_DIR}
    ${HDF5_INC}
    ${PETSC_DIR}/include
    ${PETSC_DIR}/${PETSC_ARCH}/include
)

# Collect all Source Files
file(GLOB_RECURSE PHIMATS_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

# Create the Shared Library (libphimats.so)
add_library(phimats SHARED ${PHIMATS_SOURCES})

# Link Dependencies
# Use the H5LD path you exported to find the serial HDF5 libraries
target_link_libraries(phimats 
    PRIVATE 
    ${HDF5_LIB}/libhdf5.so
    ${HDF5_LIB}/libhdf5_cpp.so
    ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so
)

# Set output directory
set_target_properties(phimats PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)