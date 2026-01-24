FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Basic Tools + GUI support libraries
RUN apt-get update && apt-get install -y \
    build-essential cmake git wget python3 python3-pip \
    libhdf5-dev gfortran \
    && rm -rf /var/lib/apt/lists/*

# Python (No Conda)
RUN pip3 install --no-cache-dir numpy scipy h5py \
    matplotlib meshio ipykernel nbformat

# Eigen 
WORKDIR /opt
RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz && \
    tar -xf eigen-3.4.0.tar.gz && mv eigen-3.4.0 eigen

# PETSs
WORKDIR /opt/petsc
RUN git clone -b release https://gitlab.com/petsc/petsc.git .
ENV PETSC_DIR=/opt/petsc
ENV PETSC_ARCH=arch-opt
RUN ./configure PETSC_ARCH=$PETSC_ARCH \
    --with-shared-libraries=1 \
    --with-fortran-bindings=0 \
    --with-debugging=0 \
    --download-f2cblaslapack \
    --download-mpich \
    --download-scalapack \
    --download-mumps \
    --download-cmake \
    COPTFLAGS='-O3 -march=native -mtune=native' \
    CXXOPTFLAGS='-O3 -march=native -mtune=native' \
    FOPTFLAGS='-O3 -march=native -mtune=native' && \
    make all -j$(nproc)

# Environments
ENV EIGEN=/opt/eigen
ENV PHIMATS_DIR=/home/phimats
ENV PHIMATSINCLUDES=$PHIMATS_DIR/src/include/
ENV PHIMATSLIBDIR=$PHIMATS_DIR/src/lib/
ENV HDF5_USE_FILE_LOCKING=FALSE
ENV H5LD=/usr/lib/x86_64-linux-gnu/hdf5/serial 
ENV H5ID=/usr/include/hdf5/serial/
ENV MPICHID=/usr/include/x86_64-linux-gnu/mpich/mpi.h
ENV MPICHLD=/usr/include/x86_64-linux-gnu/mpich/libexport 

WORKDIR $PHIMATS_DIR
ENV PYTHONPATH="/home/phimats/src/FEM_utils"