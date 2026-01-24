## **Installation Guide**

### **1️⃣ System Requirements**
- ✅ **Linux OS required** (Recommended: Ubuntu, Debian, Fedora)
- ✅ **Windows users:** Use **Windows Subsystem for Linux (WSL)**
- ⚠️ **macOS:** Potentially supported but untested 
- ❌ **No support for native Windows**  

---

### 2️⃣ Required Third-Party Software

✅ **Paraview** https://www.paraview.org  
✅ **Gmsh** https://gmsh.info/   

---

### 3️⃣ Download 

To download **PHIMATS**, clone the repository using **Git**

```sh
git clone https://github.com/ahcomat/PHIMATS.git
```

Set the project path by adding your local `PHIMATS` path in `~/.bashrc` file

```sh
export PHIMATS_DIR=/path/to/PHIMATS
```

---

### 4️⃣ Compiling 

#### 🐳 Option A: Docker (Recommended)

Uses a pre-built environment containing `PETSc`, `Eigen`, and `HDF5`, optimized for performance

1. Pull the Image:

```sh
docker pull ahcomat/phimats_dep:latest
```

2. Run this command from your PHIMATS root folder to start the container and mount your code:

```sh
docker run --rm -it \
    -v $PHIMATS_DIR:/home/phimats \
    ahcomat/phimats_dep:latest /bin/bash
```
3. Inside the container, in `PHIMATS` root folder:
```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../src && make
```

**Note** All outputs will be saved to your local drive.

#### 🛠️ Option B: Manual Installation (Advanced)

##### 1 Python Dependencies (via Conda)
 All required Python libraries are included in the **`phimats`** Conda environment.

To create the environment, run:
```sh
conda env create -f environment.yml
```

##### 2 Third-Party Libraries
✅ **C++20 compiler (GCC recommended)**  
✅ **CMake**  
✅ **Install `hdf5` library**  
✅ **Install `Eigen` library**  


##### 3 Installing PETSc
PETSc must be installed with an **optimized configuration** (`arch-opt`). Run the following command in your terminal:
```sh
./configure PETSC_ARCH=$PETSC_ARCH \
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
```

---

##### 4 Setting Up Environment Variables
Add the following lines **(after modifying the paths!)** to your `~/.bashrc` file:
```sh
# Provide your local /path/to below
export PETSC_DIR=/path/to/petsc
export EIGEN=/path/to/eigen
export PHIMATS_DIR=/path/to/PHIMATS

export PYTHONPATH=$PYTHONPATH:$PHIMATS_DIR/src/FEM_utils
export PHIMATSINCLUDES=$PHIMATS_DIR/src/include/
export PHIMATSLIBDIR=$PHIMATS_DIR/src/lib/
export PETSC_ARCH=arch-opt
export HDF5_USE_FILE_LOCKING=FALSE

# Detect HDF5 installation path automatically
export H5LD=$(h5cc -show | awk '{for(i=1;i<=NF;i++) if($i ~ "-L") print substr($i,3)}')
export H5ID=$(h5cc -show | awk '{for(i=1;i<=NF;i++) if($i ~ "-I") print substr($i,3)}')
export HDF5_USE_FILE_LOCKING=FALSE

# Detect MPICH installation path automatically
export MPICHID=$(mpicc -show | awk '{for(i=1;i<=NF;i++) if($i ~ "-I") print substr($i,3)}')
export MPICHLD=$(mpicc -show | awk '{for(i=1;i<=NF;i++) if($i ~ "-L") print substr($i,3)}')
```
After adding these, apply the changes by running:
```sh
source ~/.bashrc
```

##### 5 Compiling

In `PHIMATS` root folder:
```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../src && make
```
