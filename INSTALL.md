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

---

### 4️⃣ Compiling 

#### 🐳 Option A: Docker (Recommended)

Uses a pre-built environment containing `PETSc`, `Eigen`, and `HDF5`, optimized for performance

1. Pull the Image:

```sh
docker pull ahcomat/phimats_dep:latest
```

2. Run this command from your `PHIMATS` root folder to start the container and mount your code:

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

**Note:** All build and simulation outputs are saved directly to your local WSL drive within the project folder.

- **Persistence**: Because we "mount" your local `PHIMATS` directory into the container, these files remain on your disk even after you close the `docker` terminal or restart your computer.
- **Accessibility**: You can access these files directly through Windows File Explorer by typing `explorer.exe .` in your WSL terminal.

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

Clone `PETSc`
```
git clone -b release https://gitlab.com/petsc/petsc.git 
cd petsc
```

`PETSc` must be installed with an **optimized configuration** (`arch-opt`). Run the following command in your terminal:
```sh
./configure PETSC_ARCH=arch-opt \
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
    make PETSC_DIR=$(pwd) PETSC_ARCH=arch-opt all -j$(nproc)
```

---

##### 4 Setting Up Environment Variables

Run the `configure_env.sh` script to set up the environment variables

* Make executable

```sh
chmod +x configure_env.sh
```
* Run the script and follow the prompts

```
./configure_env.sh
```
* Apply the changes by running
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
