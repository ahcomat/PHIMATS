#!/bin/bash

# Get the absolute path of where this script is sitting
PHIMATS_PATH=$(pwd)

H5LD_PATH=$(h5cc -show | awk '{for(i=1;i<=NF;i++) if($i ~ "-L") {print substr($i,3); exit}}')
H5ID_PATH=$(h5cc -show | awk '{for(i=1;i<=NF;i++) if($i ~ "-I") print substr($i,3)}')

# --- PETSc ---
while true; do
    echo
    echo "Enter the absolute path to your PETSc folder (e.g. /home/path/to/petsc):"
    read USER_PETSC_DIR
    
    # Check if petsc.h exists in the include directory
    if [ -f "$USER_PETSC_DIR/include/petsc.h" ]; then
        echo "✅ Valid PETSc directory found!"
        break
    else
        echo "❌ Error: Could not find 'include/petsc.h' in $USER_PETSC_DIR. Please try again."
    fi
done

# --- Eigen ---
while true; do
    echo 
    echo "Enter the absolute path to your Eigen folder (e.g. /home/path/to/eigen):"
    read USER_EIGEN_DIR
    
    # Eigen has a unique signature file to identify its location
    if [ -d "$USER_EIGEN_DIR" ] && [ -f "$USER_EIGEN_DIR/signature_of_eigen3_matrix_library" ]; then
        echo "✅ Valid Eigen directory found!"
        break
    else
        echo "❌ Error: This does not look like a valid Eigen folder. Please try again."
    fi
done

echo 
echo "Setting up PHIMATS environment at: $PHIMATS_PATH"

# Create (or overwrite) the environment file
cat <<EOF > phimats_env.sh
# --- PHIMATS Environment Variables ---
export PHIMATS_DIR=$PHIMATS_PATH
export PETSC_DIR=$USER_PETSC_DIR
export EIGEN=$USER_EIGEN_DIR
export PETSC_ARCH=arch-opt

export PYTHONPATH=\$PYTHONPATH:\$PHIMATS_DIR/src/FEM_utils
export PHIMATSINCLUDES=\$PHIMATS_DIR/src/include/
export PHIMATSLIBDIR=\$PHIMATS_DIR/src/lib/
export HDF5_USE_FILE_LOCKING=FALSE

# Auto-detect HDF5 
export H5LD=$H5LD_PATH
export H5ID=$H5ID_PATH
EOF

echo
echo "Environment file 'phimats_env.sh' has been generated."

# --- FINAL STEP: Connect to .bashrc ---

# Get the absolute path to the environment file we just created
ENV_FILE_PATH="$(pwd)/phimats_env.sh"
LINE_TO_ADD="source $ENV_FILE_PATH"
BASHRC_FILE="$HOME/.bashrc"

echo
echo "Checking if PHIMATS is connected to your terminal..."
echo 

# 1. Search if the line is already in .bashrc
if ! grep -qF "$LINE_TO_ADD" "$BASHRC_FILE"; then
    # 2. If not found, use >> to APPEND it safely to the end
    echo "" >> "$BASHRC_FILE" # Add a newline for neatness
    echo "$LINE_TO_ADD" >> "$BASHRC_FILE"
    echo "✅ Successfully added PHIMATS to $BASHRC_FILE"
else
    # 3. If found, do nothing so we don't create duplicates
    echo "ℹ️ PHIMATS is already connected to your terminal."
fi

echo "🏁 Setup complete! Please run 'source ~/.bashrc' or restart your terminal."
echo