#include <iostream>
#include <cmath>

#include "Materials/Trapping/TrapHLGB.h"

using namespace std;

TrapHLGB::TrapHLGB(string dimensions, H5IO& H5File, int iSet, Logger& logger)
    : BaseTrapping(dimensions, logger) {

    string dsetName;

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/D0x1";
    D0x1 = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/D0y1";
    D0y1 = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/DQx1";
    DQx1 = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/DQy1";
    DQy1 = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zeta_HAGB";
    zeta_HAGB = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zeta_LAGB";
    zeta_LAGB = H5File.ReadScalar(dsetName);

    if (dims=="3D"){

    }
}

T_DMatx TrapHLGB::CalcDMatx(const double phi, const double T){

    // Lattice phase
    double DLx = D0x1*exp(-DQx1/(T*R));
    double DLy = D0y1*exp(-DQy1/(T*R));
    
    // Variant for storing the diffusivity (conductivity) matrix
    T_DMatx DMatx;

    if (dims=="2D"){

    Matd2x2 mat2 = Matd2x2::Zero();
    
    mat2.setZero();
    mat2(0, 0) = DLx;
    mat2(1, 1) = DLy;

    DMatx = mat2;

    } else if (dims=="3D"){
    
    }

    return DMatx;
}
