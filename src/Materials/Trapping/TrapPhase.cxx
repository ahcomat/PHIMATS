#include <iostream>
#include <cmath>

#include "Materials/Trapping/TrapPhase.h"

using namespace std;

TrapPhase::TrapPhase(string dimensions, H5IO &H5File, int iSet, Logger& logger)
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

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/D0x2";
    D0x2 = H5File.ReadScalar(dsetName);
    
    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/D0y2";
    D0y2 = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/DQy2";
    DQx2 = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/DQy2";
    DQy2 = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zeta_j";
    zeta_j = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zeta_jj";
    zeta_jj = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zeta_ij";
    zeta_ij = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zeta_ii";
    zeta_ii = H5File.ReadScalar(dsetName);

    if (dims=="3D"){

    }
}

T_DMatx TrapPhase::CalcDMatx(const double phi, const double T){

    // Lattice phase
    double DLx = D0x1*exp(-DQx1/(T*R));
    double DLy = D0y1*exp(-DQy1/(T*R));
    // GBs
    double DTx = D0x2*exp(-DQx2/(T*R));
    double DTy = D0y2*exp(-DQy2/(T*R));

    // Variant for storing the diffusivity (conductivity) matrix
    T_DMatx DMatx;

    if (dims=="2D"){

    Matd2x2 mat2 = Matd2x2::Zero();
    
    mat2.setZero();
    mat2(0, 0) = DLx*(1-phi) + phi*DTx;
    mat2(1, 1) = DLy*(1-phi) + phi*DTy;

    DMatx = mat2;

    } else if (dims=="3D"){
    
    }

    return DMatx;
}
