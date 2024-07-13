#include <iostream>
#include <cmath>

#include "Materials/Trapping/TrapGB.h"

using namespace std;

TrapGB::TrapGB(string dimensions, H5IO &H5File, int iSet, string isoType)
    : BaseTrapping(isoType, dimensions) {

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

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/kappa_GB";
    kappa_GB = H5File.ReadScalar(dsetName);

    if (dims=="3D"){

    }
}

T_DMatx TrapGB::CalcDMatx(const double gPhi, const double T){

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
    mat2(0, 0) = DLx*pow(DTx/DLx, 4*gPhi);
    mat2(1, 1) = DLy*pow(DTy/DLy, 4*gPhi);

    DMatx = mat2;

    } else if (dims=="3D"){
    
    }

    return DMatx;
}

T_DMatx TrapGB::CalcTMatx(const double gPhi, const double T){

    // Lattice phase
    double DLx = D0x1*exp(-DQx1/(T*R));
    double DLy = D0y1*exp(-DQy1/(T*R));
    // Trapping phase
    double DTx = D0x2*exp(-DQx2/(T*R));
    double DTy = D0y2*exp(-DQy2/(T*R));

    // Variant for storing the trapping matrix D*kappa_GB/(RT)
    T_DMatx TMatx;

    if (dims=="2D"){

    Matd2x2 mat2 = Matd2x2::Zero();
    
    mat2.setZero();
    mat2(0, 0) = DLx*pow(DTx/DLx, 4*gPhi)*kappa_GB/(R*T);
    mat2(1, 1) = DLy*pow(DTy/DLy, 4*gPhi)*kappa_GB/(R*T);

    TMatx = mat2;

    } else if (dims=="3D"){
    
    }

    return TMatx;
}
