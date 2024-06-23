#include <iostream>

#include "Materials/Trapping/MechTrap.h"

using namespace std;

MechTrap::MechTrap(string dimensions, H5IO &H5File, int iSet, string isoType)
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

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/Vh";
    Vh = H5File.ReadScalar(dsetName);

    if (dims=="3D"){

    }
}

T_DMatx MechTrap::CalcKMatx(const double T){

    double DLx = D0x1*exp(-DQx1/(T*R));
    double DLy = D0y1*exp(-DQy1/(T*R));

    // Variant for storing the diffusivity (conductivity) matrix
    T_DMatx KMatx;

    if (dims=="2D"){

    Matd2x2 mat2 = Matd2x2::Zero();
    
    mat2.setZero();
    mat2(0, 0) = DLx;
    mat2(1, 1) = DLy;

    KMatx = mat2;

    } else if (dims=="3D"){
    
    }

    return KMatx;
}

T_DMatx MechTrap::CalcTMatx(const double T){

    double DLx = D0x1*exp(-DQx1/(T*R));
    double DLy = D0y1*exp(-DQy1/(T*R));

    // Variant for storing the trapping matrix D*zeta/(RT)
    T_DMatx TMatx;

    if (dims=="2D"){

    Matd2x2 mat2 = Matd2x2::Zero();
    
    mat2.setZero();
    mat2(0, 0) = (DLx*Vh)/(R*T);
    mat2(1, 1) = (DLy*Vh)/(R*T);

    TMatx = mat2;

    } else if (dims=="3D"){
    
    }

    return TMatx;
}
