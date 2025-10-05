#include <iostream>

#include "Materials/Trapping/MechTrap.h"

using namespace std;

MechTrap::MechTrap(string dimensions, H5IO& H5File, int iSet, Logger& logger)
    : BaseTrapping(dimensions, logger){

    string dsetName;

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/D0x";
    D0x = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/D0y";
    D0y = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/DQx";
    DQx = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/DQy";
    DQy = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/Vh";
    Vh = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/m";
    m = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/s";
    s = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zeta_rho";
    zeta_rho = H5File.ReadScalar(dsetName);

    dsetName = "SimulationParameters/Trapping";
    string Trapping = H5File.ReadString(dsetName);

    if (Trapping=="MechTrappingPFF"){

        dsetName = "Materials/Material_"+ std::to_string(iSet)+"/Zd";
        Zd = H5File.ReadScalar(dsetName);

    } 


    if (dims=="3D"){

        dsetName = "Materials/Material_"+ std::to_string(iSet)+"/D0z";
        D0z = H5File.ReadScalar(dsetName);

        dsetName = "Materials/Material_"+ std::to_string(iSet)+"/DQz";
        DQz = H5File.ReadScalar(dsetName);

    }
}

T_DMatx MechTrap::CalcDMatx(const double phi, const double T){

    double DLx = D0x*exp(-DQx/(T*R));
    double DLy = D0y*exp(-DQy/(T*R));

    // Variant for storing the diffusivity matrix
    T_DMatx DMatx;

    if (dims=="2D"){

        Matd2x2 mat2 = Matd2x2::Zero();
    
        mat2.setZero();
        mat2(0, 0) = DLx*(1 + m*phi);
        mat2(1, 1) = DLy*(1 + m*phi);

        DMatx = mat2;

    } else if (dims=="3D"){

        double DLz = D0z*exp(-DQz/(T*R));

        Matd3x3 mat3 = Matd3x3::Zero();
    
        mat3.setZero();
        mat3(0, 0) = DLx*(1 + m*phi);
        mat3(1, 1) = DLy*(1 + m*phi);
        mat3(2, 2) = DLz*(1 + m*phi);

        DMatx = mat3;

    }

    return DMatx;
}
