#include <iostream>

#include "Materials/Transport/HeatMassTransport.h"

using namespace std;

HeatMassTransport::HeatMassTransport(string dimensions, H5IO &H5File, int iSet, string isoType)
    : BaseTransport(isoType, dimensions) {

    string dsetName;

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/rho";
    rho = H5File.ReadScalar(dsetName);
    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/c";
    c = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/xKappa";
    double xKappa = H5File.ReadScalar(dsetName);
    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/yKappa";
    double yKappa = H5File.ReadScalar(dsetName);

    if (dims=="2D"){

    Matd2x2 mat2 = Matd2x2::Zero();
    
    mat2.setZero();
    mat2(0, 0) = xKappa;
    mat2(1, 1) = yKappa;

    KMatx = mat2;

    } else if (dims=="3D"){

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/zKappa";
    double zKappa = H5File.ReadScalar(dsetName);

    Matd3x3 mat3 = Matd3x3::Zero();

    mat3.setZero();
    mat3(0, 0) = xKappa;
    mat3(1, 1) = yKappa;
    mat3(2, 2) = zKappa;

    KMatx = mat3;

    }
}

T_DMatx HeatMassTransport::getKMatx() {

    return KMatx;
}
