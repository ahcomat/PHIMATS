#include <iostream>

#include "Materials/Transport/HeatMassTransport.h"

using namespace std;

HeatMassTransport::HeatMassTransport(string dimensions, H5IO &H5File, int iSet, string isoType)
    : BaseTransport(isoType, dimensions) {

    string dsetName;

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/s";
    s = H5File.ReadScalar(dsetName);

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/Dx";
    double Dx = H5File.ReadScalar(dsetName);
    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/Dy";
    double Dy = H5File.ReadScalar(dsetName);

    if (dims=="2D"){

    Matd2x2 mat2 = Matd2x2::Zero();
    
    mat2.setZero();
    mat2(0, 0) = Dx;
    mat2(1, 1) = Dy;

    KMatx = mat2;

    } else if (dims=="3D"){

    dsetName = "Materials/Material_"+ std::to_string(iSet)+"/Dz";
    double Dz = H5File.ReadScalar(dsetName);

    Matd3x3 mat3 = Matd3x3::Zero();

    mat3.setZero();
    mat3(0, 0) = Dx;
    mat3(1, 1) = Dy;
    mat3(2, 2) = Dz;

    KMatx = mat3;

    }
}

T_DMatx HeatMassTransport::getKMatx() const {

    return KMatx;
}

double HeatMassTransport::getCapacity() const{

    return s;
}
