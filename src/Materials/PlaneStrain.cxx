#include "Materials/PlaneStrain.h"
#include <iostream>

void PlaneStrain::ReadInput(H5IO &H5FIle){

    
    std::string dsetName;
    dsetName = "SimulationParameters/Emod";
    Emod = H5FIle.ReadScalar(dsetName);

    dsetName = "SimulationParameters/nu";
    nu = H5FIle.ReadScalar(dsetName);
}

void PlaneStrain::CalcDMatx(){

    DMatx.setZero();

    double param1 = Emod*(1-nu)/((1+nu)*(1-2*nu));
    DMatx(0,0) = param1;
    DMatx(1,1) = param1;
    DMatx(0,1) = param1*nu/(1.0-nu);
    DMatx(1,0) = param1*nu/(1.0-nu);
    DMatx(2,2) = (1.0-2.0*nu)*param1/(2.0*(1.0-nu));
}

Eigen::Matrix<double, 3,3> PlaneStrain::getDMatx(){

    return DMatx;
}
