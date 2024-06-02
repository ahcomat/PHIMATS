#include <iostream>

#include "Materials/Mechanics/PlaneStrain.h"

using namespace std;

PlaneStrain::PlaneStrain(H5IO &H5File, string isoType)
    : BaseMechanics(isoType, "2D"){
    
    string dsetName; 

    if(isotropy=="Isotropic"){

        /**
         * Reads Young's modulus and Poisson's ratio
         */
        dsetName = "SimulationParameters/Emod";
        double Emod = H5File.ReadScalar(dsetName);
        dsetName = "SimulationParameters/nu";
        double nu = H5File.ReadScalar(dsetName);

        DMatx.setZero();

        double param1 = Emod*(1-nu)/((1+nu)*(1-2*nu));
        DMatx(0,0) = param1;
        DMatx(1,1) = param1;
        DMatx(0,1) = param1*nu/(1.0-nu);
        DMatx(1,0) = param1*nu/(1.0-nu);
        DMatx(2,2) = (1.0-2.0*nu)*param1/(2.0*(1.0-nu));

    } else {
        cout << "ERROR! Undefined material isotropy < " << isotropy << " >\n";
        cout << "Terminating! \n";
        exit(10);
    }
}

T_DMatx PlaneStrain::getDMatx(){

    return DMatx;
}
