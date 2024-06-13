#include <iostream>

#include "Materials/Mechanics/Elastic3D.h"

Elastic3D::Elastic3D(H5IO &H5File, int iSet, string isoType)
    : BaseMechanics(isoType, "3D") {

    string dsetName; 

    if(isotropy=="Isotropic"){

        dsetName = "Materials/Material_"+ std::to_string(iSet)+"/Emod";
        double Emod = H5File.ReadScalar(dsetName);
        dsetName = "Materials/Material_"+ std::to_string(iSet)+"/nu";
        double nu = H5File.ReadScalar(dsetName);

        DMatx.setZero();

        double ho = Emod*nu/((1+nu)*(1-2*nu));
        double uo = Emod/(2*(1+nu));

        DMatx(0,0) = ho+2*uo;
        DMatx(1,1) = ho+2*uo;
        DMatx(2,2) = ho+2*uo;

        DMatx(0,1) = ho;
        DMatx(0,2) = ho;
        DMatx(1,0) = ho;
        DMatx(1,2) = ho;
        DMatx(2,0) = ho;
        DMatx(2,1) = ho;
        
        DMatx(3,3) = 2*uo;
        DMatx(4,4) = 2*uo;
        DMatx(5,5) = 2*uo;

    } else {
        cout << "ERROR! Undefined material isotropy < " << isotropy << " >\n";
        cout << "Terminating! \n";
        exit(10);
    }
}

T_DMatx Elastic3D::getDMatx(){

    return DMatx;
}
