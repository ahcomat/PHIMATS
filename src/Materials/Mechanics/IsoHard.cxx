#include <iostream>
#include <cmath>

#include "Materials/Mechanics/IsoHard.h"

IsoHard::IsoHard(string dimensions, H5IO& H5File, int iSet)
    : Elastic(dimensions, H5File, iSet) {

    try {
        // Read plasticity and hardening law
        Platicity = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/Plasticity");
        HardLaw = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/HardeningLaw");

        // Handle supported hardening laws
        if (HardLaw == "Linear") {

            K_hard = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/K_hard");

        } else if (HardLaw == "PowerLaw") {

            K_hard = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/K_hard");
            n_pow = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/n_pow");

        } else {

            throw invalid_argument("Undefined hardening law < " + HardLaw + " >");

        }
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}

double IsoHard::R_pow(double eps_eq){
    return K_hard*pow(eps_eq, n_pow);
}

double IsoHard::dR_pow(double eps_eq){

    if (eps_eq == 0){eps_eq = 1.0e-12;}

    return K_hard*n_pow*pow(eps_eq, n_pow-1);
}

double IsoHard::Mises3D(ColVecd6& sig3D){

    double term = 0.5 * (pow(sig3D(0) - sig3D(1), 2) + 
                         pow(sig3D(1) - sig3D(2), 2) + 
                         pow(sig3D(2) - sig3D(0), 2)) + 
                  3 * (pow(sig3D(3),2) + pow(sig3D(4),2) + pow(sig3D(5),2));

    return sqrt(term);
}

// double PEEQ3D(ColVecd6& eps3D){
//     double term = 2.0/3.0 * (eps3D)
// }



T_DMatx IsoHard::getDMatx() const{

    return DMatx_ep;
}
