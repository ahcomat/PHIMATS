#include <iostream>
#include <stdexcept>
#include <cstdlib>  // Add this for EXIT_FAILURE

#include "Materials/Mechanics/Elastic.h"

Elastic::Elastic(string dimensions, H5IO& H5File, int iSet)
    : BaseMechanics(dimensions) {

    try {
        string elasticity = H5File.ReadString("Materials/Material_" + to_string(iSet)+"/Elastic/Elasticity");
        string isotropy2 = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Elastic/Isotropy");

        if (isotropy2=="Isotropic"){

            double Emod = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/Emod");
            double nu = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/nu");

            double ho = Emod*nu/((1+nu)*(1-2*nu));
            double uo = Emod/(2*(1+nu));

            InitializeIsoElasticityMatrix(elasticity, Emod, nu, ho, uo);

        } else if (isotropy2=="Cubic"){

            double C11 = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/C11");
            double C12 = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/C12");
            double C44 = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/C44");

            InitializeCubicElasticityMatrix(elasticity, C11, C12, C44);

        } else {
            throw std::invalid_argument("Undefined material isotropy: < " + isotropy2 + " >");
        }

    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}

void Elastic::InitializeIsoElasticityMatrix(const string& elasticity, double Emod, double nu, double ho, double uo) {

    if (dims != "3D")
        throw std::invalid_argument("Invalid dimension: < " + dims + " > for < " + elasticity + " > elasticity.");

    try {
        if (elasticity == "3D") {
            DMatx = Matd6x6(Matd6x6::Zero());
            auto& mat = std::get<Matd6x6>(DMatx);

            mat << ho + 2 * uo, ho, ho, 0, 0, 0,
                   ho, ho + 2 * uo, ho, 0, 0, 0,
                   ho, ho, ho + 2 * uo, 0, 0, 0,
                   0, 0, 0, 2 * uo, 0, 0,
                   0, 0, 0, 0, 2 * uo, 0,
                   0, 0, 0, 0, 0, 2 * uo;

        } else if (elasticity == "PlaneStrain" || elasticity == "PlaneStress") {

            if (dims != "2D")
                throw std::invalid_argument("Invalid dimension: < " + dims + " > for < " + elasticity + " > elasticity.");

            DMatx = Matd3x3(Matd3x3::Zero());
            auto& mat = std::get<Matd3x3>(DMatx);

            double param = Emod * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

            if (elasticity == "PlaneStrain") {
                mat << param, param * nu / (1.0 - nu), 0,
                       param * nu / (1.0 - nu), param, 0,
                       0, 0, param * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));

            } else if (elasticity == "PlaneStress") {
                param = Emod / (1.0 - nu * nu);
                mat << param, param * nu, 0,
                       param * nu, param, 0,
                       0, 0, param * (1.0 - nu) / 2.0;
            }
        } else {
            throw std::invalid_argument("Undefined material elasticity: < " + elasticity + " >");
        }
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}

void Elastic::InitializeCubicElasticityMatrix(const string& elasticity, double C11, double C12, double C44) {
    try {
        if (elasticity == "3D") {
            DMatx = Matd6x6(Matd6x6::Zero());
            auto& mat = std::get<Matd6x6>(DMatx);

            mat << C11, C12, C12, 0, 0, 0,
                   C12, C11, C12, 0, 0, 0,
                   C12, C12, C11, 0, 0, 0,
                   0, 0, 0, C44, 0, 0,
                   0, 0, 0, 0, C44, 0,
                   0, 0, 0, 0, 0, C44;

        } else if (elasticity == "PlaneStrain" || elasticity == "PlaneStress") {
            DMatx = Matd3x3(Matd3x3::Zero());
            auto& mat = std::get<Matd3x3>(DMatx);

            double param = C11 - C12;
            mat << C11, C12, 0,
                   C12, C11, 0,
                   0, 0, param / 2;
        } else {
            throw std::invalid_argument("Undefined material elasticity: < " + elasticity + " >");
        }
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}


T_DMatx Elastic::getDMatx() const{

    return DMatx;
}
