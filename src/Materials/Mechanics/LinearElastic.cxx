#include <iostream>
#include <stdexcept>
#include <cstdlib>  // Add this for EXIT_FAILURE

#include "Materials/Mechanics/LinearElastic.h"

LinearElastic::LinearElastic(string dimensions, H5IO& H5File, int iSet, Logger& logger)
    : BaseMechanics(dimensions, logger) {

    try {
        analysisType = H5File.ReadString("Materials/Material_" + to_string(iSet)+"/Elastic/AnalysisType");
        string isotropy = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Elastic/Isotropy");

        if (isotropy=="Isotropic"){

            Emod = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/Emod");
            nu = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/nu");

            ho = Emod*nu/((1+nu)*(1-2*nu));
            uo = Emod/(2*(1+nu));

            InitializeIsoElasticityMatrix(analysisType, Emod, nu, ho, uo);

        } else if (isotropy=="Cubic"){

            C11 = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/C11");
            C12 = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/C12");
            C44 = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/Elastic/C44");

            InitializeCubicElasticityMatrix(analysisType, C11, C12, C44);

        } else {
            throw std::invalid_argument("Undefined material isotropy: < " + isotropy + " >");
        }

    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}

void LinearElastic::InitializeIsoElasticityMatrix(const string& analysisType, double Emod, double nu, double ho, double uo) {

    try {
        if (analysisType == "3D") {

            if (dims != "3D"){
                throw std::invalid_argument("Invalid dimension: < " + dims + " > for < " + analysisType + " > analysis.");
                }
                
            DMatx_e = Matd6x6(Matd6x6::Zero());
            auto& mat = std::get<Matd6x6>(DMatx_e);

            mat << ho + 2 * uo, ho, ho, 0, 0, 0,
                   ho, ho + 2 * uo, ho, 0, 0, 0,
                   ho, ho, ho + 2 * uo, 0, 0, 0,
                   0, 0, 0, uo, 0, 0,
                   0, 0, 0, 0, uo, 0,
                   0, 0, 0, 0, 0, uo;

        } else if (analysisType == "PlaneStrain" || analysisType == "PlaneStress") {

            if (dims != "2D")
                throw std::invalid_argument("Invalid dimension: < " + dims + " > for < " + analysisType + " > analysis.");

            DMatx_e = Matd3x3(Matd3x3::Zero());
            auto& mat = std::get<Matd3x3>(DMatx_e);

            double param = Emod * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

            if (analysisType == "PlaneStrain") {
                mat << param, param * nu / (1.0 - nu), 0,
                       param * nu / (1.0 - nu), param, 0,
                       0, 0, param * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));

            } else if (analysisType == "PlaneStress") {
                param = Emod / (1.0 - nu * nu);
                mat << param, param * nu, 0,
                       param * nu, param, 0,
                       0, 0, param * (1.0 - nu) / 2.0;
            }
        } else {
            throw std::invalid_argument("Undefined material analysis type: < " + analysisType + " >");
        }
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}

void LinearElastic::InitializeCubicElasticityMatrix(const string& analysisType, double C11, double C12, double C44) {
    try {
        if (analysisType == "3D") {
            DMatx_e = Matd6x6(Matd6x6::Zero());
            auto& mat = std::get<Matd6x6>(DMatx_e);

            mat << C11, C12, C12, 0, 0, 0,
                   C12, C11, C12, 0, 0, 0,
                   C12, C12, C11, 0, 0, 0,
                   0, 0, 0, C44, 0, 0,
                   0, 0, 0, 0, C44, 0,
                   0, 0, 0, 0, 0, C44;

        } else if (analysisType == "PlaneStrain" || analysisType == "PlaneStress") {
            DMatx_e = Matd3x3(Matd3x3::Zero());
            auto& mat = std::get<Matd3x3>(DMatx_e);

            double param = C11 - C12;
            mat << C11, C12, 0,
                   C12, C11, 0,
                   0, 0, param / 2;
        } else {
            throw std::invalid_argument("Undefined material analysis type: < " + analysisType + " >");
        }
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }
}

T_DMatx LinearElastic::getDMatx() const{

    return DMatx_e;
}
