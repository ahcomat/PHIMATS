#include <iostream>
#include <cmath>

#include "Materials/Mechanics/IsoHard.h"
#include "Materials/Mechanics/IsoHard.tpp"  // Include template implementation

IsoHard::IsoHard(string dimensions, H5IO& H5File, int iSet, Logger& logger)
    : LinearElastic(dimensions, H5File, iSet, logger) {

    try {
        // Read plasticity and hardening law
        Platicity = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/Plasticity");
        hardLaw = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/HardeningLaw");
        sig_y0 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/sig_y0");

        // Handle supported hardening laws
        if (hardLaw == "PowerLaw") {

            K_hard = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/K_hard");
            n_pow = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/n_pow");
            hardening = HardeningLaw::PowerLaw;

        } else if (hardLaw == "Voce") {

            hardening = HardeningLaw::Voce;

        } else {

            throw invalid_argument("Undefined hardening law < " + hardLaw + " >");

        }
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }

    if (dims == "3D"){

        if (hardLaw == "PowerLaw") {
            // Initialize the function pointer
            selectedRM3D = &IsoHard::RM3D<PowerLaw>;
        } else if (hardLaw == "Voce") {
            selectedRM3D = &IsoHard::RM3D<Voce>;
        }

        DMatx_ep = Matd6x6(Matd6x6::Zero());
        double dummy;

    } else if (dims == "2D") {

        if (analysisType == "PlaneStrain") {
            analysis2D = AnalysisType::PlaneStrain;
            if (hardLaw == "PowerLaw") {
                // Initialize the function pointer
                selectedRM2D = &IsoHard::RM2D<PlaneStrain, PowerLaw>;
            } else if (hardLaw == "Voce") {
                selectedRM2D = &IsoHard::RM2D<PlaneStrain, Voce>;
            }
        } else if (analysisType == "PlaneStress") {
            analysis2D = AnalysisType::PlaneStress;
            if (hardLaw == "PowerLaw") {
                // Initialize the function pointer
                selectedRM2D = &IsoHard::RM2D<PlaneStress, PowerLaw>;
            } else if (hardLaw == "Voce") {
                selectedRM2D = &IsoHard::RM2D<PlaneStress, Voce>;
            }
        } else {
            throw std::invalid_argument("Invalid 2D analysis type: " + analysisType);
        }

         DMatx_ep =  Matd3x3(Matd3x3::Zero());

    } else {

        throw std::invalid_argument("Invalid dimension: < " + dims + " > for < " + analysisType + " > analysis.");

    }

    // Initalize to elastic values for initial stiffness matrix.
    DMatx_ep = DMatx_e;
}

double IsoHard::Mises3D(const ColVecd6& sig3D){

    double sx = sig3D(0); 
    double sy = sig3D(1); 
    double sz = sig3D(2); 
    double txy = sig3D(3); 
    double tyz = sig3D(4); 
    double tzx = sig3D(5); 

    // Compute von Mises stress term
    double term = 0.5 * (pow(sx - sy, 2) +
                         pow(sy - sz, 2) +
                         pow(sz - sx, 2)) +
                  3 * (pow(txy, 2) + pow(tyz, 2) + pow(tzx, 2));

    return sqrt(term);
}

void IsoHard::ReturnMapping3D(ColVecd6& deps, ColVecd6& sig, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, const ColVecd6& eps_e_old, const ColVecd6& eps_p_old, const double& eps_eq_old, const int iStep){

    // Ensure selectedRM3D is valid
    if (!selectedRM3D) {
        throw std::runtime_error("ReturnMapping3D function pointer is not set.");
    }

    (this->*selectedRM3D)(deps, sig, eps_e, eps_p, eps_eq, sig_eq, eps_e_old, eps_p_old, eps_eq_old, iStep);
}

/// @brief Select appropriate template specialization 
void IsoHard::ReturnMapping2D(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep){

        // Ensure selectedRM3D is valid
    if (!selectedRM2D) {
        throw std::runtime_error("ReturnMapping2D function pointer is not set.");
    }

    (this->*selectedRM2D)(deps, sig, eps_e, eps_p, eps_eq, sig_eq, eps_e_old, eps_p_old, eps_eq_old, iStep);
}

T_DMatx IsoHard::getDMatx() const{

    return DMatx_ep;
}


