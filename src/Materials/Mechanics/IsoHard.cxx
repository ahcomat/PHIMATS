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

            rho_0 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/rho_0");
            M = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/M");
            b = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/b");
            alpha = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/alpha");
            k1 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/k1");
            k2 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/k2");
            rho_s = pow(k1/k2, 2.0);
            C_prime = (k1/k2) - sqrt(rho_0);

        } else if (hardLaw == "Voce") {

            hardening = HardeningLaw::Voce;

        } else if (hardLaw == "KME") {

            rho_0 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/rho_0");
            M = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/M");
            b = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/b");
            alpha = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/alpha");
            k1 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/k1");
            k2 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/k2");
            rho_s = pow(k1/k2, 2.0);
            C_prime = (k1/k2) - sqrt(rho_0);
            hardening = HardeningLaw::KME;

        } else {

            throw invalid_argument("Undefined hardening law < " + hardLaw + " >");

        }
    } catch (const std::runtime_error& e) {
        logger.log("\nException caught in IsoHard constructor\n", "Error", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
    }

    if (dims == "3D"){

        if (hardLaw == "PowerLaw") {
            // Initialize the function pointer
            selectedRM3D = &IsoHard::RM3D<PowerLaw>;
        } else if (hardLaw == "Voce") {
            selectedRM3D = &IsoHard::RM3D<Voce>;
        } else if(hardLaw == "KME") {
            selectedRM3D = &IsoHard::RM3D<KME>;
        }

        CMatx_ep = Matd6x6(Matd6x6::Zero());

    } else if (dims == "2D") {

        if (analysisType == "PlaneStrain") {
            analysis2D = AnalysisType::PlaneStrain;
            if (hardLaw == "PowerLaw") {
                // Initialize the function pointer
                selectedRM2D = &IsoHard::RM2D<PlaneStrain, PowerLaw>;
            } else if (hardLaw == "Voce") {
                selectedRM2D = &IsoHard::RM2D<PlaneStrain, Voce>;
            } else if (hardLaw == "KME") {
                selectedRM2D = &IsoHard::RM2D<PlaneStrain, KME>;
            }
        } else if (analysisType == "PlaneStrainPFF") {
            analysis2D = AnalysisType::PlaneStress;
            if (hardLaw == "PowerLaw") {
                // Initialize the function pointer
                selectedRM2DPFF = &IsoHard::RM2DPFF<PlaneStrain, PowerLaw>;
            } else if (hardLaw == "Voce") {
                selectedRM2DPFF = &IsoHard::RM2DPFF<PlaneStrain, Voce>;
            } else if (hardLaw == "KME") {
                selectedRM2DPFF = &IsoHard::RM2DPFF<PlaneStrain, KME>;
            }
        } else if (analysisType == "AxiSymmetric") {
            analysis2D = AnalysisType::AxiSymmetric;
            if (hardLaw == "PowerLaw") {
                selectedRMAxi = &IsoHard::RMAxi<PowerLaw>;
            } else if (hardLaw == "Voce") {
                selectedRMAxi = &IsoHard::RMAxi<Voce>;
            } else if (hardLaw == "KME") {
                selectedRMAxi = &IsoHard::RMAxi<KME>;
            }
            CMatx_ep = Matd4x4(Matd4x4::Zero()); // Initialize as 4x4
        } else if (analysisType == "PlaneStress") {
            analysis2D = AnalysisType::PlaneStress;
            if (hardLaw == "PowerLaw") {
                // Initialize the function pointer
                selectedRM2D = &IsoHard::RM2D<PlaneStress, PowerLaw>;
            } else if (hardLaw == "Voce") {
                selectedRM2D = &IsoHard::RM2D<PlaneStress, Voce>;
            } else if (hardLaw == "KME") {
                selectedRM2D = &IsoHard::RM2D<PlaneStress, KME>;
            }
        } else {
            throw std::invalid_argument("Invalid 2D analysis type: " + analysisType);
        }

         CMatx_ep =  Matd3x3(Matd3x3::Zero());

    } else {

        throw std::invalid_argument("Invalid dimension: < " + dims + " > for < " + analysisType + " > analysis.");

    }

    // Initalize to elastic values for initial stiffness matrix.
    CMatx_ep = CMatx_e;
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

double IsoHard::MisesAxi(const ColVecd4& sig) {
    
    double s_rr = sig(0); 
    double s_zz = sig(1); 
    double s_hoop = sig(2); 
    double t_rz = sig(3); 

    double term = 0.5 * (std::pow(s_rr - s_zz, 2) + 
                         std::pow(s_zz - s_hoop, 2) + 
                         std::pow(s_hoop - s_rr, 2)) + 
                  3.0 * std::pow(t_rz, 2);

    return std::sqrt(term);
}

void IsoHard::ReturnMapping3D(ColVecd6& deps, ColVecd6& sig, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd6& eps_e_old, const ColVecd6& eps_p_old, const double& eps_eq_old, const int iStep){

    // Ensure selectedRM3D is valid
    if (!selectedRM3D) {
        throw std::runtime_error("ReturnMapping3D function pointer is not set.");
    }

    (this->*selectedRM3D)(deps, sig, eps_e, eps_p, eps_eq, sig_eq, sig_h, rho, eps_e_old, eps_p_old, eps_eq_old, iStep);
}

void IsoHard::ReturnMapping2D(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep){

    // Ensure selectedRM2D is valid
    if (!selectedRM2D) {
        throw std::runtime_error("ReturnMapping2D function pointer is not set. Make sure you are not using a PFF material model.");
    }

    (this->*selectedRM2D)(deps, sig, eps_e, eps_p, eps_eq, sig_eq, sig_h, rho, eps_e_old, eps_p_old, eps_eq_old, iStep);
}

void IsoHard::ReturnMapping2D_PFF(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep, const double gPhi_d, const double& wp_old, double& wp){

    // Ensure selectedRM2DPFF is valid
    if (!selectedRM2DPFF) {
        throw std::runtime_error("ReturnMapping2DPFF function pointer is not set.");
    }

    (this->*selectedRM2DPFF)(deps, sig, eps_e, eps_p, eps_eq, sig_eq, sig_h, rho, eps_e_old, eps_p_old, eps_eq_old, iStep, gPhi_d, wp_old, wp);
}

void IsoHard::ReturnMappingAxi(ColVecd4& deps, ColVecd4& sig, ColVecd4& eps_e, ColVecd4& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd4& eps_e_old, const ColVecd4& eps_p_old, const double& eps_eq_old, const int iStep){

    // Ensure selectedRM3D is valid
    if (!selectedRMAxi) {
        throw std::runtime_error("ReturnMappingAxi function pointer is not set. Make sure you are not using a PFF material model.");
    }

    (this->*selectedRMAxi)(deps, sig, eps_e, eps_p, eps_eq, sig_eq, sig_h, rho, eps_e_old, eps_p_old, eps_eq_old, iStep);
}

T_DMatx IsoHard::getCMatx() const{

    return CMatx_ep;
}


