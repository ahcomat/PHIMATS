#include <iostream>
#include <cmath>

#include "Materials/Mechanics/IsoHard.h"
#include "Materials/Mechanics/IsoHard.tpp"  // Include template implementation

IsoHard::IsoHard(string dimensions, H5IO& H5File, int iSet)
    : LinearElastic(dimensions, H5File, iSet) {

    try {
        // Read plasticity and hardening law
        Platicity = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/Plasticity");
        hardLaw = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/HardeningLaw");
        sig_y0 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/sig_y0");

        // Handle supported hardening laws
        if (hardLaw == "PowerLaw") {

            K_hard = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/K_hard");
            n_pow = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/n_pow");

        } else {

            throw invalid_argument("Undefined hardening law < " + hardLaw + " >");

        }
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        cerr << "Terminating!" << endl;
        exit(EXIT_FAILURE);
    }

    if (dims == "3D"){

        DMatx_ep = Matd6x6(Matd6x6::Zero());

    } else if (dims == "2D") {

        if (analysisType == "PlaneStrain") {
            analysis2D = AnalysisType::PlaneStrain;
        } else if (analysisType == "PlaneStress") {
            analysis2D = AnalysisType::PlaneStress;
        } else {
            throw std::invalid_argument("Invalid 2D analysis type: " + analysisType);
        }

         DMatx_ep =  Matd3x3(Matd3x3::Zero());

    } else {

        throw std::invalid_argument("Invalid dimension: < " + dims + " > for < " + analysisType + " > analysis.");

    }
}

double IsoHard::R_pow(const double& eps_eq){
    return K_hard*pow(eps_eq, n_pow);
}

double IsoHard::dR_pow(const double& eps_eq){

    double eps = std::max(eps_eq, 1.0e-12);

    return K_hard*n_pow*pow(eps, n_pow-1);
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

void IsoHard::ReturnMapping3D(ColVecd6& sig, ColVecd6& eps, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, const int iStep){

    // Elastic strain
    eps_e = eps - eps_p;

    // Trial stress
    ColVecd6 sig_trial = std::get<Matd6x6>(DMatx_e)*eps_e;

    //Mises stress
    double sig_trial_eq = Mises3D(sig_trial);

    // Yield function 
    double f_yield = sig_trial_eq - R_pow(eps_eq) - sig_y0;

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        std::get<Matd6x6>(DMatx_ep) = std::get<Matd6x6>(DMatx_e);

    } else { // --> Plastic step

        // Iteration counter
        int nIter_RM = 0;

        // Initialize plastic multiplier
        double p = eps_eq;
        double delta_p = 0;

        while(abs(f_yield > tol)){

            // Update Iteration counter
            nIter_RM++;

            // Update plastic strain increment 
            delta_p += f_yield/(3*uo + dR_pow(p));
            p = eps_eq + delta_p;
            f_yield = sig_trial_eq - 3*uo*delta_p - R_pow(p) - sig_y0;

            if(nIter_RM > max_iter){

                std::ostringstream oss;
                oss << "\nReturn mapping did not converge at step: " << iStep << "\n"
                    << "Reached maximum iterations: " << nIter_RM << "\n"
                    << "Final yield function value: " << f_yield << "\n"
                    << "Plastic strain increment: " << delta_p << "\n"
                    << "Equivalent plastic strain: " << p << "\n";
                
                throw std::runtime_error(oss.str());
            }
        }

        // Update variables
        ColVecd6 sig_trial_dev = sig_trial - (1.0/3.0)*sig_trial.segment<3>(0).sum()*I6; // Deviatoric stress
        ColVecd6 N_tr = (3.0/2.0)*sig_trial_dev/sig_trial_eq; // Plastic flow direction

        eps_eq = p;             // Equivalent plastic strain
        eps_p += delta_p*N_tr;  // Plastic strain tensor
        eps_e = eps - eps_p;    // Elastic strain tensor
        
        sig = std::get<Matd6x6>(DMatx_e)*eps_e;  // Stress tensor
        sig_eq = Mises3D(sig);  // Von Mises stress

        // Tangent stiffness matrix
        double Hmod =  dR_pow(eps_eq);
        std::get<Matd6x6>(DMatx_ep) = std::get<Matd6x6>(DMatx_e) - 
                                      (2*uo/(1+Hmod/(3*uo)))*N_tr*N_tr.transpose();
    }
}

/// @brief Select appropriate template specialization 
void IsoHard::ReturnMapping2D(ColVecd3& sig, ColVecd3& eps, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, const int iStep){

    if (analysis2D == AnalysisType::PlaneStrain) {
        return ReturnMapping2D<PlaneStrain>(sig, eps, eps_e, eps_p, eps_eq, sig_eq, iStep);
    } else if (analysis2D == AnalysisType::PlaneStress) {
        return ReturnMapping2D<PlaneStress>(sig, eps, eps_e, eps_p, eps_eq, sig_eq, iStep);
    } else {
        throw std::logic_error("Unhandled analysis type.");
    }

}

T_DMatx IsoHard::getDMatx() const{

    return DMatx_ep;
}


