#include <iostream>
#include <cmath>

#include "Materials/Mechanics/IsoHard.h"

IsoHard::IsoHard(string dimensions, H5IO& H5File, int iSet)
    : LinearElastic(dimensions, H5File, iSet) {

    try {
        // Read plasticity and hardening law
        Platicity = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/Plasticity");
        HardLaw = H5File.ReadString("Materials/Material_" + to_string(iSet) + "/Plastic/HardeningLaw");
        sig_y0 = H5File.ReadScalar("Materials/Material_" + to_string(iSet) + "/Plastic/sig_y0");

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

    if (dims == "3D"){

        DMatx_ep = Matd6x6(Matd6x6::Zero());

    } else if (dims == "2D") {

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

    double term = 0.5 * (pow(sig3D(0) - sig3D(1), 2) + 
                         pow(sig3D(1) - sig3D(2), 2) + 
                         pow(sig3D(2) - sig3D(0), 2)) + 
                  3 * (pow(sig3D(3),2) + pow(sig3D(4),2) + pow(sig3D(5),2));

    return sqrt(term);
}

// double PEEQ3D(ColVecd6& eps3D){
//     double term = 2.0/3.0 * (eps3D)
// }

void IsoHard::ReturnMapping3D(ColVecd6& sig, ColVecd6& eps, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, int iStep){

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
            nIter_RM += 1;

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
        ColVecd6 sig_trial_dev = sig_trial - (1/3)*sig_trial.segment<3>(0).sum()*I6; // Deviatoric stress
        ColVecd6 N_tr = (3/2)*sig_trial_dev/sig_trial_eq; // Plastic flow direction

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

T_DMatx IsoHard::getDMatx() const{

    return DMatx_ep;
}
