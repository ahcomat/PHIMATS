/// @brief Specialization for plane strain.
template <>
inline double IsoHard::Mises2D<PlaneStrain>(const ColVecd3& sig2D){

    double sx = sig2D(0); 
    double sy = sig2D(1); 
    double txy = sig2D(2);  

    // Poisson's ratio from Lamé constants
    double nu = ho / (2 * (ho + uo));

    double hydrostatic_term = nu * (sx + sy);
    double term = sx * sx - sx * sy + sy * sy + 3 * txy * txy - hydrostatic_term * hydrostatic_term;

    return sqrt(term);
}

/// @brief Specialization for plane stress.
template <>
inline double IsoHard::Mises2D<PlaneStress>(const ColVecd3& sig2D){

    double sx = sig2D(0); 
    double sy = sig2D(1); 
    double txy = sig2D(2); 

    double term = sx * sx - sx * sy + sy * sy + 3 * txy * txy;

    return sqrt(term);
}

// /// @brief Specialization 
template <typename AnalysisType>
void IsoHard::ReturnMapping2D(ColVecd3& sig, ColVecd3& eps, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, int iStep){

    // Elastic strain
    eps_e = eps - eps_p;

    // Trial stress
    ColVecd3 sig_trial = std::get<Matd3x3>(DMatx_e)*eps_e;

    double sig_trial_eq = Mises2D<AnalysisType>(sig_trial);

    // Yield function 
    double f_yield = sig_trial_eq - R_pow(eps_eq) - sig_y0;

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        std::get<Matd3x3>(DMatx_ep) = std::get<Matd3x3>(DMatx_e);

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
        ColVecd3 sig_trial_dev = sig_trial - (1.0/3.0)*sig_trial.segment<2>(0).sum()*I3; // Deviatoric stress
        ColVecd3 N_tr = (3.0/2.0)*sig_trial_dev/sig_trial_eq; // Plastic flow direction

        eps_eq = p;             // Equivalent plastic strain
        eps_p += delta_p*N_tr;  // Plastic strain tensor
        eps_e = eps - eps_p;    // Elastic strain tensor
        
        sig = std::get<Matd3x3>(DMatx_e)*eps_e;  // Stress tensor
        sig_eq = Mises2D<AnalysisType>(sig);  // Von Mises stress

        // Tangent stiffness matrix
        double Hmod =  dR_pow(eps_eq);
        std::get<Matd3x3>(DMatx_ep) = std::get<Matd3x3>(DMatx_e) - 
                                      (2*uo/(1+Hmod/(3*uo)))*N_tr*N_tr.transpose();
    }
}