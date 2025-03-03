
template <>
inline void IsoHard::UHard<PowerLaw>(const double& eqpl, double& syield, double& hard){

    syield = sig_y0 +  K_hard * std::pow(eqpl, n_pow);
    double eps = std::max(eqpl, 1.0e-12);
    hard = K_hard*n_pow*pow(eps, n_pow-1);
}

template <>
inline void IsoHard::UHard<Voce>(const double& eqpl, double& syield, double& hard){

    // double sigma0 = 200.0;
    // double sigmaS = 300.0;
    // double h = 50.0;
    // syield = sigma0 + (sigmaS - sigma0) * (1 - std::exp(-h * eqpl));
    // hard = h * (sigmaS - sigma0) * std::exp(-h * eqpl);
}

template <typename HardeningLaw>
void IsoHard::RM3D(ColVecd6& deps, ColVecd6& sig, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd6& eps_e_old, const ColVecd6& eps_p_old, const double& eps_eq_old, const int iStep){

    double Ebulk3 = Emod/(1.0 - 2.0*nu);

    // Elastic strain
    eps_e = eps_e_old + deps;

    // Trial stress
    ColVecd6 sig_trial = std::get<Matd6x6>(DMatx_e)*eps_e;

    //Mises stress
    double sig_trial_eq = Mises3D(sig_trial);

    // Initialize plastic multiplier and increment
    double p = eps_eq_old;
    double deqpl = 0;

    // Current yield stress and hardening modulus
    double sYield0, sYield, hard; 
    UHard<HardeningLaw>(p, sYield0, hard); 
    sYield = sYield0;

    // Yield function 
    double f_yield = sig_trial_eq - sYield0*(1.0 + 1.0e-6);

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        sig_h = (1.0/3.0)*sig_trial.segment<3>(0).sum();
        std::get<Matd6x6>(DMatx_ep) = std::get<Matd6x6>(DMatx_e);

    } else { // --> Plastic step

        // Iteration counter
        int nIter_RM = 0;

        // Deviatoric stress
        ColVecd6 sig_trial_dev = sig_trial - (1.0/3.0)*sig_trial.segment<3>(0).sum()*I6; 
        ColVecd6 N_tr = (3.0/2.0)*sig_trial_dev/sig_trial_eq; // Plastic flow direction

        while(abs(f_yield > tol)){

            // Update Iteration counter
            nIter_RM++;

            // Update plastic strain increment 
            deqpl += f_yield/(3*uo + hard);
            p = eps_eq_old + deqpl;
            UHard<HardeningLaw>(p, sYield, hard);
            f_yield = sig_trial_eq - 3*uo*deqpl - sYield;

            if(nIter_RM > max_iter){

                std::ostringstream oss;
                oss << "\nReturn mapping did not converge at step: " << iStep << "\n"
                    << "Reached maximum iterations: " << nIter_RM << "\n"
                    << "Final yield function value: " << f_yield << "\n"
                    << "Plastic strain increment: " << deqpl << "\n"
                    << "Equivalent plastic strain: " << p << "\n";
                
                throw std::runtime_error(oss.str());
            }
        }

        // --------- CONVERGED ---> Update variables

        eps_eq = p;             // Equivalent plastic strain
        eps_p += deqpl*N_tr;    // Plastic strain tensor
        eps_e -= deqpl*N_tr;    // Elastic strain tensor
        
        sig = std::get<Matd6x6>(DMatx_e)*eps_e;  // Stress tensor
        sig_eq = Mises3D(sig);  // Von Mises stress
        sig_h = (1.0/3.0)*sig_trial.segment<3>(0).sum(); // Hydrostatic stress

        // Tangent stiffness matrix

        double effg = uo * sYield/sig_trial_eq;
        double efflam = (Ebulk3 - 2.0*effg)/3.0;
        double effhard = 3.0*uo * hard/(3.0*uo + hard) - (3*effg);

        std::get<Matd6x6>(DMatx_ep).setZero();

        std::get<Matd6x6>(DMatx_ep).topLeftCorner<3, 3>().setConstant(efflam);

        for (int i = 0; i < 3; ++i) {
            std::get<Matd6x6>(DMatx_ep)(i, i) += 2.0*effg;
            std::get<Matd6x6>(DMatx_ep)(i+3, i+3) += effg;
        }

        std::get<Matd6x6>(DMatx_ep) += effhard*(2.0/3.0)*N_tr*(2.0/3.0)*N_tr.transpose();

    }
}

/// @brief Specialization for plane strain.
template <>
inline double IsoHard::Mises2D<PlaneStrain>(const ColVecd3& sig2D){

    double sx = sig2D(0); 
    double sy = sig2D(1); 
    double txy = sig2D(2);  
    double sz = nu * (sx + sy);

    double term = 0.5 * (pow(sx - sy, 2) +
                         pow(sy - sz, 2) +
                         pow(sz - sx, 2)) + 3 * txy * txy;

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

/// @brief Specialization for plane strain.
template <>
inline double IsoHard::Shydro2D<PlaneStrain>(const ColVecd3& sig2D){

    double sx = sig2D(0); 
    double sy = sig2D(1); 
    double sz = nu * (sx + sy);

    double hydrostatic = 1.0/3.0* (sx + sy + sz);

    return hydrostatic;
}

/// @brief Specialization for plane stress.
template <>
inline double IsoHard::Shydro2D<PlaneStress>(const ColVecd3& sig2D){

    double sx = sig2D(0); 
    double sy = sig2D(1); 

    double hydrostatic = 1.0/3.0* (sx + sy);

    return hydrostatic;
}

/// @brief Specialization 
template <typename AnalysisType, typename HardeningLaw>
void IsoHard::RM2D(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const int iStep){

    double Ebulk3 = Emod/(1.0 - 2.0*nu);

    // Elastic strain
    eps_e = eps_e_old + deps;

    // Trial stress
    ColVecd3 sig_trial = std::get<Matd3x3>(DMatx_e)*eps_e;

    //Mises stress
    double sig_trial_eq = Mises2D<AnalysisType>(sig_trial);

    // Initialize plastic multiplier and increment
    double p = eps_eq_old;
    double deqpl = 0;

    // Current yield stress and hardening modulus
    double sYield0, sYield, hard; 
    UHard<HardeningLaw>(p, sYield0, hard); 
    sYield = sYield0;

    // Yield function 
    double f_yield = sig_trial_eq - sYield0*(1.0 + 1.0e-6);

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        sig_h = Shydro2D<AnalysisType>(sig_trial);
        std::get<Matd3x3>(DMatx_ep) = std::get<Matd3x3>(DMatx_e);

    } else { // --> Plastic step

        // Iteration counter
        int nIter_RM = 0;

        // Deviatoric stress

        ColVecd3 sig_trial_dev = sig_trial - Shydro2D<AnalysisType>(sig_trial)*I3; 
        ColVecd3 N_tr = (3.0/2.0)*sig_trial_dev/sig_trial_eq; // Plastic flow direction

        while(abs(f_yield > tol)){

            // Update Iteration counter
            nIter_RM++;

            // Update plastic strain increment 
            deqpl += f_yield/(3*uo + hard);
            p = eps_eq_old + deqpl;
            UHard<HardeningLaw>(p, sYield, hard);
            f_yield = sig_trial_eq - 3*uo*deqpl - sYield;

            if(nIter_RM > max_iter){

                std::ostringstream oss;
                oss << "\nReturn mapping did not converge at step: " << iStep << "\n"
                    << "Reached maximum iterations: " << nIter_RM << "\n"
                    << "Final yield function value: " << f_yield << "\n"
                    << "Plastic strain increment: " << deqpl << "\n"
                    << "Equivalent plastic strain: " << p << "\n";
                
                throw std::runtime_error(oss.str());
            }
        }

        // --------- CONVERGED ---> Update variables

        eps_eq = p;             // Equivalent plastic strain
        eps_p += deqpl*N_tr;    // Plastic strain tensor
        eps_e -= deqpl*N_tr;    // Elastic strain tensor
        
        sig = std::get<Matd3x3>(DMatx_e)*eps_e;  // Stress tensor
        sig_eq = Mises2D<AnalysisType>(sig);  // Von Mises stress
        sig_h = Shydro2D<AnalysisType>(sig_trial); // Hydostatic stress

        // Tangent stiffness matrix

        double effg = uo * sYield/sig_trial_eq;
        double efflam = (Ebulk3 - 2.0*effg)/3.0;
        double effhard = 3.0*uo * hard/(3.0*uo + hard) - (3*effg);

        std::get<Matd3x3>(DMatx_ep).setZero();

        std::get<Matd3x3>(DMatx_ep).topLeftCorner<2, 2>().setConstant(efflam);

        for (int i = 0; i < 2; ++i) {
            std::get<Matd3x3>(DMatx_ep)(i, i) += 2.0*effg;
            std::get<Matd3x3>(DMatx_ep)(i+1, i+1) += effg;
        }

        std::get<Matd3x3>(DMatx_ep) += effhard*(2.0/3.0)*N_tr*(2.0/3.0)*N_tr.transpose();

    }
}