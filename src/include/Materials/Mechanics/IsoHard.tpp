
template <>
inline void IsoHard::UHard<PowerLaw>(const double& eqpl, double& syield, double& rho, double& hard){

    syield = K_hard * std::pow(eps_0 + eqpl, n_pow);
    hard = K_hard * n_pow * std::pow(eps_0 + eqpl, n_pow - 1.0);

    // Passive evolution of dislocation density, i.e. does not affect hardening
    double param = (k1/k2) - C_prime*exp(-(M*k2/2)*eqpl);
    rho = pow(param, 2);
    rho = rho/rho_s; // Normalizing
}

template <>
inline void IsoHard::UHard<Voce>(const double& eqpl, double& syield, double& rho, double& hard){

    // double sigma0 = 200.0;
    // double sigmaS = 300.0;
    // double h = 50.0;
    // syield = sigma0 + (sigmaS - sigma0) * (1 - std::exp(-h * eqpl));
    // hard = h * (sigmaS - sigma0) * std::exp(-h * eqpl);
}

template <>
inline void IsoHard::UHard<KME>(const double& eqpl, double& syield, double& rho, double& hard){


    /* 
    Note that here calculate the dislocation density, use it in evaluating
    equations then normalized it at the end
    */
    double param = (k1/k2) - C_prime*exp(-(M*k2/2)*eqpl);
    rho = pow(param, 2);
    syield = sig_y0 +  M*alpha*uo*b*sqrt(rho); 
    hard = pow(M, 2)*alpha*uo*b/2*(k1 - k2*sqrt(rho));
    rho = rho/rho_s; // Normalizing
}

template <typename HardeningLaw>
void IsoHard::RM3D(ColVecd6& deps, ColVecd6& sig, ColVecd6& eps_e, ColVecd6& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd6& eps_e_old, const ColVecd6& eps_p_old, const double& eps_eq_old, const int iStep){

    double Ebulk3 = Emod/(1.0 - 2.0*nu);

    // Elastic strain
    eps_e = eps_e_old + deps;

    // Plastic strain
    eps_p = eps_p_old;

    // Trial stress
    ColVecd6 sig_trial = std::get<Matd6x6>(CMatx_e)*eps_e;

    //Mises stress
    double sig_trial_eq = Mises3D(sig_trial);

    // Initialize plastic multiplier and increment
    double p = eps_eq_old;
    double deqpl = 0;

    // Current yield stress and hardening modulus
    double sYield0, sYield, hard; 
    UHard<HardeningLaw>(p, sYield0, rho, hard); 
    sYield = sYield0;

    // Yield function 
    double f_yield = sig_trial_eq - sYield0*(1.0 + 1.0e-6);

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        sig_h = (1.0/3.0)*sig_trial.segment<3>(0).sum();
        std::get<Matd6x6>(CMatx_ep) = std::get<Matd6x6>(CMatx_e);

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
            UHard<HardeningLaw>(p, sYield, rho, hard);
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
        
        sig = std::get<Matd6x6>(CMatx_e)*eps_e;  // Stress tensor
        sig_eq = Mises3D(sig);  // Von Mises stress
        sig_h = (1.0/3.0)*sig_trial.segment<3>(0).sum(); // Hydrostatic stress

        // Tangent stiffness matrix

        double effg = uo * sYield/sig_trial_eq;
        double efflam = (Ebulk3 - 2.0*effg)/3.0;
        double effhard = 3.0*uo * hard/(3.0*uo + hard) - (3*effg);

        std::get<Matd6x6>(CMatx_ep).setZero();

        std::get<Matd6x6>(CMatx_ep).topLeftCorner<3, 3>().setConstant(efflam);

        for (int i = 0; i < 3; ++i) {
            std::get<Matd6x6>(CMatx_ep)(i, i) += 2.0*effg;
            std::get<Matd6x6>(CMatx_ep)(i+3, i+3) += effg;
        }

        std::get<Matd6x6>(CMatx_ep) += effhard*(2.0/3.0)*N_tr*(2.0/3.0)*N_tr.transpose();

    }
}

/// @brief Specialization for plane strain.
template <>
inline double IsoHard::Mises2D<PlaneStrain>(const ColVecd3& sig2D, const double& sig_z){

    // Use the components directly from the arguments
    const double sx  = sig2D(0); 
    const double sy  = sig2D(1); 
    const double txy = sig2D(2);  
    const double sz  = sig_z; 

    // Manual squaring for speed in the NR loop
    const double s11_22 = sx - sy;
    const double s22_33 = sy - sz;
    const double s33_11 = sz - sx;

    const double term = 0.5 * (s11_22 * s11_22 + 
                               s22_33 * s22_33 + 
                               s33_11 * s33_11) + 3.0 * (txy * txy);

    return std::sqrt(term);
}

/// @brief Specialization for plane stress.
template <>
inline double IsoHard::Mises2D<PlaneStress>(const ColVecd3& sig2D, const double& sig_z){

    double sx = sig2D(0); 
    double sy = sig2D(1); 
    double txy = sig2D(2);  

    double term = sx * sx - sx * sy + sy * sy + 3 * txy * txy;

    return sqrt(term);
}

/// @brief Specialization 
template <typename AnalysisType, typename HardeningLaw>
void IsoHard::RM2D(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& sig_z, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const double& sig_z_old, const int iStep, const Matd3x3& Ce, Matd3x3& Cep){

    // Elastic strain
    eps_e = eps_e_old + deps;

    // Plastic strain
    eps_p = eps_p_old;

    // Sigma zz
    double sig_z_trial = sig_z_old + ho * (deps(0) + deps(1));

    // Trial stress
    ColVecd3 sig_trial = Ce * eps_e;

    //Mises stress
    double sig_trial_eq = Mises2D<AnalysisType>(sig_trial, sig_z_trial);

    // Initialize plastic multiplier and increment
    double p = eps_eq_old;
    double deqpl = 0;

    // Current yield stress and hardening modulus
    double sYield0, sYield, hard; 
    UHard<HardeningLaw>(p, sYield0, rho, hard); 
    sYield = sYield0;

    // Yield function 
    double f_yield = sig_trial_eq - sYield0*(1.0 + 1.0e-6);

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        sig_z = sig_z_trial;
        sig_h = (sig_trial(0) + sig_trial(1) + sig_z_trial) / 3.0;
        Cep = Ce;

    } else { // --> Plastic step

        // Iteration counter
        int nIter_RM = 0;

        // Deviatoric stress
        double sig_h_tr = (sig_trial(0) + sig_trial(1) + sig_z_trial) / 3.0;

        ColVecd4 sig_trial4;
        sig_trial4 << sig_trial(0),  // Index 0: xx
                      sig_trial(1),  // Index 1: yy
                      sig_z_trial,   // Index 2: zz
                      sig_trial(2);  // Index 3: xy (Shear);

        ColVecd4 sig_trial4_dev = sig_trial4 - sig_h_tr*I4; 

        // Plastic flow direction accounting for zz
        ColVecd4 N_tr4 = (3.0/2.0)*sig_trial4_dev/sig_trial_eq; 

        // Flow direction in 2D 
        ColVecd3 N_tr;
        N_tr << N_tr4(0),      // Index 0: xx
                N_tr4(1),      // Index 1: yy
                N_tr4(3);      // Index 3: xy (Tensor Shear);

        while(abs(f_yield) > tol){

            // Update Iteration counter
            nIter_RM++;

            // Update plastic strain increment 
            deqpl += f_yield / (three_uo + hard);
            p = eps_eq_old + deqpl;
            UHard<HardeningLaw>(p, sYield, rho, hard);
            f_yield = sig_trial_eq - three_uo * deqpl - sYield;

            if(nIter_RM > max_iter){

                std::ostringstream oss;
                oss << "\nReturn mapping did not converge at step: " << iStep << "\n"
                    << "Reached maximum iterations: " << nIter_RM << "\n"
                    << "Final yield function value: " << f_yield << "\n"
                    << "Plastic strain increment: " << deqpl << "\n"
                    << "Plastic strain increment: " << deqpl << "\n"
                    << "Equivalent plastic strain: " << p << "\n";
                
                throw std::runtime_error(oss.str());
            }
        }

        // --------- CONVERGED ---> Update variables

        eps_eq = p;             // Equivalent plastic strain

        ColVecd3 d_eps_p;
        d_eps_p << deqpl * N_tr(0),      // xx
                   deqpl * N_tr(1),      // yy
                   2 * deqpl * N_tr(2);  // xy (Engineering shear)

        eps_p = eps_p_old + d_eps_p;
        eps_e = (eps_e_old + deps) - d_eps_p;
        
        sig = Ce*eps_e;  // Stress tensor

        // 1. Calculate the elastic z-increment (correctly using N_tr4 for zz flow)
        double d_eps_z_e = - (deqpl * N_tr4(2));

        // 2. Calculate the stress INCREMENT (d_sig_z) 
        // Use deps(0) and deps(1) because they are the total increments for this step
        double d_sig_z = ho * (deps(0) + deps(1) + d_eps_z_e) + 2.0 * uo * d_eps_z_e;

        // 3. Update the TOTAL stress by adding the increment to the old state
        sig_z = sig_z_old + d_sig_z;

        sig_eq = Mises2D<AnalysisType>(sig, sig_z);  // Von Mises stress
        sig_h = (sig(0) + sig(1) + sig_z) / 3.0;     // Hydostatic stress

        // Tangent stiffness matrix
        
        RowVecd3 Ce_N = Ce*N_tr;  // Note that we use the Tensor shear in N_tr
        double N_Ce_N = N_tr.dot(Ce_N);
        double denom = (2.0 / 3.0) * hard + N_Ce_N;
        Cep = (Ce - (Ce_N.transpose() * Ce_N)/denom);

    }
}

/// @brief Specialization 
template <typename AnalysisType, typename HardeningLaw>
void IsoHard::RM2DPFF(ColVecd3& deps, ColVecd3& sig, ColVecd3& eps_e, ColVecd3& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& sig_z, double& rho, const ColVecd3& eps_e_old, const ColVecd3& eps_p_old, const double& eps_eq_old, const double& sig_z_old, const int iStep, const double gPhi_d, const double& wp_old, double& wp, const Matd3x3& Ce, Matd3x3& Cep){

    // Elastic strain
    eps_e = eps_e_old + deps;

    // Plastic strain
    eps_p = eps_p_old;
    wp = wp_old;

    // Sigma zz trial (damaged!)
    double sig_z_trial = (sig_z_old + ho * (deps(0) + deps(1))) * gPhi_d;

    // Trial stress (damaged!)
    ColVecd3 sig_trial = Ce * eps_e * gPhi_d;

    // Mises stress (damaged!)
    double sig_trial_eq = Mises2D<AnalysisType>(sig_trial, sig_z_trial);

    // Initialize plastic multiplier and increment
    double p = eps_eq_old;
    double deqpl = 0;

    // Current yield stress and hardening modulus
    double sYield0, sYield, hard; 
    UHard<HardeningLaw>(p, sYield0, rho, hard); 
    sYield = sYield0;

    // Yield function 
    double f_yield = sig_trial_eq - sYield0*(1.0 + 1.0e-6);

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        sig_z = sig_z_old + ho * (deps(0) + deps(1));  // We store the undamaged value 
        sig_h = (sig_trial(0) + sig_trial(1) + sig_z_trial) / 3.0;
        Cep = Ce*gPhi_d;

    } else { // --> Plastic step

        // Iteration counter
        int nIter_RM = 0;

        // Deviatoric stress

        // Deviatoric stress
        double sig_h_tr = (sig_trial(0) + sig_trial(1) + sig_z_trial) / 3.0;

        ColVecd4 sig_trial4;
        sig_trial4 << sig_trial(0),  // Index 0: xx
                      sig_trial(1),  // Index 1: yy
                      sig_z_trial,   // Index 2: zz
                      sig_trial(2);  // Index 3: xy (Shear);

        ColVecd4 sig_trial4_dev = sig_trial4 - sig_h_tr*I4; 

        // Plastic flow direction accounting for zz
        ColVecd4 N_tr4 = (3.0/2.0)*sig_trial4_dev/sig_trial_eq; 

        // Flow direction in 2D  
        ColVecd3 N_tr;
        N_tr << N_tr4(0),    // Index 0: xx
                N_tr4(1),    // Index 1: yy
                N_tr4(3);    // Index 3: xy (Tensor Shear);

        while(abs(f_yield) > tol_PFF){

            // Update Iteration counter
            nIter_RM++;

            // Update plastic strain increment 
            deqpl += f_yield/(three_uo*gPhi_d + hard);
            p = eps_eq_old + deqpl;
            UHard<HardeningLaw>(p, sYield, rho, hard);
            f_yield = sig_trial_eq - gPhi_d*three_uo*deqpl - sYield;

            if(nIter_RM > max_iter){

                std::ostringstream oss;
                oss << "\nReturn mapping did not converge at step: " << iStep << "\n"
                    << "Reached maximum iterations: " << nIter_RM << "\n"
                    << "Final yield function value: " << f_yield << "\n"
                    << "Plastic strain increment: " << deqpl << "\n"
                    << "Damage: " << gPhi_d << "\n"
                    << "Equivalent plastic strain: " << p << "\n";
                
                throw std::runtime_error(oss.str());
            }
        }

        // --------- CONVERGED ---> Update variables

        eps_eq = p;             // Equivalent plastic strain

        ColVecd3 d_eps_p;
        d_eps_p << deqpl * N_tr(0),      // xx
                   deqpl * N_tr(1),      // yy
                   2 * deqpl * N_tr(2);  // xy (Engineering shear)

        eps_p = eps_p_old + d_eps_p;
        eps_e = (eps_e_old + deps) - d_eps_p;
        
        // Undamaged values
        sig = Ce*eps_e;  // Stress tensor

        double d_eps_z_e = - (deqpl * N_tr4(2)); // sig_z
        double d_sig_z = ho * (deps(0) + deps(1) + d_eps_z_e) + 2.0 * uo * d_eps_z_e;
        sig_z = sig_z_old + d_sig_z;  // We store undamaged sig_z

        // Plastic work density using undamaged stress tensor
        // 1. Calculate the plastic strain increment in Z
        // (Since total eps_z = 0 in Plane Strain, deps_p_z = -deps_e_z)
        double d_eps_p_z = - d_eps_z_e; 

        // 2. Calculate the plastic work increment
        double d_wp = sig(0) * d_eps_p(0) +    // sigma_xx * deps_p_xx
                      sig(1) * d_eps_p(1) +    // sigma_yy * deps_p_yy
                      sig_z  * d_eps_p_z +     // sigma_zz * deps_p_zz (The missing piece!)
                      sig(2) * d_eps_p(2);     // tau_xy   * dgamma_p_xy (Engineering shear)

        // 3. Update total plastic work density
        wp = wp_old + d_wp;

        // Apply damage to stresses for the final state
        sig = sig*gPhi_d;   
        sig_eq = Mises2D<AnalysisType>(sig, sig_z*gPhi_d);  // Von Mises stress
        sig_h = (sig(0) + sig(1) + sig_z*gPhi_d) / 3.0;     // Hydostatic stress

        // Tangent stiffness matrix

        RowVecd3 Ce_N = Ce*N_tr;
        double N_Ce_N = N_tr.dot(Ce_N);
        double denom = (2.0 / 3.0) * hard + N_Ce_N;
        Cep = (Ce - (Ce_N.transpose() * Ce_N)/denom)*gPhi_d;

    }
}

inline double IsoHard::MisesAxi(const ColVecd4& sig) {
    const double s_rr = sig(0); 
    const double s_zz = sig(1); 
    const double s_hoop = sig(2); 
    const double t_rz = sig(3); 

    const double d1 = s_rr - s_zz;
    const double d2 = s_zz - s_hoop;
    const double d3 = s_hoop - s_rr;

    // Manual multiplication is typically faster than std::pow for squares
    double term = 0.5 * (d1 * d1 + d2 * d2 + d3 * d3) + 3.0 * (t_rz * t_rz);

    return std::sqrt(term);
}

/// @brief Specialization 
template <typename HardeningLaw>
void IsoHard::RMAxi(ColVecd4& deps, ColVecd4& sig, ColVecd4& eps_e, ColVecd4& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd4& eps_e_old, const ColVecd4& eps_p_old, const double& eps_eq_old, const int iStep, const Matd4x4& Ce, Matd4x4& Cep){

    // Elastic strain
    eps_e = eps_e_old + deps;

    // Plastic strain
    eps_p = eps_p_old;

    // Trial stress
    ColVecd4 sig_trial = Ce*eps_e;

    //Mises stress
    double sig_trial_eq = MisesAxi(sig_trial);

    // Initialize plastic multiplier and increment
    double p = eps_eq_old;
    double deqpl = 0;

    // Current yield stress and hardening modulus
    double sYield0, sYield, hard; 
    UHard<HardeningLaw>(p, sYield0, rho, hard); 
    sYield = sYield0;

    // Yield function 
    double f_yield = sig_trial_eq - sYield0*(1.0 + 1.0e-6);

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        sig_h = sig_trial.dot(I4) / 3.0;
        Cep = Ce; // Zero-overhead copy;

    } else { // --> Plastic step

        // Iteration counter
        int nIter_RM = 0;

        // Deviatoric stress

        ColVecd4 sig_trial_dev = sig_trial - (sig_trial.dot(I4) / 3.0) * I4; 
        ColVecd4 N_tr = (3.0/2.0)*sig_trial_dev/sig_trial_eq; // Plastic flow direction

        while(abs(f_yield) > tol){

            // Update Iteration counter
            nIter_RM++;

            // Update plastic strain increment, Use the pre-calculated constant
            deqpl += f_yield / (three_uo + hard);
            p = eps_eq_old + deqpl;
            UHard<HardeningLaw>(p, sYield, rho, hard);
            f_yield = sig_trial_eq - three_uo * deqpl - sYield;

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

        ColVecd4 d_eps_p;
        d_eps_p << deqpl * N_tr(0),      // rr
                   deqpl * N_tr(1),      // zz
                   deqpl * N_tr(2),      // θθ
                   2 * deqpl * N_tr(3);  // rz (Engineering shear)

        eps_p = eps_p_old + d_eps_p;
        eps_e = (eps_e_old + deps) - d_eps_p;
        
        sig = Ce*eps_e;  // Stress tensor
        sig_eq = MisesAxi(sig);;  // Von Mises stress
        sig_h = sig.dot(I4) / 3.0; // Hydostatic stress

        // Tangent stiffness matrix

        RowVecd4 Ce_N = Ce*N_tr;
        double N_Ce_N = N_tr.dot(Ce_N);
        double denom = (2.0 / 3.0) * hard + N_Ce_N;
        Cep.noalias() = Ce - (Ce_N.transpose() * Ce_N)/denom;

    }
}

/// @brief Specialization 
template <typename HardeningLaw>
void IsoHard::RMAxiPFF(ColVecd4& deps, ColVecd4& sig, ColVecd4& eps_e, ColVecd4& eps_p, double& eps_eq, double& sig_eq, double& sig_h, double& rho, const ColVecd4& eps_e_old, const ColVecd4& eps_p_old, const double& eps_eq_old, const int iStep, const double gPhi_d, const double& wp_old, double& wp, const Matd4x4& Ce, Matd4x4& Cep){

    // Elastic strain
    eps_e = eps_e_old + deps;

    // Plastic strain
    eps_p = eps_p_old;

    // Trial stress
    ColVecd4 sig_trial = Ce*eps_e * gPhi_d;

    //Mises stress
    double sig_trial_eq = MisesAxi(sig_trial);

    // Initialize plastic multiplier and increment
    double p = eps_eq_old;
    double deqpl = 0;

    // Current yield stress and hardening modulus
    double sYield0, sYield, hard; 
    UHard<HardeningLaw>(p, sYield0, rho, hard); 
    sYield = sYield0;

    // Yield function 
    double f_yield = sig_trial_eq - sYield0*(1.0 + 1.0e-6);

    // Check yielding 
    if (f_yield <= 0){  // --> Elastic step

        // Update 
        sig = sig_trial;
        sig_eq = sig_trial_eq;
        sig_h = sig_trial.dot(I4) / 3.0;
        Cep = Ce * gPhi_d; // Zero-overhead copy;

    } else { // --> Plastic step

        // Iteration counter
        int nIter_RM = 0;

        // Deviatoric stress

        ColVecd4 sig_trial_dev = sig_trial - (sig_trial.dot(I4) / 3.0) * I4; 
        ColVecd4 N_tr = (3.0/2.0)*sig_trial_dev/sig_trial_eq; // Plastic flow direction

        while(abs(f_yield) > tol_PFF){

            // Update Iteration counter
            nIter_RM++;

            // Update plastic strain increment, Use the pre-calculated constant
            deqpl += f_yield / (three_uo*gPhi_d + hard);
            p = eps_eq_old + deqpl;
            UHard<HardeningLaw>(p, sYield, rho, hard);
            f_yield = sig_trial_eq - three_uo * gPhi_d * deqpl - sYield;

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

        ColVecd4 d_eps_p;
        d_eps_p << deqpl * N_tr(0),      // rr
                   deqpl * N_tr(1),      // zz
                   deqpl * N_tr(2),      // θθ
                   2 * deqpl * N_tr(3);  // rz (Engineering shear)

        eps_p = eps_p_old + d_eps_p;
        eps_e = (eps_e_old + deps) - d_eps_p;

        sig = Ce*eps_e;    // Undamaged stress tensor

        // sig here is currently undamaged (sig = Ce * eps_e)
        double d_wp = sig(0) * d_eps_p(0) + // sigma_rr * deps_p_rr
                    sig(1) * d_eps_p(1) + // sigma_zz * deps_p_zz
                    sig(2) * d_eps_p(2) + // sigma_hoop * deps_p_hoop
                    sig(3) * d_eps_p(3);   // tau_rz * dgamma_p_rz (Engineering shear)

        wp = wp_old + d_wp;
        
        sig = sig * gPhi_d;    // Damaged stress tensor
        sig_eq = MisesAxi(sig);    // Von Mises stress
        sig_h = sig.dot(I4) / 3.0; // Hydostatic stress

        // Tangent stiffness matrix

        RowVecd4 Ce_N = Ce*N_tr;
        double N_Ce_N = N_tr.dot(Ce_N);
        double denom = (2.0 / 3.0) * hard + N_Ce_N;
        Cep.noalias() = Ce - (Ce_N.transpose() * Ce_N)/denom;
        Cep = Cep*gPhi_d;

    }
}