/**
 * @file BaseElemPFF.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for phase-field fracture elements. 
 * @date 2025-09-21
 * 
 * @copyright Copyright (C) 2025 Abdelrahman Hussein
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  
 */

#ifndef BASEELEMPFF_H
#define BASEELEMPFF_H

#include "H5IO.h"
#include "FiniteElements/BaseElements.h"
#include "Materials/PFF/BasePFF.h"


class BaseElemPFF: public BaseElements{

public:

BaseElemPFF(int nElDim, int nElNodes, int nElGauss, int nElPhiDofs, Logger& logger)
    : BaseElements(nElDim, nElNodes, nElGauss, logger), nElPhiDofs(nElPhiDofs) {};

/**
 * @brief Get dimensions of the element. 
 * 
 * @return int 
 */
int get_nDim() const { return nElDim; };

/**
 * @brief Get the number of element Gauss points. 
 * 
 * @return int 
 */
int get_nGauss() const { return nElGauss; };

/**
 * @brief Get the number of element phi DOFs. Same as number of nodes. 
 * 
 * @return int 
 */
int get_nElPhiDofs() const { return nElPhiDofs; };

/**
 * @brief Return a const reference to the `elemPhiDof`. 
 * 
 * @return const vector<vector<int>>& 
 */
const vector<vector<int>>& get_elemPhiDof() const { return elemPhiDof; };

/**
 * @brief Reads the data `nElements`, `nElementSets` and `elemNodeConn` from hdf5 file.
 * 
 * @param H5File_in  
 */
void ReadElementsData(H5IO &H5File_in, H5IO &H5File_mesh, int iSet);

/**
 * @brief Converts a 2D vector in Voigt notation to tensor notation.
 * 
 * @param v 
 * @return Eigen::Matrix2d 
 */
inline Eigen::Matrix2d voigtToTensor2D(const ColVecd3& v) {

    Eigen::Matrix2d t;
    t << v(0), 0.5 * v(2),
         0.5 * v(2), v(1);
    return t;

}

/**
 * @brief Converts 2D tensor to 2D Voigt notation.
 * 
 * @param t 
 * @return ColVecd3 
 */
inline ColVecd3 tensorToVoigt2D(const Eigen::Matrix2d& t) {

    ColVecd3 v;
    v << t(0,0), t(1,1), 2.0 * t(0,1);
    return v;
    
}

/**
 * @brief Set the const wc object
 * 
 * @param const_wc 
 */
inline void set_const_wc(double const_wc){

    for (auto& row : elem_wc)
    std::fill(row.begin(), row.end(), const_wc);

}

/**
 * @brief Set a decaying wc as a function of concentration. 
 * 
 * @param elCon 
 * @param const_wc 
 * @param wc_min 
 * @param beta 
 */
inline void set_decay_wc(const std::vector<std::vector<double>>& elCon, 
                         double wc_0, double wc_min, double beta, double c_crit) {

    for (int iElem = 0; iElem < nElements; iElem++) {
        for (int iGauss = 0; iGauss < nElGauss; iGauss++) {
            
            // local concentration at the integration point
            double c_hat = accessVec(elCon, iElem, iGauss)/c_crit;
            
            // Formula: w_c(c) = w_min + (w_0 - w_min) * exp(-beta * c_hat)
            accessVec(elem_wc, iElem, iGauss) = wc_min + (wc_0 - wc_min) * std::exp(-beta * c_hat);
        }
    }
}

/**
 * @brief Calculates the elastic crack driving force using `psi_plus` and updates the history parameter `elemH`.
 * 
 */
inline void CalcDrivForcE() {
    for (size_t iElem = 0; iElem < elemH.size(); ++iElem) {
        for (size_t iGauss = 0; iGauss < elemH[iElem].size(); ++iGauss) {
            double drivForce = accessVec(psi_plus, iElem, iGauss) / accessVec(elem_wc, iElem, iGauss);
            accessVec(elemH, iElem, iGauss) = std::max(drivForce, accessVec(elemH, iElem, iGauss));
        }
    }
}

 /**
  * @brief Calculates elastoplastic crack driving force with a threshold based on `psi_plus` and `el_wp`.
  *       For details, see Miehe et al. CMAME (2016) PI and PII.
  * 
  * @param el_wp_ptr Plastic word density.
  * @param eta Contribution percentage H = (ηψ⁺ + (1-η)wₚ) / w꜀.
  */
inline void CalcDrivForcEP(const std::vector<std::vector<double>>* el_wp_ptr, const double eta) {
    for (size_t iElem = 0; iElem < elemH.size(); ++iElem) {
        for (size_t iGauss = 0; iGauss < elemH[iElem].size(); ++iGauss) {

            double psiElastic = accessVec(psi_plus, iElem, iGauss);
            double wpPlastic  = accessVec(*el_wp_ptr, iElem, iGauss);
            double wcLocal    = accessVec(elem_wc, iElem, iGauss);

            double drivForce = (psiElastic*eta + wpPlastic*(1.0-eta)) / wcLocal;

            accessVec(elemH, iElem, iGauss) = std::max(drivForce, accessVec(elemH, iElem, iGauss));
        }
    }
}

/**
 * @brief Calculates elastoplastic crack driving force with a threshold based on `psi_plus` and `el_wp`.
 *        For details, see Miehe et al. CMAME (2016) PI and PII.
 * 
 * @param el_wp_ptr 
 * @param zeta Parameter the controls the post initiation behavoir. Default = 1. 
 * @param eta Contribution percentage H = (ηψ⁺ + (1-η)wₚ) / w꜀.
 */
inline void CalcDrivForcEP_TH(const std::vector<std::vector<double>>* el_wp_ptr, const double zeta, const double eta) {
    for (size_t iElem = 0; iElem < elemH.size(); ++iElem) {
        for (size_t iGauss = 0; iGauss < elemH[iElem].size(); ++iGauss) {

            double psiElastic = accessVec(psi_plus, iElem, iGauss);
            double wpPlastic  = accessVec(*el_wp_ptr, iElem, iGauss);
            double wcLocal    = accessVec(elem_wc, iElem, iGauss);

            // Raw driving force
            double rawForce = (psiElastic*eta + wpPlastic*(1-eta)) / wcLocal;

            // Apply thresholding (subtract 1, clamp to zero)
            double drivForce = zeta*std::max(rawForce - 1.0, 0.0);

            accessVec(elemH, iElem, iGauss) = std::max(drivForce, accessVec(elemH, iElem, iGauss));
        }
    }
}

/**
 * @brief Calculates a state-dependent ductile driving force using a triaxiality-gated plastic work 
 *        contribution and a damage-scaled energy transition.
 * 
 * @param el_wp_ptr Plastic work density pointer.
 * @param elTriax_ptr Triaxiality pointer.
 * @param zeta Parameter the controls the post initiation behavoir. Default = 1. 
 * @param eta Contribution percentage.
 */
inline void CalcDrivForcHybridDuctile_TH(const std::vector<std::vector<double>>* el_wp_ptr, const std::vector<std::vector<double>>* elTriax_ptr, const double zeta, const double eta) {
    for (size_t iElem = 0; iElem < elemH.size(); ++iElem) {
        for (size_t iGauss = 0; iGauss < elemH[iElem].size(); ++iGauss) {

            double psiElastic = accessVec(psi_plus, iElem, iGauss);
            double wpPlastic  = accessVec(*el_wp_ptr, iElem, iGauss);
            double triax  = accessVec(*elTriax_ptr, iElem, iGauss);
            double wcLocal    = accessVec(elem_wc, iElem, iGauss);

            double triaxiality_gate = std::tanh(4.0 * std::max(0.0, triax)); 

            // Combined Plastic Driving Force Contribution
            double wp_contribution = triaxiality_gate * (1.0 - eta);

            double rawForce = (psiElastic * eta + wpPlastic * wp_contribution) / wcLocal;
            accessVec(elemH, iElem, iGauss) = std::max(zeta * std::max(rawForce - 1.0, 0.0), accessVec(elemH, iElem, iGauss));
        }
    }
}

/**
 * @brief Get a constant reference to `el_gPhi_d`.
 * 
 * @return const std::vector<std::vector<double>>& 
 */
const std::vector<std::vector<double>>& getEl_gPhi_d() const ;

/**
 * @brief Get a constant reference to `elPhi_ptr`.
 * 
 * @return const std::vector<std::vector<double>>& 
 */
const std::vector<std::vector<double>>& getElphi() const ;

/**
 * @brief Performs spectral decomposition of the elastic strain tensor and calculates the positive
 * and negative parts of the elastic strain energy density.
 * 
 * @param lam Lame constant lambda.
 * @param Gmod Shear modulus. 
 * @param elStrain_e 
 */
virtual void CalcPsiSpectral(double lam, double Gmod, const T_elStres& elStrain_e) = 0;

/**
 * @brief Calculates the element stiffness matrix.
 * 
 */
virtual void CalcElemStiffMatx() = 0;

/**
 * @brief Calculates the RHS vector.
 * 
 * @param Fint Buffer for vecFp.
 */
virtual void CalcFH(double* FH) = 0;

/**
 * @brief Maps int-point values to nodes for output.
 * 
 * @param nodH Crack driving force.
 * @param nodPsi_plus Positive part of the strain energy density.
 * @param nod_wp plastic work density. 
 * @param nodCount 
 */
virtual void CalcNodVals(vector<double>& nodH, vector<double>& nodPsi_plus, vector<double>& nod_wp, const std::vector<std::vector<double>>* el_wp_ptr, vector<double>& nodCount) = 0;

/**
 * @brief Calculated the integration point phi. 
 * 
 * @param globalBuffer 
 */
virtual void Calc_gPhi_d(const double* globalBuffer) = 0;

protected:      

/// @brief Number of element phi dofs.
const int nElPhiDofs;   

/// @brief PFF internal length scale.
double const_ell = 0.0;

/// @brief Element concentration (temperature) dofs. In this case, it is identical to `elemNodeConn`.
vector<vector<int>> elemPhiDof;  

/// @brief Int-pt phi [nElGauss]. 
vector<vector<double>> elPhi;

/// @brief Int-pt gPhi_d [nElGauss]. 
vector<vector<double>> el_gPhi_d;   

/// @brief Int-pt psi_plus [nElGauss]. 
vector<vector<double>> psi_plus; 

/// @brief Int-pt psi_minus [nElGauss]. 
vector<vector<double>> psi_minus; 

/// @brief Int-pt driving force [nElGauss]. 
vector<vector<double>> elemH;  

/// @brief Int-pt critical energy density [nElGauss]. 
vector<vector<double>> elem_wc; 

/// @brief Pointer of int-pt phi [nElGauss]. 
vector<vector<double>>* elPhi_ptr;

/// @brief Int-pt gPhi_d [nElGauss]. 
vector<vector<double>>* el_gPhi_d_ptr;

};
#endif