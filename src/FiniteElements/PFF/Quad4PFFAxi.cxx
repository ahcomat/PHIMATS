#include<iostream>
#include<algorithm>

#include "FiniteElements/PFF/Quad4PFFAxi.h"

Quad4PFFAxi::Quad4PFFAxi(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, double ell, Logger& logger)
    : Quad4PFF(H5File_in, H5File_mesh, Nodes, iSet, ell, logger) {}

void Quad4PFFAxi::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const RowVecd4& shFunc, 
                            const double& wt, double& intVol, Matd2x4& cartDeriv) {
    // Kinematics from Quad4Axi
    Matd2x2 jacMat = sFuncDeriv * elNodCoord;
    
    // Radius r = sum(N_i * r_i) where col 0 is r
    double radius = shFunc.dot(elNodCoord.col(0));

    // Axisymmetric volume: 2 * PI * r * det(J) * w
    intVol = 2.0 * M_PI * radius * jacMat.determinant() * wt;

    // Standard Cartesian derivatives
    cartDeriv = jacMat.inverse() * sFuncDeriv;
}

void Quad4PFFAxi::CalcPsiSpectral(double lam, double Gmod, const T_elStres& elStrain_e) {

    Eigen::Vector3d eigvals;
    Eigen::Matrix3d eigvecs;

    // Positive and negative parts of the strain tensor (3x3 for Axisymmetry)
    Eigen::Matrix3d E_plus  = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d E_minus = Eigen::Matrix3d::Zero();

    // Axisymmetric strain has 4 components: [err, ezz, ett, erz]
    auto& elStran = *std::get<std::vector<std::vector<ColVecd4>>*>(elStrain_e);
    
    for(int iElem=0; iElem<nElements; iElem++){
        for (int iGauss=0; iGauss<nElGauss; iGauss++){
            const ColVecd4& stran_v = elStran.at(iElem).at(iGauss);

            // Construct 3x3 tensor to include hoop strain (index 2)
            Eigen::Matrix3d stran_t = Eigen::Matrix3d::Zero();
            stran_t(0,0) = stran_v(0); // r
            stran_t(1,1) = stran_v(1); // z
            stran_t(2,2) = stran_v(2); // theta (hoop)
            stran_t(0,1) = stran_t(1,0) = stran_v(3) * 0.5; // rz shear
            
            // Calculate the trace
            double tr_eps = stran_t.trace();
            double tr_plus  = std::max(tr_eps, 0.0);
            double tr_minus = std::min(tr_eps, 0.0);

            // Perform spectral decomposition on 3D tensor
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(stran_t);

            eigvals  = eigensolver.eigenvalues();  // Principal strains
            eigvecs = eigensolver.eigenvectors();  // Principal directions

            // Positive and negative parts of the strain tensor
            E_plus.setZero(); E_minus.setZero();

            for (int i = 0; i < 3; ++i) {
                double eigval = eigvals(i);
                Eigen::Vector3d eigvec = eigvecs.col(i);

                E_plus  += std::max(eigval, 0.0) * eigvec * eigvec.transpose();
                E_minus += std::min(eigval, 0.0) * eigvec * eigvec.transpose();
            }

            psi_plus.at(iElem).at(iGauss)  = 0.5 * lam * tr_plus  * tr_plus
                         + Gmod * (E_plus * E_plus).trace();

            psi_minus.at(iElem).at(iGauss) = 0.5 * lam * tr_minus * tr_minus
                                    + Gmod * (E_minus * E_minus).trace();

        }
    }
}