#include<iostream>
#include<algorithm>

#include "FiniteElements/Trapping/Quad4THAxi.h"

Quad4THAxi::Quad4THAxi(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, Logger& logger, H5IO* H5File_rve)
    : Quad4TH(H5File_in, H5File_mesh, Nodes, iSet, logger, H5File_rve) {}

void Quad4THAxi::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const RowVecd4& shFunc, 
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