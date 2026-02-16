#include<iostream>
#include<algorithm>

#include "FiniteElements/PFF/Quad4PFF.h"

/*
 Unifying indices:
 i -> number of nodes.
 j -> spatial dimensions.
 k -> number of stresses.
 l -> total displacement dofs.
*/

Quad4PFF::Quad4PFF(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, double ell, Logger& logger)
    : BaseElemPFF(2, 4, 4, 4, logger){ // nElDim, nElNodes, nElGauss, nElConDofs 

    InitShapeFunc();
    ReadElementsData(H5File_in, H5File_mesh, iSet);
    const_ell = ell;
    InitializeElements(Nodes, H5File_in);
}   

Quad4PFF::~Quad4PFF(){}

void Quad4PFF::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    vector<double> ip = {-0.57735027, 0.57735027};
    vector<double> dummy(nElDim);
    for(std::size_t iGauss=0; iGauss<ip.size(); iGauss++){
        for(std::size_t jGauss=0; jGauss<ip.size(); jGauss++){
            dummy.at(0) = ip.at(iGauss);
            dummy.at(1) = ip.at(jGauss);
            gaussPts.push_back(dummy);
        }
    }

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nElGauss);
    shapeFuncDeriv.resize(nElGauss);

    for(int i=0; i<nElGauss; i++){
        shapeFunc.at(i) =  CalcShapeFunc(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
        shapeFuncDeriv.at(i) = CalcShapeFuncDeriv(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
    }
}

RowVecd4 Quad4PFF:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd4 shape;

    shape(0) = (1.0-eta-xi+xi*eta)*0.25;
    shape(1) = (1.0-eta+xi-xi*eta)*0.25;
    shape(2) = (1.0+eta+xi+xi*eta)*0.25;
    shape(3) = (1.0+eta-xi-xi*eta)*0.25;

    return shape;
}

Matd2x4 Quad4PFF::CalcShapeFuncDeriv(double xi, double eta){

    // dN_ji
    Matd2x4 shapeDeriv;

    shapeDeriv(0,0) = (-1.0+eta)*0.25;
    shapeDeriv(0,1) = (1.0-eta)*0.25;
    shapeDeriv(0,2) = (1.0+eta)*0.25;
    shapeDeriv(0,3) = (-1.0-eta)*0.25;

    shapeDeriv(1,0) = (-1.0+xi)*0.25;
    shapeDeriv(1,1) = (-1.0-xi)*0.25;
    shapeDeriv(1,2) = (1.0+xi)*0.25;
    shapeDeriv(1,3) = (1.0-xi)*0.25;
             
    return shapeDeriv;
}

void Quad4PFF::InitializeElements(Nodes &Nodes, H5IO &H5File_in){

    // Initialize the storages 
    elStiffMatx.resize(nElements);
    elPhi.resize(nElements);
    el_gPhi_d.resize(nElements);
    psi_plus.resize(nElements);
    psi_minus.resize(nElements);
    elemH.resize(nElements);
    elem_wc.resize(nElements);

    // After allocation
    elStiffMatxVariant = &elStiffMatx;
    el_gPhi_d_ptr = &el_gPhi_d;
    elPhi_ptr = &elPhi;

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd4x2 dummyElNodCoord; // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    try{

        // Loop through elements.
        for(int iElem=0; iElem<nElements; iElem++){

            accessVec(elPhi, iElem).resize(nElGauss);
            accessVec(el_gPhi_d, iElem).resize(nElGauss);
            accessVec(psi_plus, iElem).resize(nElGauss);
            accessVec(psi_minus, iElem).resize(nElGauss);
            accessVec(elemH, iElem).resize(nElGauss);
            accessVec(elem_wc, iElem).resize(nElGauss);

            accessVec(gaussPtCart, iElem).resize(nElGauss);
            accessVec(BMat, iElem).resize(nElGauss);

            // Loop through nodes to get coordinates.
            for(int iNod=0; iNod<nElNodes; iNod++){

                dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(0);
                dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(1);
            }

            elemNodCoord.at(iElem) = dummyElNodCoord;

            // Loop through integration points.
            for(int iGauss=0; iGauss<nElGauss; iGauss++){
            
                // Cart coord of iGauss point.
                dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
                // Derivatives and int-pt volume
                CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), shapeFunc.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss));

            }

            gaussPtCart.at(iElem) = dummyElemGauss;
            intPtVol.at(iElem) = dummyIntVol;
        }

    } catch (const std::runtime_error& e) {

        logger.log("\nException caught in Quad4PFF::InitializeElements:\n", "Error");
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

    }
}

RowVecd2 Quad4PFF::getGaussCart(RowVecd4& sFunc, Matd4x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Quad4PFF::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const RowVecd4& shFunc, const double& wt, double& intVol, Matd2x4& cartDeriv){

    // Calculates the jacobian matrix J_jj = dN_ji x_ij
    Matd2x2 jacMat = sFuncDeriv*elNodCoord;

    // Jacobian determinant.
    intVol = jacMat.determinant()*wt;
        
#ifdef DEBUG

    if(intVol<=0){

        cerr << "Error: Negative jacobian determinant resulting in negative area.\n"
             << "Terminating!\n\n";
        exit(10);
    }

#endif

    // Cartesian derivatives for the current integration point dx_ji = [J^-1]_jj dN_ji
    cartDeriv = jacMat.inverse()*sFuncDeriv;
}

void Quad4PFF::CalcPsiSpectral(double lam, double Gmod, const T_elStres& elStrain_e){

    Eigen::Vector2d eigvals;
    Eigen::Matrix2d eigvecs;

    // Positive and negative parts of the strain tensor
    Eigen::Matrix2d E_plus  = Eigen::Matrix2d::Zero();
    Eigen::Matrix2d E_minus = Eigen::Matrix2d::Zero();

    // Access and dereference the pointer
    std::vector<std::vector<ColVecd3>>& elStran = *std::get<std::vector<std::vector<ColVecd3>>*>(elStrain_e);
    
    for(int iElem=0; iElem<nElements; iElem++){
        for (int iGauss=0; iGauss<nElGauss; iGauss++){
            
            // Get ip Voigt strain reference
            const ColVecd3& stran_v =  elStran.at(iElem).at(iGauss);

            // Convert to tensor
            Eigen::Matrix2d stran_t = voigtToTensor2D(stran_v);

            // Calculate the trace
            double tr_eps = stran_t.trace();
            double tr_plus  = std::max(tr_eps, 0.0);
            double tr_minus = std::min(tr_eps, 0.0);
            
            // Spectral decomposition
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(stran_t);

            eigvals  = eigensolver.eigenvalues();  // Principal strains
            eigvecs = eigensolver.eigenvectors();  // Principal directions

            // Positive and negative parts of the strain tensor
            E_plus.setZero(); E_minus.setZero();

            for (int i = 0; i < 2; ++i) {
                double eigval = eigvals(i);
                Eigen::Vector2d eigvec = eigvecs.col(i);

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

void Quad4PFF::CalcElemStiffMatx(){

    try{

            // Loop through all elements.
            for(int iElem=0; iElem<nElements; iElem++){

                // MUST BE POPULATED WITH ZEROS    
                accessVec(elStiffMatx, iElem).setZero();        

                // Integration over all Gauss points.
                for (int iGauss=0; iGauss<nElGauss; iGauss++){

                    const Matd2x4& dummyBMat = accessVec(BMat, iElem, iGauss); // derivative matrix for the given gauss point.
                    const RowVecd4& dummyShFunc = accessVec(shapeFunc, iGauss);
                    const double& dummydVol = accessVec(intPtVol, iElem, iGauss);  // Volume of the current int-pt 

                    accessVec(elStiffMatx, iElem).noalias() += (
                        (accessVec(elemH, iElem, iGauss) + 1)*(dummyShFunc.transpose()*dummyShFunc) +
                        (const_ell*const_ell*dummyBMat.transpose()*dummyBMat)
                    ) * dummydVol;
                }
            } 

    } catch (const std::runtime_error& e) {

        logger.log("\nException caught in Quad4PFF::CalcElemStiffMatx:\n", "", false);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

    }

    // // TODO: For debug!
    // for (auto& iStifMat : elStiffMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n\n";

}

void Quad4PFF::CalcFH(double* FH){

    ColVecd4 dummyFp; // for element nodal RHS.

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGauss=0; iGauss<nElGauss; iGauss++){

            dummyFp = accessVec(shapeFunc, iGauss).transpose()*accessVec(elemH, iElem, iGauss)*accessVec(intPtVol, iElem, iGauss);

            for(int pom=0; pom<nElPhiDofs; pom++){
                FH[accessVec(elemPhiDof, iElem, pom)] += dummyFp(pom);
            }
        }
    }
}

void Quad4PFF::Calc_gPhi_d(const double* globalBuffer){

    ColVecd4 dummyPhi;
    double phi = 0;

    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElPhiDofs; iDof++){
            dummyPhi(iDof) = globalBuffer[accessVec(elemPhiDof, iElem, iDof)];
        }

        for(int iGauss=0; iGauss<nElGauss; iGauss++){

            double rawPhi = accessVec(shapeFunc, iGauss) * dummyPhi;
    
            double clampedPhi = std::max(0.0, std::min(1.0, rawPhi));            
            accessVec(elPhi, iElem, iGauss) = clampedPhi;
            
            double phi = 1.0 - clampedPhi;
            accessVec(el_gPhi_d, iElem, iGauss) = (phi * phi) + 1e-6;
        }
    }
}