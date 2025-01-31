#include<iostream>
#include<algorithm>

#include "FiniteElements/Trapping/Quad4TH.h"
#include "Materials/Trapping/TrapGB.h"

#ifndef DEBUG
#define at(x) operator[](x)
#endif

/*
 Unifying indices:
 i -> number of nodes.
 j -> spatial dimensions.
 k -> number of stresses.
 l -> total displacement dofs.
*/

Quad4TH::Quad4TH(H5IO &H5File_in, Nodes &Nodes, int iSet, Logger& logger)
    : BaseElemTrap(2, 4, 4, 4, logger){ // nElDim, nElNodes, nElGauss, nElConDofs 

    InitShapeFunc();
    ReadElementsData(H5File_in, iSet);
    InitializeElements(Nodes, H5File_in);
}   

Quad4TH::~Quad4TH(){

    // Exit message
    cout << "Quad4TH elements exited correctly" << "\n";
}

void Quad4TH::InitShapeFunc(){

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

RowVecd4 Quad4TH:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd4 shape;

    shape(0) = (1.0-eta-xi+xi*eta)*0.25;
    shape(1) = (1.0-eta+xi-xi*eta)*0.25;
    shape(2) = (1.0+eta+xi+xi*eta)*0.25;
    shape(3) = (1.0+eta-xi-xi*eta)*0.25;

    return shape;
}

Matd2x4 Quad4TH::CalcShapeFuncDeriv(double xi, double eta){

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

void Quad4TH::InitializeElements(Nodes &Nodes, H5IO &H5File_in){

    // Initialize the storages for int-pt flux and phi
    elFlux.resize(nElements);
    elPhi.resize(nElements);

    // Read nodal values of phi
    nodPhi.resize(nNodes);
    string dsetName = "Phi";
    H5File_in.ReadField1D(dsetName, nodPhi);      

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd4x2 dummyElNodCoord; // For node coordinates.
    ColVecd4 dummyElNodPhi;  // For element nodal values of phi.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    // Loop through elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elFlux.at(iElem).resize(nElGauss);
        elPhi.at(iElem).resize(nElGauss);
        gaussPtCart.at(iElem).resize(nElGauss);
        BMat.at(iElem).resize(nElGauss);

        // Loop through nodes to get coordinates.
        for(int iNod=0; iNod<nElNodes; iNod++){

            dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(0);
            dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(1);

            dummyElNodPhi[iNod] = nodPhi.at(elemNodeConn.at(iElem).at(iNod));
        }

        elemNodCoord.at(iElem) = dummyElNodCoord;

        // Loop through integration points.
        for(int iGauss=0; iGauss<nElGauss; iGauss++){
        
            // Cart coord of iGauss point.
            dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
            // Derivatives and int-pt volume
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss));
            // Int-pt phi
            elPhi.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNodPhi;
        }
        gaussPtCart.at(iElem) = dummyElemGauss;
        intPtVol.at(iElem) = dummyIntVol;
    }
}

RowVecd2 Quad4TH::getGaussCart(RowVecd4& sFunc, Matd4x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Quad4TH::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv){

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

void Quad4TH::getInPtCoords(T_nodStres& glIntPtCoords){
    
}

void Quad4TH::CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi){

    // Int-pt gradients of phi.
    vector<vector<ColVecd2>> elGradPhi(nElements);// For element nodal values of phi.
    ColVecd4 dummyElNodPhi;        // For element nodal values of phi.
    ColVecd4 dummyElNodLapPhi;  // For element nodal values of phi laplacian.
    int iNode;  // counter for the number of nodes.

    for(int iElem=0; iElem<nElements; iElem++){

        elGradPhi.at(iElem).resize(nElGauss);

        // Loop through element nodes to get nodal values.
        for(int iNod=0; iNod<nElNodes; iNod++){
            dummyElNodPhi[iNod] = nodPhi.at(elemNodeConn.at(iElem).at(iNod));
        }

        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            const Matd2x4& dummyBMat = BMat.at(iElem).at(iGaus); // derivative matrix for the given gauss point.
            elGradPhi.at(iElem).at(iGaus) = dummyBMat*dummyElNodPhi;
            dummyElNodLapPhi = dummyBMat.transpose()*dummyBMat*dummyElNodPhi;

            iNode = 0;
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd2>>(nodGrad).at(*iNod2) += elGradPhi.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                nodCount.at(*iNod2) += shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                nodLapPhi[*iNod2]   += dummyElNodLapPhi[iNode]; 
                iNode += 1;
            }
        }
    }
}

void Quad4TH::CalcElemStiffMatx(BaseTrapping* mat, const double T){

    // Matd2x2 KMat; 
    // Matd2x2 TMat; 

    // elStiffMatx.resize(nElements);
    // elCapMatx.resize(nElements);
    // elMKTMatx.resize(nElements);
    // vector<Matd4x4> elKDMatx(nElements); 
    // vector<Matd4x4> elKTMatx(nElements); 

    // double dummydVol;        // dummy for int-pt volume.
    // ColVecd4 dummyElNodPhi;  // For element nodal values of phi.
    // double phi;              // dummy for int-pt phi.

    // // Loop through all elements.
    // for(int iElem=0; iElem<nElements; iElem++){

    //     // MUST BE POPULATED WITH ZEROS    
    //     elStiffMatx.at(iElem).setZero();
    //     elKDMatx.at(iElem).setZero();  
    //     elKTMatx.at(iElem).setZero();  
    //     elCapMatx.at(iElem).setZero(); 

    //     // Loop through element nodes to get nodal values.
    //     for(int iNod=0; iNod<nElNodes; iNod++){
    //         dummyElNodPhi[iNod] = nodPhi.at(elemNodeConn.at(iElem).at(iNod));
    //     }              

    //     // Integration over all Gauss points.
    //     for (int iGauss=0; iGauss<nElGauss; iGauss++){

    //         phi = elPhi.at(iElem).at(iGauss);  // Phi of the current int-pt

    //         KMat = std::get<Matd2x2>(dynamic_cast<PhaseTrap*>(mat)->CalcKMatx((1-phi),  phi, T));
    //         TMat = std::get<Matd2x2>(dynamic_cast<PhaseTrap*>(mat)->CalcTMatx((1-phi),  phi, T));

    //         const Matd2x4& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
    //         const RowVecd4& dummyShFunc = shapeFunc.at(iGauss);
    //         dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current int-pt 
    //         // [B_ji]^T k_jj B_ji
    //         elKDMatx.at(iElem).noalias() += dummyBMat.transpose()*KMat*dummyBMat*dummydVol;
    //         // [B_ji]^T k_jj B_ji
    //         // elKTMatx.at(iElem).noalias() += dummyBMat.transpose()*TMat*dummyBMat*dummyElNodPhi*dummyShFunc*dummydVol; 
    //         elKTMatx.at(iElem).noalias() += dummyBMat.transpose()*TMat*dummyBMat*phi*dummydVol; 

    //         // [N_i]^T s N_i
    //         elCapMatx.at(iElem).noalias() += (dummyShFunc.transpose()*dummyShFunc)*dummydVol;
    //     }

    //     // cout << elKTMatx.at(iElem) << "\n\n";
    //     elStiffMatx.at(iElem) = dt*(elKDMatx.at(iElem) - elKTMatx.at(iElem)) + elCapMatx.at(iElem);
    //     elMKTMatx.at(iElem) = elCapMatx.at(iElem);
    //     // elMKTMatx.at(iElem) = elCapMatx.at(iElem);
    // }

    // // // TODO: For debug!
    // // for (auto& iStifMat : elKDMatx)
    // //     cout << iStifMat << "\n\n";
    // // cout << elKDMatx.at(0) << "\n";

    // // Pointer to the vector, not the vector itself.
    // elStiffMatxVariant = &elStiffMatx;
    // // elCapMatxVariant = &elCapMatx;
    // elMKTMatxVariant = &elMKTMatx;

}

void Quad4TH::CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T){

    // ColVecd4 dummyCon; // for element nodal concentration.
    // int iNode;  // counter for the number of nodes.

    // // Integration point values.
    // for(int iElem=0; iElem<nElements; iElem++){

    //     // Get element nodal concentration from the solution vector. 
    //     for(int iDof=0; iDof<nElConDofs; iDof++){
    //         dummyCon(iDof) = globalBuffer[elemConDof.at(iElem).at(iDof)];
    //     }

    //     // Gauss points
    //     for(int iGaus=0; iGaus<nElGauss; iGaus++){

    //         // Int pt flux
    //         elFlux.at(iElem).at(iGaus) = -(std::get<Matd2x2>(KMatx)*BMat.at(iElem).at(iGaus))*dummyCon;

    //         // Nodal values
    //         iNode = 0;
    //         for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

    //             std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += elFlux.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
    //             nodCount.at(*iNod2) += shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
    //             iNode += 1;
    //         }
    //     }
    // }

    // // // TODO: For debug!
    // // cout << elFlux.at(0).at(0) << "\n\n";
}
