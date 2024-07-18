#include<iostream>
#include<algorithm>

#include "FiniteElements/Trapping/Tri3TH.h"
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

Tri3TH::Tri3TH(H5IO &H5File_in, Nodes &Nodes, int iSet)
    : BaseElemTrap(2, 3, 3, 3){ // nElDim, nElNodes, nElConDofs, nElGauss 

    InitShapeFunc();
    ReadElementsData(H5File_in, iSet);
    InitializeElements(Nodes, H5File_in);
}   

Tri3TH::~Tri3TH(){

    // Exit message
    cout << "Tri3TH elements exited correctly" << "\n";
}

void Tri3TH::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    gaussPts.push_back({1.0/2.0, 1.0/2.0});
    gaussPts.push_back({0.0, 1.0/2.0});
    gaussPts.push_back({1.0/2, 0.0});

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nElGauss);
    shapeFuncDeriv.resize(nElGauss);

    for(int i=0; i<nElGauss; i++){
        shapeFunc.at(i) =  CalcShapeFunc(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
        shapeFuncDeriv.at(i) = CalcShapeFuncDeriv(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
    }
}

RowVecd3 Tri3TH:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd3 shape;

    shape(0) = xi;
    shape(1) = eta;
    shape(2) = 1 - xi - eta;

    return shape;
}

Matd2x3 Tri3TH::CalcShapeFuncDeriv(double xi, double eta){

    // dN_ji
    Matd2x3 shapeDeriv;

    shapeDeriv(0,0) = 1.0;
    shapeDeriv(0,1) = 0.0;
    shapeDeriv(0,2) = -1.0;

    shapeDeriv(1,0) = 0.0;
    shapeDeriv(1,1) = 1.0;
    shapeDeriv(1,2) = -1.0;
             
    return shapeDeriv;
}

void Tri3TH::InitializeElements(Nodes &Nodes, H5IO &H5File_in){

    // Initialize the storages
    elFlux.resize(nElements);
    elStiffMatx.resize(nElements);
    elCapMatx.resize(nElements);

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd3x2 dummyElNodCoord;   // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    if (Trapping==1){       // GB

        ColVecd3 dummyElNod_gPhi;  // For element nodal values of phi.

        // Initialize the storages for int-pt gPhi
        el_gPhi.resize(nElements);
        // Read nodal values of gPhi
        nod_gPhi.resize(nNodes);
        string dsetName = "gPhi";
        H5File_in.ReadFieldDoub1D(dsetName, nod_gPhi);

        // Loop through elements.
        for(int iElem=0; iElem<nElements; iElem++){

            elFlux.at(iElem).resize(nElGauss);
            el_gPhi.at(iElem).resize(nElGauss);


            gaussPtCart.at(iElem).resize(nElGauss);
            BMat.at(iElem).resize(nElGauss);

            // Loop through nodes to get coordinates.
            for(int iNod=0; iNod<nElNodes; iNod++){

                dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(0);
                dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(1);

                dummyElNod_gPhi[iNod] = nod_gPhi.at(elemNodeConn.at(iElem).at(iNod));
            }

            elemNodCoord.at(iElem) = dummyElNodCoord;

            // Loop through integration points.
            for(int iGauss=0; iGauss<nElGauss; iGauss++){
            
                // Cart coord of iGauss point.
                dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
                // Derivatives and int-pt volume
                CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss));
                // Int-pt phi
                el_gPhi.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_gPhi;
            }
            gaussPtCart.at(iElem) = dummyElemGauss;
            intPtVol.at(iElem) = dummyIntVol;
        }   

    } else if (Trapping==2) {       // Phase 

        ColVecd3 dummyElNod_martensite; ColVecd3 dummyElNod_gPhiMM;
        ColVecd3 dummyElNod_gPhiff; ColVecd3 dummyElNod_gPhifM;

        // Initialize the storages for int-pt traps
        el_martensite.resize(nElements); el_gPhiff.resize(nElements);
        el_gPhiMM.resize(nElements); el_gPhifM.resize(nElements);

        // Read nodal values of traps
        nod_martensite.resize(nNodes);
        string dsetName = "martensite";
        H5File_in.ReadFieldDoub1D(dsetName, nod_martensite);

        nod_gPhiff.resize(nNodes);
        dsetName = "gPhi_ff";
        H5File_in.ReadFieldDoub1D(dsetName, nod_gPhiff);

        nod_gPhifM.resize(nNodes);
        dsetName = "gPhi_fM";
        H5File_in.ReadFieldDoub1D(dsetName, nod_gPhifM);

        nod_gPhiMM.resize(nNodes);
        dsetName = "gPhi_MM";
        H5File_in.ReadFieldDoub1D(dsetName, nod_gPhiMM);

        // Loop through elements.
        for(int iElem=0; iElem<nElements; iElem++){

            elFlux.at(iElem).resize(nElGauss);

            el_martensite.at(iElem).resize(nElGauss);
            el_gPhiff.at(iElem).resize(nElGauss);
            el_gPhifM.at(iElem).resize(nElGauss);
            el_gPhiMM.at(iElem).resize(nElGauss);

            gaussPtCart.at(iElem).resize(nElGauss);
            BMat.at(iElem).resize(nElGauss);

            // Loop through nodes to get coordinates.
            for(int iNod=0; iNod<nElNodes; iNod++){

                dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(0);
                dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(1);

                dummyElNod_martensite[iNod] = nod_martensite.at(elemNodeConn.at(iElem).at(iNod));
                dummyElNod_gPhiff[iNod] = nod_gPhiff.at(elemNodeConn.at(iElem).at(iNod));
                dummyElNod_gPhifM[iNod] = nod_gPhifM.at(elemNodeConn.at(iElem).at(iNod));
                dummyElNod_gPhiMM[iNod] = nod_gPhiMM.at(elemNodeConn.at(iElem).at(iNod));
            }

            elemNodCoord.at(iElem) = dummyElNodCoord;

            // Loop through integration points.
            for(int iGauss=0; iGauss<nElGauss; iGauss++){
            
                // Cart coord of iGauss point.
                dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
                // Derivatives and int-pt volume
                CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss));
                // Int-pt phi
                el_martensite.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_martensite;
                el_gPhiff.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_gPhiff;
                el_gPhifM.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_gPhifM;
                el_gPhiMM.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_gPhiMM;

            }
            gaussPtCart.at(iElem) = dummyElemGauss;
            intPtVol.at(iElem) = dummyIntVol;
        }   

    }  else {

        cerr << "Error trapping flag: " << Trapping << "\n";
        cerr << "Terminating!\n\n";
        exit(10);
    }
}

RowVecd2 Tri3TH::getGaussCart(RowVecd3& sFunc, Matd3x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Tri3TH::CalcCartDeriv(Matd3x2& elNodCoord, Matd2x3& sFuncDeriv, const double& wt, double& intVol, Matd2x3& cartDeriv){

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

void Tri3TH::CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi){

    // Int-pt gradients of phi.
    vector<vector<ColVecd2>> elGrad_gPhi(nElements);    // For element nodal values of gPhi.
    ColVecd3 dummyElNodPhi;        // For element nodal values of phi.
    ColVecd3 dummyElNodLapPhi;     // For element nodal values of phi laplacian.
    int iNode;                     // counter for the number of nodes.

    for(int iElem=0; iElem<nElements; iElem++){

        elGrad_gPhi.at(iElem).resize(nElGauss);

        // Loop through element nodes to get nodal values.
        for(int iNod=0; iNod<nElNodes; iNod++){
            dummyElNodPhi[iNod] = nod_gPhi.at(elemNodeConn.at(iElem).at(iNod));
        }

        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            const Matd2x3& dummyBMat = BMat.at(iElem).at(iGaus); // derivative matrix for the given gauss point.
            elGrad_gPhi.at(iElem).at(iGaus) = dummyBMat*dummyElNodPhi;
            dummyElNodLapPhi = dummyBMat.transpose()*dummyBMat*dummyElNodPhi;

            iNode = 0;
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd2>>(nodGrad).at(*iNod2) += elGrad_gPhi.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += 1;
                nodLapPhi[*iNod2]   += dummyElNodLapPhi[iNode]; 
                iNode += 1;
            }
        }
    }
}

void Tri3TH::CalcElemStiffMatx(BaseTrapping* mat, const double T){

    Matd2x2 DMat; 
    Matd2x2 TMat; 

    elStiffMatx.resize(nElements);
    elCapMatx.resize(nElements);
    vector<Matd3x3> elKDMatx(nElements); 
    vector<Matd3x3> elKTMatx(nElements); 

    double dummydVol;          // dummy for int-pt volume.
    ColVecd3 dummyElNod_gPhi;  // For element nodal values of phi.
    double gPhi;               // dummy for int-pt phi.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        // MUST BE POPULATED WITH ZEROS    
        elStiffMatx.at(iElem).setZero();
        elKDMatx.at(iElem).setZero();  
        elKTMatx.at(iElem).setZero();  
        elCapMatx.at(iElem).setZero(); 

        // Loop through element nodes to get nodal values.
        for(int iNod=0; iNod<nElNodes; iNod++){
            dummyElNod_gPhi[iNod] = nod_gPhi.at(elemNodeConn.at(iElem).at(iNod));
        }              

        // Integration over all Gauss points.
        for (int iGauss=0; iGauss<nElGauss; iGauss++){

            gPhi = el_gPhi.at(iElem).at(iGauss);  // gPhi of the current int-pt

            DMat = std::get<Matd2x2>(dynamic_cast<TrapGB*>(mat)->CalcDMatx(gPhi, T));
            TMat = std::get<Matd2x2>(dynamic_cast<TrapGB*>(mat)->CalcTMatx(gPhi, T));

            const Matd2x3& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
            const RowVecd3& dummyShFunc = shapeFunc.at(iGauss);
            dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current int-pt 

            // [B_ji]^T k_jj B_ji
            elKDMatx.at(iElem).noalias() += dummyBMat.transpose()*DMat*dummyBMat*dummydVol;
            // [B_ji]^T k_jj B_ji
            elKTMatx.at(iElem).noalias() += dummyBMat.transpose()*TMat*dummyBMat*dummyElNod_gPhi*dummyShFunc*dummydVol; 
            // [N_i]^T N_i
            elCapMatx.at(iElem).noalias() += (dummyShFunc.transpose()*dummyShFunc)*dummydVol;
        }

        elStiffMatx.at(iElem) = dt*elKDMatx.at(iElem) - dt*elKTMatx.at(iElem) + elCapMatx.at(iElem);
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elKDMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
    elCapMatxVariant = &elCapMatx;
}

void Tri3TH::CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, vector<double>& nodCount, const double T){

    Matd2x2 DMat; 
    Matd2x2 TMat; 
    double gPhi;
    double IntPtCon;   // Integration point concnetration
    ColVecd3 dummyCon; // for element nodal concentration.
    ColVecd3 dummyElNod_gPhi; // for element nodal concentration.

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal concentration and gPhi from the solution vector. 
        for(int iDof=0; iDof<nElConDofs; iDof++){
            dummyCon[iDof] = globalBuffer[elemConDof.at(iElem).at(iDof)];
            dummyElNod_gPhi[iDof] = nod_gPhi.at(elemConDof.at(iElem).at(iDof));
        }

        // Gauss points
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            gPhi = el_gPhi.at(iElem).at(iGaus);
            IntPtCon = shapeFunc.at(iGaus)*dummyCon;

            DMat = std::get<Matd2x2>(dynamic_cast<TrapGB*>(mat)->CalcDMatx(gPhi, T));
            TMat = std::get<Matd2x2>(dynamic_cast<TrapGB*>(mat)->CalcTMatx(gPhi, T));

            // Int pt flux
            elFlux.at(iElem).at(iGaus) = - DMat*BMat.at(iElem).at(iGaus)*dummyCon 
                                         + TMat*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_gPhi;

            // elFlux.at(iElem).at(iGaus) = - DMat*BMat.at(iElem).at(iGaus)*dummyCon;

            // elFlux.at(iElem).at(iGaus) = TMat*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_gPhi;

            // Nodal values
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){
                std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += elFlux.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += 1;
            }
        }
    }

    // // TODO: For debug!
    // cout << elFlux.at(0).at(0) << "\n\n";
}
