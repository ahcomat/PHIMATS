#include<iostream>
#include<algorithm>

#include "FiniteElements/Trapping/Quad4TH.h"
#include "Materials/Trapping/TrapGB.h"
#include "Materials/Trapping/TrapPhase.h"
#include "Materials/Trapping/MechTrap.h"

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
    elStiffMatx.resize(nElements);
    elCapMatx.resize(nElements);     
    
    elCon.resize(nElements);
    elCon_ptr = &elCon;

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd4x2 dummyElNodCoord; // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    try{

        if (Trapping=="GBTrapping"){    // GB

            ColVecd4 dummyElNod_gPhi;   // For element nodal values of phi.

            // Initialize the storages for int-pt gPhi
            el_gPhi.resize(nElements);
            // Read nodal values of gPhi
            nod_gPhi.resize(nNodes);
            string dsetName = "gPhi";
            H5File_in.ReadField1D(dsetName, nod_gPhi);

            // Loop through elements.
            for(int iElem=0; iElem<nElements; iElem++){

                accessVec(elFlux, iElem).resize(nElGauss);
                accessVec(elCon, iElem).resize(nElGauss);

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

        } else if (Trapping=="2PhaseTrapping") {       // 2Phase 

            ColVecd4 dummyElNod_phi_j; ColVecd4 dummyElNod_gPhi_jj;
            ColVecd4 dummyElNod_gPhi_ii; ColVecd4 dummyElNod_gPhi_ij;

            // Initialize the storages for int-pt traps
            el_phi_j.resize(nElements); el_gPhi_ii.resize(nElements);
            el_gPhi_jj.resize(nElements); el_gPhi_ij.resize(nElements);

            // Read nodal values of traps
            nod_phi_j.resize(nNodes);
            string dsetName = "phi_j";
            H5File_in.ReadField1D(dsetName, nod_phi_j);

            nod_gPhi_ii.resize(nNodes);
            dsetName = "gPhi_ii";
            H5File_in.ReadField1D(dsetName, nod_gPhi_ii);

            nod_gPhi_ij.resize(nNodes);
            dsetName = "gPhi_ij";
            H5File_in.ReadField1D(dsetName, nod_gPhi_ij);

            nod_gPhi_jj.resize(nNodes);
            dsetName = "gPhi_jj";
            H5File_in.ReadField1D(dsetName, nod_gPhi_jj);

            // Loop through elements.
            for(int iElem=0; iElem<nElements; iElem++){

                accessVec(elFlux, iElem).resize(nElGauss);
                accessVec(elCon, iElem).resize(nElGauss);

                el_phi_j.at(iElem).resize(nElGauss);
                el_gPhi_ii.at(iElem).resize(nElGauss);
                el_gPhi_ij.at(iElem).resize(nElGauss);
                el_gPhi_jj.at(iElem).resize(nElGauss);

                gaussPtCart.at(iElem).resize(nElGauss);
                BMat.at(iElem).resize(nElGauss);

                // Loop through nodes to get coordinates.
                for(int iNod=0; iNod<nElNodes; iNod++){

                    dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(0);
                    dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(1);

                    dummyElNod_phi_j[iNod] = nod_phi_j.at(elemNodeConn.at(iElem).at(iNod));
                    dummyElNod_gPhi_ii[iNod] = nod_gPhi_ii.at(elemNodeConn.at(iElem).at(iNod));
                    dummyElNod_gPhi_ij[iNod] = nod_gPhi_ij.at(elemNodeConn.at(iElem).at(iNod));
                    dummyElNod_gPhi_jj[iNod] = nod_gPhi_jj.at(elemNodeConn.at(iElem).at(iNod));
                }

                elemNodCoord.at(iElem) = dummyElNodCoord;

                // Loop through integration points.
                for(int iGauss=0; iGauss<nElGauss; iGauss++){
                
                    // Cart coord of iGauss point.
                    dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
                    // Derivatives and int-pt volume
                    CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss));
                    // Int-pt phi
                    el_phi_j.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_phi_j;
                    el_gPhi_ii.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_gPhi_ii;
                    el_gPhi_ij.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_gPhi_ij;
                    el_gPhi_jj.at(iElem).at(iGauss)  = shapeFunc.at(iGauss)*dummyElNod_gPhi_jj;

                }
                gaussPtCart.at(iElem) = dummyElemGauss;
                intPtVol.at(iElem) = dummyIntVol;
            }

        } else if (Trapping=="MechTrapping" || Trapping=="MechTrappingPFF") {         // Stresses and dislocations 
            
            nod_sigma_h.resize(nNodes);
            nod_rho.resize(nNodes);

            // Loop through elements.
            for(int iElem=0; iElem<nElements; iElem++){

                accessVec(elFlux, iElem).resize(nElGauss);
                accessVec(elCon, iElem).resize(nElGauss);

                gaussPtCart.at(iElem).resize(nElGauss);
                BMat.at(iElem).resize(nElGauss);

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
                    CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss));

                }
                gaussPtCart.at(iElem) = dummyElemGauss;
                intPtVol.at(iElem) = dummyIntVol;
        
            }

        } else {
            throw std::runtime_error("Undefined trapping type < " + Trapping + " >");
        }

    } catch (const std::runtime_error& e) {

        logger.log("\nException caught in Quad4TH::InitializeElements:\n", "", false);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

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

    for(int iElem=0; iElem<nElements; iElem++){

        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            std::get<std::vector<ColVecd4>>(glIntPtCoords).at(elemIDs.at(iElem)*nElGauss+iGaus)[0] =  gaussPtCart.at(iElem).at(iGaus)[0];
            std::get<std::vector<ColVecd4>>(glIntPtCoords).at(elemIDs.at(iElem)*nElGauss+iGaus)[1] =  gaussPtCart.at(iElem).at(iGaus)[1];
            std::get<std::vector<ColVecd4>>(glIntPtCoords).at(elemIDs.at(iElem)*nElGauss+iGaus)[2] =  0;

        }
    }
    
}

// void Quad4TH::CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi){

//     // Int-pt gradients of phi.
//     vector<vector<ColVecd2>> elGradPhi(nElements);// For element nodal values of phi.
//     ColVecd4 dummyElNodPhi;        // For element nodal values of phi.
//     ColVecd4 dummyElNodLapPhi;  // For element nodal values of phi laplacian.
//     int iNode;  // counter for the number of nodes.

//     for(int iElem=0; iElem<nElements; iElem++){

//         elGradPhi.at(iElem).resize(nElGauss);

//         for(int iGaus=0; iGaus<nElGauss; iGaus++){

//             const Matd2x4& dummyBMat = BMat.at(iElem).at(iGaus); // derivative matrix for the given gauss point.
//             elGradPhi.at(iElem).at(iGaus) = dummyBMat*dummyElNodPhi;
//             dummyElNodLapPhi = dummyBMat.transpose()*dummyBMat*dummyElNodPhi;

//             iNode = 0;
//             for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

//                 std::get<std::vector<ColVecd2>>(nodGrad).at(*iNod2) += elGradPhi.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
//                 nodCount.at(*iNod2) += shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
//                 nodLapPhi[*iNod2]   += dummyElNodLapPhi[iNode]; 
//                 iNode += 1;
//             }
//         }
//     }
// }

void Quad4TH::CalcElemStiffMatx(BaseTrapping* mat, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){

    Matd2x2 DMat; 

    vector<Matd4x4> elKDMatx(nElements); 
    vector<Matd4x4> elKTMatx(nElements);
    vector<Matd4x4> elKZMatx(nElements);

    Matd4x4 BDB;

    double dummydVol;       // dummy for int-pt volume.

    try{

        if (Trapping=="GBTrapping"){                   // GB

            Matd2x2 TMat; 

            ColVecd4 dummyElNod_gPhi;  // For element nodal values of phi.
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

                    DMat = std::get<Matd2x2>(mat->CalcDMatx(gPhi, T));
                    TMat = std::get<Matd2x2>(dynamic_cast<TrapGB*>(mat)->CalcTMatx(gPhi, T));

                    const Matd2x4& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
                    const RowVecd4& dummyShFunc = shapeFunc.at(iGauss);
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

        } else if (Trapping=="2PhaseTrapping") {       // 2Phase 

            ColVecd4 dummyElNod_phi_j, dummyElNod_gPhi_ij, dummyElNod_gPhi_ii,
            dummyElNod_gPhi_jj;
            double phi_j, gPhi_ii, gPhi_ij, gPhi_jj;              

            vector<int> NodeConn;

            double zeta_j  = dynamic_cast<TrapPhase*>(mat)->get_zeta_j();
            double zeta_ij = dynamic_cast<TrapPhase*>(mat)->get_zeta_ij();
            double zeta_jj = dynamic_cast<TrapPhase*>(mat)->get_zeta_jj();
            double zeta_ii = dynamic_cast<TrapPhase*>(mat)->get_zeta_ii();

            // Loop through all elements.
            for(int iElem=0; iElem<nElements; iElem++){

                NodeConn = elemNodeConn.at(iElem);

                // MUST BE POPULATED WITH ZEROS    
                elStiffMatx.at(iElem).setZero();
                elKDMatx.at(iElem).setZero();  
                elKTMatx.at(iElem).setZero();  
                elCapMatx.at(iElem).setZero(); 

                // Loop through element nodes to get nodal values.
                for(int iNod=0; iNod<nElNodes; iNod++){
                    dummyElNod_phi_j[iNod] = nod_phi_j.at(NodeConn.at(iNod));
                    dummyElNod_gPhi_ii[iNod] = nod_gPhi_ii.at(NodeConn.at(iNod));
                    dummyElNod_gPhi_ij[iNod] = nod_gPhi_ij.at(NodeConn.at(iNod));
                    dummyElNod_gPhi_jj[iNod] = nod_gPhi_jj.at(NodeConn.at(iNod));
                }              

                // Integration over all Gauss points.
                for (int iGauss=0; iGauss<nElGauss; iGauss++){

                    phi_j   = el_phi_j.at(iElem).at(iGauss);  // values of the current int-pt
                    gPhi_ii = el_gPhi_ii.at(iElem).at(iGauss);
                    gPhi_ij = el_gPhi_ij.at(iElem).at(iGauss);
                    gPhi_jj = el_gPhi_jj.at(iElem).at(iGauss);

                    DMat = std::get<Matd2x2>(mat->CalcDMatx(phi_j, T));

                    const Matd2x4& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
                    const RowVecd4& dummyShFunc = shapeFunc.at(iGauss);
                    dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current int-pt 

                    // [B_ji]^T k_jj B_ji
                    elKDMatx.at(iElem).noalias() += dummyBMat.transpose()*DMat*dummyBMat*dummydVol;
                    // [B_ji]^T k_jj B_ji
                    elKTMatx.at(iElem).noalias() += dummyBMat.transpose()*DMat*zeta_j/(R*T)*dummyBMat*dummyElNod_phi_j*dummyShFunc*dummydVol + 
                    dummyBMat.transpose()*DMat*zeta_ii/(R*T)*dummyBMat*dummyElNod_gPhi_ii*dummyShFunc*dummydVol + 
                    dummyBMat.transpose()*DMat*zeta_ij/(R*T)*dummyBMat*dummyElNod_gPhi_ij*dummyShFunc*dummydVol +
                    dummyBMat.transpose()*DMat*zeta_jj/(R*T)*dummyBMat*dummyElNod_gPhi_jj*dummyShFunc*dummydVol; 
                    // [N_i]^T N_i
                    elCapMatx.at(iElem).noalias() += (dummyShFunc.transpose()*dummyShFunc)*dummydVol;
                }

                elStiffMatx.at(iElem) = dt*elKDMatx.at(iElem) - dt*elKTMatx.at(iElem) + elCapMatx.at(iElem);
            }
        
        } else if (Trapping=="MechTrapping") {         // Stresses and dislocations 

            MechTrap* mechTrapMat = dynamic_cast<MechTrap*>(mat);
            double s = mechTrapMat->get_s();
            double Vh = mechTrapMat->get_Vh();
            double zeta_rho = mechTrapMat->get_zeta_rho();

            ColVecd4 dummyElNod_sigma_h, dummyElNod_rho;
            double rho;

            vector<int> NodeConn;

            // Loop through all elements.
            for(int iElem=0; iElem<nElements; iElem++){

                NodeConn = elemNodeConn.at(iElem);

                // MUST BE POPULATED WITH ZEROS    
                accessVec(elStiffMatx, iElem).setZero();
                accessVec(elKDMatx, iElem).setZero();
                accessVec(elKTMatx, iElem).setZero();
                accessVec(elCapMatx, iElem).setZero();
                BDB.setZero();

                // Loop through element nodes to get nodal values.
                for(int iNod=0; iNod<nElNodes; iNod++){
                    dummyElNod_sigma_h[iNod] = nod_sigma_h.at(NodeConn.at(iNod));
                    dummyElNod_rho[iNod] = nod_rho.at(NodeConn.at(iNod));
                }       
                
                // Integration over all Gauss points.
                for (int iGauss=0; iGauss<nElGauss; iGauss++){

                    rho = accessVec(shapeFunc, iGauss)*dummyElNod_rho;
                    DMat = std::get<Matd2x2>(mechTrapMat->CalcDMatx(rho, T)); 

                    const Matd2x4& dummyBMat = accessVec(BMat, iElem, iGauss); // derivative matrix for the given gauss point.
                    const RowVecd4& dummyShFunc = accessVec(shapeFunc, iGauss);
                    dummydVol = accessVec(intPtVol, iElem, iGauss);  // Volume of the current int-pt

                    BDB = dummyBMat.transpose()*DMat*dummyBMat*dummydVol;

                    // [B_ji]^T k_jj B_ji
                    accessVec(elKDMatx, iElem).noalias() += BDB;
                    // [B_ji]^T k_jj B_ji
                    accessVec(elKTMatx, iElem).noalias() += BDB*(Vh/(R*T)*dummyElNod_sigma_h*dummyShFunc
                                                                 + zeta_rho/(R*T)*dummyElNod_rho*dummyShFunc);
                    // [N_i]^T N_i
                    elCapMatx.at(iElem).noalias() += s*(dummyShFunc.transpose()*dummyShFunc)*dummydVol;

                }

                accessVec(elStiffMatx, iElem) = dt*accessVec(elKDMatx, iElem) - dt*accessVec(elKTMatx, iElem) + accessVec(elCapMatx, iElem);
            }

        } else if (Trapping=="MechTrappingPFF") {         // PFF 

            MechTrap* mechTrapMat = dynamic_cast<MechTrap*>(mat);
            double s = mechTrapMat->get_s();
            double Vh = mechTrapMat->get_Vh();
            double zeta_rho = mechTrapMat->get_zeta_rho();
            double Zd = mechTrapMat->get_Zd();

            ColVecd4 dummyElNod_sigma_h, dummyElNod_rho;
            double intPtRho;
            double phi2;

            vector<int> NodeConn;

            // Loop through all elements.
            for(int iElem=0; iElem<nElements; iElem++){

                NodeConn = elemNodeConn.at(iElem);

                // MUST BE POPULATED WITH ZEROS    
                accessVec(elStiffMatx, iElem).setZero();
                accessVec(elKDMatx, iElem).setZero();
                accessVec(elKTMatx, iElem).setZero();
                accessVec(elKZMatx, iElem).setZero();
                accessVec(elCapMatx, iElem).setZero();
                BDB.setZero();

                // Loop through element nodes to get nodal values.
                for(int iNod=0; iNod<nElNodes; iNod++){
                    dummyElNod_sigma_h[iNod] = nod_sigma_h.at(NodeConn.at(iNod));
                    dummyElNod_rho[iNod] = nod_rho.at(NodeConn.at(iNod));
                }       
                
                // Integration over all Gauss points.
                for (int iGauss=0; iGauss<nElGauss; iGauss++){

                    const Matd2x4& dummyBMat = accessVec(BMat, iElem, iGauss); // derivative matrix for the given gauss point.
                    const RowVecd4& dummyShFunc = accessVec(shapeFunc, iGauss);
                    dummydVol = accessVec(intPtVol, iElem, iGauss);  // Volume of the current int-pt

                    intPtRho = dummyShFunc*dummyElNod_rho;
                    phi2 = accessVec(*elPhi_d_ptr, iElem, iGauss)*accessVec(*elPhi_d_ptr, iElem, iGauss);
                    DMat = std::get<Matd2x2>(mechTrapMat->CalcDMatx(intPtRho, T))*(1 - 0.99*phi2); 
                    BDB = dummyBMat.transpose()*DMat*dummyBMat*dummydVol;

                    // [B_ji]^T k_jj B_ji
                    accessVec(elKDMatx, iElem).noalias() += BDB;
                    // [B_ji]^T k_jj B_ji
                    accessVec(elKTMatx, iElem).noalias() += BDB*(Vh/(R*T)*dummyElNod_sigma_h*dummyShFunc
                                                                 + zeta_rho/(R*T)*dummyElNod_rho*dummyShFunc);

                    accessVec(elKZMatx, iElem).noalias() += Zd*phi2*(dummyShFunc.transpose()*dummyShFunc)*dummydVol;

                    // [N_i]^T N_i
                    elCapMatx.at(iElem).noalias() += s*(dummyShFunc.transpose()*dummyShFunc)*dummydVol;

                }

                accessVec(elStiffMatx, iElem) = dt*accessVec(elKDMatx, iElem) - dt*accessVec(elKTMatx, iElem)
                                              + dt*accessVec(elKZMatx, iElem) + accessVec(elCapMatx, iElem);
            }
        } else {
            throw std::runtime_error("Undefined trapping type < " + Trapping + " >");
        }

    } catch (const std::runtime_error& e) {

        logger.log("\nException caught in Quad4TH::CalcElemStiffMatx:\n", "", false);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

    }

    // // TODO: For debug!
    // for (auto& iStifMat : elKDMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
    elCapMatxVariant = &elCapMatx;
}

double Quad4TH::CalcAvCon(const double* globalBuffer){

    double IntPtCon, AvCon = 0, TotVol = 0;
    ColVecd4 dummyCon; 

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal concentration from the solution vector. 
        for(int iDof=0; iDof<nElConDofs; iDof++){
            dummyCon[iDof] = globalBuffer[elemConDof.at(iElem).at(iDof)];
        }

        // Gauss points
        for(int iGaus=0; iGaus<nElGauss; iGaus++){
            IntPtCon = shapeFunc.at(iGaus)*dummyCon;
            TotVol += intPtVol.at(iElem).at(iGaus);
            AvCon += IntPtCon*intPtVol.at(iElem).at(iGaus);
        }
    }

    return AvCon/TotVol;
}

void Quad4TH::CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){

    Matd2x2 DMat; 
    double IntPtCon;   // Integration point concnetration
    ColVecd4 dummyCon; // for element nodal concentration.
    int iNode;         // counter for the number of nodes.
    Matd2x4 DB;

    try{
        if (Trapping=="MechTrapping"){       // hydrostatic stresses an dislocations

            MechTrap* mechTrapMat = dynamic_cast<MechTrap*>(mat);
            double Vh = mechTrapMat->get_Vh();
            double zeta_rho = mechTrapMat->get_zeta_rho();
            ColVecd4 dummyElNod_sigma_h, dummyElNod_rho;
            double intPtRho;

            // Integration point values.
            for(int iElem=0; iElem<nElements; iElem++){

                // Loop through element nodes to get nodal values.
                for(int iDof=0; iDof<nElConDofs; iDof++){
                    dummyCon[iDof] = globalBuffer[accessVec(elemConDof, iElem, iDof)];
                    dummyElNod_sigma_h[iDof] = accessVec(nod_sigma_h, accessVec(elemConDof, iElem, iDof));
                    dummyElNod_rho[iDof] = accessVec(nod_rho, accessVec(elemConDof, iElem, iDof));
                }

                // Gauss points
                for(int iGauss=0; iGauss<nElGauss; iGauss++){

                    const RowVecd4& dummyShFunc = accessVec(shapeFunc, iGauss);

                    IntPtCon = dummyShFunc*dummyCon;
                    intPtRho = dummyShFunc*dummyElNod_rho;
                    DMat = std::get<Matd2x2>(mechTrapMat->CalcDMatx(intPtRho, T));
                    DB = DMat*accessVec(BMat, iElem, iGauss);

                    // Int pt flux
                    accessVec(elFlux, iElem, iGauss) = DB*(- dummyCon + Vh/(R*T)*IntPtCon*dummyElNod_sigma_h
                                                                      + zeta_rho/(R*T)*dummyElNod_rho);

                    std::get<std::vector<ColVecd3>>(intPtFlux).at(accessVec(elemIDs, iElem)*nElGauss+iGauss)[0] =  accessVec(elFlux, iElem, iGauss)[0];
                    std::get<std::vector<ColVecd3>>(intPtFlux).at(accessVec(elemIDs, iElem)*nElGauss+iGauss)[1] =  accessVec(elFlux, iElem, iGauss)[1];
                    std::get<std::vector<ColVecd3>>(intPtFlux).at(accessVec(elemIDs, iElem)*nElGauss+iGauss)[2] =  0;

                    // Nodal values
                    iNode = 0;
                    for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                        std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += accessVec(elFlux, iElem, iGauss)*dummyShFunc[iNode]*accessVec(wts, iGauss);
                        accessVec(nodCount, *iNod2) += dummyShFunc[iNode]*accessVec(wts, iGauss);
                        iNode += 1;
                    }
                }
            }

        } else if (Trapping == "MechTrappingPFF"){       // PFF

            MechTrap* mechTrapMat = dynamic_cast<MechTrap*>(mat);

            double Vh = mechTrapMat->get_Vh();
            double zeta_rho = mechTrapMat->get_zeta_rho();

            ColVecd4 dummyElNod_sigma_h, dummyElNod_rho;
            double intPtRho;

            // Integration point values.
            for(int iElem=0; iElem<nElements; iElem++){

                // Loop through element nodes to get nodal values.
                for(int iDof=0; iDof<nElConDofs; iDof++){
                    dummyCon[iDof] = globalBuffer[accessVec(elemConDof, iElem, iDof)];
                    dummyElNod_sigma_h[iDof] = accessVec( nod_sigma_h, accessVec(elemConDof, iElem, iDof));
                    dummyElNod_rho[iDof] = accessVec( nod_rho, accessVec(elemConDof, iElem, iDof));
                }

                // Gauss points
                for(int iGauss=0; iGauss<nElGauss; iGauss++){

                    const RowVecd4& N_i = accessVec(shapeFunc, iGauss);

                    IntPtCon = N_i*dummyCon;
                    double phi2 = accessVec(*elPhi_d_ptr, iElem, iGauss)*accessVec(*elPhi_d_ptr, iElem, iGauss);
                    intPtRho = N_i*dummyElNod_rho;
                    DMat = std::get<Matd2x2>(mechTrapMat->CalcDMatx(intPtRho, T))*(1 - 0.99*phi2);
                    DB = DMat*accessVec(BMat, iElem, iGauss);

                    double dummydVol = accessVec(intPtVol, iElem, iGauss);

                    // Int pt flux
                    accessVec(elFlux, iElem, iGauss) = DB*( -dummyCon + Vh/(R*T)*IntPtCon*dummyElNod_sigma_h 
                                                                      + zeta_rho/(R*T)*dummyElNod_rho);

                    std::get<std::vector<ColVecd3>>(intPtFlux).at(accessVec(elemIDs, iElem)*nElGauss+iGauss)[0] =  accessVec(elFlux, iElem, iGauss)[0];
                    std::get<std::vector<ColVecd3>>(intPtFlux).at(accessVec(elemIDs, iElem)*nElGauss+iGauss)[1] =  accessVec(elFlux, iElem, iGauss)[1];
                    std::get<std::vector<ColVecd3>>(intPtFlux).at(accessVec(elemIDs, iElem)*nElGauss+iGauss)[2] =  0;

                    // Nodal values
                    iNode = 0;
                    for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                        std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += accessVec(elFlux, iElem, iGauss)*N_i[iNode]*dummydVol;
                        accessVec(nodCount, *iNod2) += N_i[iNode]*dummydVol;
                        iNode += 1;
                    }
                }
            }

        } else {
            throw std::runtime_error("Undefined trapping type < " + Trapping + " >");
        }

    } catch (const std::runtime_error& e) {

        logger.log("\nException caught in Quad4TH::CalcFlux:\n", "ERROR", false);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

    }

    // // // TODO: For debug!
    // // cout << elFlux.at(0).at(0) << "\n\n";
}

void Quad4TH::CalcFsrc(const double conB, BaseTrapping* mat, double* FsrcBuffer, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){

    MechTrap* mechTrapMat = dynamic_cast<MechTrap*>(mat);

    double Vh = mechTrapMat->get_Vh();
    double Zd = mechTrapMat->get_Zd();
    double zeta_rho = mechTrapMat->get_zeta_rho();
    ColVecd4 dummyElNod_sigma_h, dummyElNod_rho, Ceq, dummyFsrc;

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){

        // Loop through element nodes to get nodal values.
        for(int iDof=0; iDof<nElConDofs; iDof++){
            dummyElNod_sigma_h[iDof] = accessVec( nod_sigma_h, accessVec(elemConDof, iElem, iDof));
            dummyElNod_rho[iDof] = accessVec( nod_rho, accessVec(elemConDof, iElem, iDof));
        }

        // Gauss points
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            // Equilibrium conentration
            Ceq = conB * ((dummyElNod_sigma_h * Vh + dummyElNod_rho * zeta_rho ) / (R * T)).array().exp().matrix();

            const RowVecd4& N_i = accessVec(shapeFunc, iGaus);
            const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  

            double phi2 = accessVec(*elPhi_d_ptr, iElem, iGaus)*accessVec(*elPhi_d_ptr, iElem, iGaus);
            Matd4x4 NTN = N_i.transpose() * N_i;
            dummyFsrc = dt*Zd*phi2*NTN*Ceq*dummydVol;

            for(int iNod2 = 0; iNod2 < nElConDofs; ++iNod2) {
                int gNod = elemNodeConn[iElem][iNod2];
                FsrcBuffer[gNod] += dummyFsrc[iNod2];
            }
        }
    }
}

void Quad4TH::CalcElCon(const double* globalBuffer){

    ColVecd4 nodCon;

    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElConDofs; iDof++){
            nodCon(iDof) = globalBuffer[accessVec(elemConDof, iElem, iDof)];
        }

        for(int iGauss=0; iGauss<nElGauss; iGauss++){

            const RowVecd4& N_i = accessVec(shapeFunc, iGauss);
            accessVec(elCon, iElem, iGauss) = N_i*nodCon;

        }
    }
}