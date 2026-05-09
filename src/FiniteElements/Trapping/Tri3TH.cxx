#include<iostream>
#include<algorithm>

#include "FiniteElements/Trapping/Tri3TH.h"
#include "Materials/Trapping/TrapGB.h"
#include "Materials/Trapping/TrapHLGB.h"
#include "Materials/Trapping/TrapPhase.h"

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

Tri3TH::Tri3TH(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, Logger& logger,  H5IO* H5File_rve)
    : BaseElemTrap(2, 3, 3, 3, logger){ // nElDim, nElNodes, nElGauss, nElConDofs

    InitShapeFunc();
    ReadElementsData(H5File_in, H5File_mesh, iSet);
    InitializeElements(Nodes, H5File_rve);
}   

Tri3TH::~Tri3TH(){}

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

void Tri3TH::InitializeElements(Nodes &Nodes, H5IO* H5File_rve){

    // Initialize the storages
    elFlux.resize(nElements);
    elStiffMatx.resize(nElements);
    elCapMatx.resize(nElements);

    elCon.resize(nElements);
    elConTrap.resize(nElements); elDeltaTrap.resize(nElements);

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd3x2 dummyElNodCoord;   // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    if (Trapping != "MechTrapping" && H5File_rve == nullptr) {
    logger.log("Error in Tri3TH::InitializeElements\n", "ERROR", true);
    logger.log("    Physics Category '" + Trapping + "' requires external RVE data,", "", false);
    logger.log("    but the RVE HDF5 file pointer is null.", "", false);
    logger.log("    Check if the RVE file exists and is correctly passed.", "", false);
    logger.log("    Critical error encountered. Terminating!\n", "", false);
    exit(EXIT_FAILURE);
    }

    try{

        if (Trapping=="GBTrapping"){    // GB

            ColVecd3 dummyElNod_gPhi;   // For element nodal values of phi.

            // Initialize the storages for int-pt gPhi
            el_gPhi.resize(nElements);
            // Read nodal values of gPhi
            nod_gPhi.resize(nNodes);
            string dsetName = "gPhi";
            H5File_rve->ReadField1D(dsetName, nod_gPhi);

            // Loop through elements.
            for(int iElem=0; iElem<nElements; iElem++){

                elFlux.at(iElem).resize(nElGauss);
                accessVec(elCon, iElem).resize(nElGauss);
                accessVec(elConTrap, iElem).resize(nElGauss);
                accessVec(elDeltaTrap, iElem).resize(nElGauss);

                el_gPhi.at(iElem).resize(nElGauss);

                gaussPtCart.at(iElem).resize(nElGauss);
                BMat.at(iElem).resize(nElGauss);

                // Loop through nodes to get coordinates.
                for(int iNod=0; iNod<nElNodes; iNod++){

                    dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(0);
                    dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(1);

                    dummyElNod_gPhi[iNod] = accessVec(nod_gPhi, accessVec(elemNodeConn, iElem, iNod));
                }

                elemNodCoord.at(iElem) = dummyElNodCoord;

                // Loop through integration points.
                for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                    // Cart coord of iGaus point.
                    dummyElemGauss.at(iGaus) = getGaussCart(accessVec(shapeFunc, iGaus), dummyElNodCoord);
                    // Derivatives and int-pt volume
                    CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGaus), wts.at(iGaus), dummyIntVol.at(iGaus), accessVec(BMat, iElem, iGaus));
                    // Int-pt phi
                    accessVec(el_gPhi, iElem, iGaus)  = accessVec(shapeFunc, iGaus)*dummyElNod_gPhi;
                }
                gaussPtCart.at(iElem) = dummyElemGauss;
                intPtVol.at(iElem) = dummyIntVol;
            }  

        } else if (Trapping=="HLGBTrapping"){           //HLGB

            ColVecd3 dummyElNod_gPhi_HAGB, dummyElNod_gPhi_LAGB;

            // Initialize the storages for int-pt gPhi
            el_gPhi_HAGB.resize(nElements); el_gPhi_LAGB.resize(nElements);
            // Initialize the storages for nodal values of gPhi
            nod_gPhi_HAGB.resize(nNodes); nod_gPhi_LAGB.resize(nNodes);

            // Read gPhis
            string dsetName = "gPhi_HAGB";
            H5File_rve->ReadField1D(dsetName, nod_gPhi_HAGB);
            dsetName = "gPhi_LAGB";
            H5File_rve->ReadField1D(dsetName, nod_gPhi_LAGB);

            // Loop through elements.
            for(int iElem=0; iElem<nElements; iElem++){

                elFlux.at(iElem).resize(nElGauss);
                el_gPhi_HAGB.at(iElem).resize(nElGauss);
                el_gPhi_LAGB.at(iElem).resize(nElGauss);

                gaussPtCart.at(iElem).resize(nElGauss);
                BMat.at(iElem).resize(nElGauss);

                // Loop through nodes to get coordinates.
                for(int iNod=0; iNod<nElNodes; iNod++){

                    dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(0);
                    dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(1);

                    dummyElNod_gPhi_HAGB[iNod] = accessVec(nod_gPhi_HAGB, accessVec(elemNodeConn, iElem, iNod));
                    dummyElNod_gPhi_LAGB[iNod] = accessVec(nod_gPhi_LAGB, accessVec(elemNodeConn, iElem, iNod));
                }

                elemNodCoord.at(iElem) = dummyElNodCoord;

                // Loop through integration points.
                for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                    // Cart coord of iGaus point.
                    dummyElemGauss.at(iGaus) = getGaussCart(accessVec(shapeFunc, iGaus), dummyElNodCoord);
                    // Derivatives and int-pt volume
                    CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGaus), wts.at(iGaus), dummyIntVol.at(iGaus), accessVec(BMat, iElem, iGaus));
                    // Int-pt phi
                    accessVec(el_gPhi_HAGB, iElem, iGaus)  = accessVec(shapeFunc, iGaus)*dummyElNod_gPhi_HAGB;
                    accessVec(el_gPhi_HAGB, iElem, iGaus)  = accessVec(shapeFunc, iGaus)*dummyElNod_gPhi_LAGB;
                }
                gaussPtCart.at(iElem) = dummyElemGauss;
                intPtVol.at(iElem) = dummyIntVol;
            }   
                   
        } else if (Trapping=="2PhaseTrapping") {       // 2Phase 

            ColVecd3 dummyElNod_phi_j; ColVecd3 dummyElNod_gPhi_jj;
            ColVecd3 dummyElNod_gPhi_ii; ColVecd3 dummyElNod_gPhi_ij;

            // Initialize the storages for int-pt traps
            el_phi_j.resize(nElements); el_gPhi_ii.resize(nElements);
            el_gPhi_jj.resize(nElements); el_gPhi_ij.resize(nElements);

            // Read nodal values of traps
            nod_phi_j.resize(nNodes);
            string dsetName = "phi_j";
            H5File_rve->ReadField1D(dsetName, nod_phi_j);

            nod_gPhi_ii.resize(nNodes);
            dsetName = "gPhi_ii";
            H5File_rve->ReadField1D(dsetName, nod_gPhi_ii);

            nod_gPhi_ij.resize(nNodes);
            dsetName = "gPhi_ij";
            H5File_rve->ReadField1D(dsetName, nod_gPhi_ij);

            nod_gPhi_jj.resize(nNodes);
            dsetName = "gPhi_jj";
            H5File_rve->ReadField1D(dsetName, nod_gPhi_jj);

            // Loop through elements.
            for(int iElem=0; iElem<nElements; iElem++){

                elFlux.at(iElem).resize(nElGauss);

                accessVec(elCon, iElem).resize(nElGauss);
                accessVec(elConTrap, iElem).resize(nElGauss);
                accessVec(elDeltaTrap, iElem).resize(nElGauss);

                el_phi_j.at(iElem).resize(nElGauss);
                el_gPhi_ii.at(iElem).resize(nElGauss);
                el_gPhi_ij.at(iElem).resize(nElGauss);
                el_gPhi_jj.at(iElem).resize(nElGauss);

                gaussPtCart.at(iElem).resize(nElGauss);
                BMat.at(iElem).resize(nElGauss);

                // Loop through nodes to get coordinates.
                for(int iNod=0; iNod<nElNodes; iNod++){

                    dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(0);
                    dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(1);

                    dummyElNod_phi_j[iNod] = nod_phi_j.at(accessVec(elemNodeConn, iElem, iNod));
                    dummyElNod_gPhi_ii[iNod] = nod_gPhi_ii.at(accessVec(elemNodeConn, iElem, iNod));
                    dummyElNod_gPhi_ij[iNod] = nod_gPhi_ij.at(accessVec(elemNodeConn, iElem, iNod));
                    dummyElNod_gPhi_jj[iNod] = nod_gPhi_jj.at(accessVec(elemNodeConn, iElem, iNod));
                }

                elemNodCoord.at(iElem) = dummyElNodCoord;

                // Loop through integration points.
                for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                    // Cart coord of iGaus point.
                    dummyElemGauss.at(iGaus) = getGaussCart(accessVec(shapeFunc, iGaus), dummyElNodCoord);
                    // Derivatives and int-pt volume
                    CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGaus), wts.at(iGaus), dummyIntVol.at(iGaus), accessVec(BMat, iElem, iGaus));
                    // Int-pt phi
                    accessVec(el_phi_j, iElem, iGaus)  = accessVec(shapeFunc, iGaus)*dummyElNod_phi_j;
                    accessVec(el_gPhi_ii, iElem, iGaus)  = accessVec(shapeFunc, iGaus)*dummyElNod_gPhi_ii;
                    accessVec(el_gPhi_ij, iElem, iGaus)  = accessVec(shapeFunc, iGaus)*dummyElNod_gPhi_ij;
                    accessVec(el_gPhi_jj, iElem, iGaus)  = accessVec(shapeFunc, iGaus)*dummyElNod_gPhi_jj;

                }
                gaussPtCart.at(iElem) = dummyElemGauss;
                intPtVol.at(iElem) = dummyIntVol;
            }
        } else if (Trapping=="MechTrapping") {         // Stresses and dislocations 

            nod_sigma_h.resize(nNodes);

            // Loop through elements.
            for(int iElem=0; iElem<nElements; iElem++){

                elFlux.at(iElem).resize(nElGauss);

                gaussPtCart.at(iElem).resize(nElGauss);
                BMat.at(iElem).resize(nElGauss);

                // Loop through nodes to get coordinates.
                for(int iNod=0; iNod<nElNodes; iNod++){

                    dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(0);
                    dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(accessVec(elemNodeConn, iElem, iNod)).at(1);

                }

                elemNodCoord.at(iElem) = dummyElNodCoord;

                // Loop through integration points.
                for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                    // Cart coord of iGaus point.
                    dummyElemGauss.at(iGaus) = getGaussCart(accessVec(shapeFunc, iGaus), dummyElNodCoord);
                    // Derivatives and int-pt volume
                    CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGaus), wts.at(iGaus), dummyIntVol.at(iGaus), accessVec(BMat, iElem, iGaus));

                }
                gaussPtCart.at(iElem) = dummyElemGauss;
                intPtVol.at(iElem) = dummyIntVol;
        
            }
        }
    } catch (const std::runtime_error& e) {
        logger.log("\nException caught in Tri3TH::InitializeElements:\n", "", false);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("\nCritical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
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

void Tri3TH::getInPtCoords(T_nodStres& glIntPtCoords){

    for(int iElem=0; iElem<nElements; iElem++){

        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            std::get<std::vector<ColVecd3>>(glIntPtCoords).at(accessVec(elemIDs,iElem)*nElGauss+iGaus)[0] =  gaussPtCart.at(iElem).at(iGaus)[0];
            std::get<std::vector<ColVecd3>>(glIntPtCoords).at(accessVec(elemIDs,iElem)*nElGauss+iGaus)[1] =  gaussPtCart.at(iElem).at(iGaus)[1];
            std::get<std::vector<ColVecd3>>(glIntPtCoords).at(accessVec(elemIDs,iElem)*nElGauss+iGaus)[2] =  0;

        }
    }
}

// void Tri3TH::CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi){

//     // Int-pt gradients of phi.
//     vector<vector<ColVecd2>> elGrad_gPhi(nElements);    // For element nodal values of gPhi.
//     ColVecd3 dummyElNodPhi;        // For element nodal values of phi.
//     ColVecd3 dummyElNodLapPhi;     // For element nodal values of phi laplacian.
//     int iNode;                     // counter for the number of nodes.

//     for(int iElem=0; iElem<nElements; iElem++){

//         elGrad_gPhi.at(iElem).resize(nElGauss);

//         // Loop through element nodes to get nodal values.
//         for(int iNod=0; iNod<nElNodes; iNod++){
//             dummyElNodPhi[iNod] = accessVec(nod_gPhi, accessVec(elemNodeConn, iElem, iNod));
//         }

//         for(int iGaus=0; iGaus<nElGauss; iGaus++){

//             const Matd2x3& dummyBMat = accessVec(BMat, iElem, iGaus); // derivative matrix for the given gauss point.
//             elGrad_gPhi.at(iElem).at(iGaus) = dummyBMat*dummyElNodPhi;
//             dummyElNodLapPhi = dummyBMat.transpose()*dummyBMat*dummyElNodPhi;

//             iNode = 0;
//             for(auto iNod2=accessVec(elemNodeConn,iElem).begin(); iNod2!=accessVec(elemNodeConn,iElem).end(); iNod2++){

//                 std::get<std::vector<ColVecd2>>(nodGrad).at(*iNod2) += elGrad_gPhi.at(iElem).at(iGaus);
//                 nodCount.at(*iNod2) += 1;
//                 nodLapPhi[*iNod2]   += dummyElNodLapPhi[iNode]; 
//                 iNode += 1;
//             }
//         }
//     }
// }

void Tri3TH::CalcElemStiffMatx(BaseTrapping* mat, const double T, const std::vector<std::vector<double>>* elPhi_d_pt){

    Matd2x2 DMat; 
    Matd3x3 BDB, BTB;
    vector<Matd3x3> elKDMatx(nElements); 
    vector<Matd3x3> elKTMatx(nElements); 

    // Precompute thermal constants to avoid division inside any loop
    const double invRT = 1.0 / (R * T); 

    if (Trapping=="GBTrapping"){       // GB

        TrapGB* gbTrap = dynamic_cast<TrapGB*>(mat);
        Matd2x2 TMat; 

        ColVecd3 dummyElNod_gPhi;  // For element nodal values of phi.

        // Loop through all elements.
        for(int iElem=0; iElem<nElements; iElem++){

            // MUST BE POPULATED WITH ZEROS    
            accessVec(elStiffMatx, iElem).setZero();
            accessVec(elKDMatx, iElem).setZero();  
            accessVec(elKTMatx, iElem).setZero();  
            accessVec(elCapMatx, iElem).setZero(); 

            // Loop through element nodes to get nodal values.
            for(int iNod=0; iNod<nElNodes; iNod++){
                dummyElNod_gPhi[iNod] = accessVec(nod_gPhi, accessVec(elemNodeConn, iElem, iNod));
            }              

            // Integration over all Gauss points.
            for (int iGaus=0; iGaus<nElGauss; iGaus++){

                const double& gPhi = accessVec(el_gPhi, iElem, iGaus);  // gPhi of the current int-pt

                DMat = std::get<Matd2x2>(gbTrap->CalcDMatx(gPhi, T));
                TMat = std::get<Matd2x2>(gbTrap->CalcTMatx(gPhi, T));

                const Matd2x3& dummyBMat = accessVec(BMat, iElem, iGaus); // derivative matrix for the given gauss point.
                const RowVecd3& dummyShFunc = accessVec(shapeFunc, iGaus);
                const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  // Volume of the current int-pt 
                BDB = dummyBMat.transpose()*DMat*dummyBMat*dummydVol;
                BTB = dummyBMat.transpose()*TMat*dummyBMat*dummydVol;

                // [B_ji]^T k_jj B_ji
                accessVec(elKDMatx, iElem).noalias() += BDB;
                // [B_ji]^T k_jj B_ji
                accessVec(elKTMatx, iElem).noalias() += (BTB*dummyElNod_gPhi)*dummyShFunc; 
                // [N_i]^T N_i
                accessVec(elCapMatx, iElem).noalias() += (dummyShFunc.transpose()*dummyShFunc)*dummydVol;
            }

            accessVec(elStiffMatx, iElem) = dt*accessVec(elKDMatx, iElem) - dt*accessVec(elKTMatx, iElem) + accessVec(elCapMatx, iElem);
        }

    } else if (Trapping=="HLGBTrapping"){

        TrapHLGB* hlgbTrap = dynamic_cast<TrapHLGB*>(mat);

        ColVecd3 dummyElNod_gPhi_HAGB, dummyElNod_gPhi_LAGB;  // For element nodal values of gPhi.

        // Pre-scale material properties
        const double c_zeta_HAGB = hlgbTrap->get_zeta_HAGB() * invRT;
        const double c_zeta_LAGB = hlgbTrap->get_zeta_LAGB() * invRT;

        // Loop through all elements.
        for(int iElem=0; iElem<nElements; iElem++){

            // MUST BE POPULATED WITH ZEROS    
            accessVec(elStiffMatx, iElem).setZero();
            accessVec(elKDMatx, iElem).setZero();  
            accessVec(elKTMatx, iElem).setZero();  
            accessVec(elCapMatx, iElem).setZero(); 
            BDB.setZero();

            // Loop through element nodes to get nodal values.
            for(int iNod=0; iNod<nElNodes; iNod++){
                dummyElNod_gPhi_HAGB[iNod] = accessVec(nod_gPhi_HAGB, accessVec(elemNodeConn, iElem, iNod));
                dummyElNod_gPhi_LAGB[iNod] = accessVec(nod_gPhi_LAGB, accessVec(elemNodeConn, iElem, iNod));
            }              

            // Integration over all Gauss points.
            for (int iGaus=0; iGaus<nElGauss; iGaus++){

                DMat = std::get<Matd2x2>(hlgbTrap->CalcDMatx(1.0, T));
                
                // Useless for now, but can be used for varying diffusivity
                const double& gPhi_HAGB = accessVec(el_gPhi_HAGB, iElem, iGaus);
                const double& gPhi_LAGB = accessVec(el_gPhi_LAGB, iElem, iGaus);

                const Matd2x3& dummyBMat = accessVec(BMat, iElem, iGaus); // derivative matrix for the given gauss point.
                const RowVecd3& dummyShFunc = accessVec(shapeFunc, iGaus);
                const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  // Volume of the current int-pt 
                BDB = dummyBMat.transpose()*DMat*dummyBMat*dummydVol;

                // [B_ji]^T k_jj B_ji
                accessVec(elKDMatx, iElem).noalias() += BDB;
                // [B_ji]^T k_jj B_ji
                accessVec(elKTMatx, iElem).noalias() += (BDB*dummyElNod_gPhi_HAGB)*(c_zeta_HAGB*dummyShFunc) + 
                    (BDB*dummyElNod_gPhi_LAGB)*(c_zeta_LAGB*dummyShFunc);

                // [N_i]^T N_i
                accessVec(elCapMatx, iElem).noalias() += (dummyShFunc.transpose()*dummyShFunc)*dummydVol;
            }

            accessVec(elStiffMatx, iElem) = dt*accessVec(elKDMatx, iElem) - dt*accessVec(elKTMatx, iElem) + accessVec(elCapMatx, iElem);
        }

    } else if (Trapping=="2PhaseTrapping"){        // 2Phase

        TrapPhase* phaseTrap = dynamic_cast<TrapPhase*>(mat);

        ColVecd3 dummyElNod_phi_j, dummyElNod_gPhi_ij, dummyElNod_gPhi_ii, dummyElNod_gPhi_jj;

        const double c_zeta_j  = phaseTrap->get_zeta_j() * invRT;
        const double c_zeta_ij = phaseTrap->get_zeta_ij() * invRT;
        const double c_zeta_jj = phaseTrap->get_zeta_jj() * invRT;
        const double c_zeta_ii = phaseTrap->get_zeta_ii() * invRT;             

        vector<int> NodeConn;

        // Loop through all elements.
        for(int iElem=0; iElem<nElements; iElem++){

            NodeConn = accessVec(elemNodeConn,iElem);

            // MUST BE POPULATED WITH ZEROS    
            accessVec(elStiffMatx, iElem).setZero();
            accessVec(elKDMatx, iElem).setZero();  
            accessVec(elKTMatx, iElem).setZero();  
            accessVec(elCapMatx, iElem).setZero(); 
            BDB.setZero();

            // Loop through element nodes to get nodal values.
            for(int iNod=0; iNod<nElNodes; iNod++){
                dummyElNod_phi_j[iNod] = nod_phi_j.at(NodeConn.at(iNod));
                dummyElNod_gPhi_ii[iNod] = nod_gPhi_ii.at(NodeConn.at(iNod));
                dummyElNod_gPhi_ij[iNod] = nod_gPhi_ij.at(NodeConn.at(iNod));
                dummyElNod_gPhi_jj[iNod] = nod_gPhi_jj.at(NodeConn.at(iNod));
            }              

            // Integration over all Gauss points.
            for (int iGaus=0; iGaus<nElGauss; iGaus++){

                const double& phi_j = accessVec(el_phi_j, iElem, iGaus);  // values of the current int-pt
                const double& gPhi_ii = accessVec(el_gPhi_ii, iElem, iGaus);
                const double& gPhi_ij = accessVec(el_gPhi_ij, iElem, iGaus);
                const double& gPhi_jj = accessVec(el_gPhi_jj, iElem, iGaus);

                DMat = std::get<Matd2x2>(phaseTrap->CalcDMatx(phi_j, T));

                const Matd2x3& dummyBMat = accessVec(BMat, iElem, iGaus); // derivative matrix for the given gauss point.
                const RowVecd3& dummyShFunc = accessVec(shapeFunc, iGaus);
                const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  // Volume of the current int-pt 

                BDB = dummyBMat.transpose()*DMat*dummyBMat*dummydVol;

                // [B_ji]^T k_jj B_ji
                accessVec(elKDMatx, iElem).noalias() += BDB;
                // [B_ji]^T k_jj B_ji
                accessVec(elKTMatx, iElem).noalias() += (BDB*dummyElNod_phi_j)*(c_zeta_j*dummyShFunc) + 
                    (BDB*dummyElNod_gPhi_ii)*(c_zeta_ii*dummyShFunc) + 
                    (BDB*dummyElNod_gPhi_ij)*(c_zeta_ij*dummyShFunc) +
                    (BDB*dummyElNod_gPhi_jj)*(c_zeta_jj*dummyShFunc); 
                // [N_i]^T N_i
                accessVec(elCapMatx, iElem).noalias() += (dummyShFunc.transpose()*dummyShFunc)*dummydVol;
            }

            accessVec(elStiffMatx, iElem) = dt*accessVec(elKDMatx, iElem) - dt*accessVec(elKTMatx, iElem) + accessVec(elCapMatx, iElem);
        }
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elKDMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
    elCapMatxVariant = &elCapMatx;
}

double Tri3TH::CalcAvCon(const double* globalBuffer){

    double IntPtCon, AvCon = 0, TotVol = 0;
    ColVecd3 dummyCon; 

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal concentration from the solution vector. 
        for(int iDof=0; iDof<nElConDofs; iDof++){
            dummyCon[iDof] = globalBuffer[accessVec(elemConDof,iElem,iDof)];
        }

        // Gauss points
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            const RowVecd3& N_i = accessVec(shapeFunc, iGaus);
            const double& dummydVol = accessVec(intPtVol, iElem, iGaus);

            IntPtCon = N_i*dummyCon;
            TotVol += dummydVol;
            AvCon += IntPtCon*dummydVol;
        }
    }

    return AvCon/TotVol;
}

double Tri3TH::CalcAvTrapCon(){

    double conTrap, AvCon = 0, TotVol = 0;

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            const double& dummydVol = accessVec(intPtVol, iElem, iGaus);

            conTrap = accessVec(elConTrap, iElem, iGaus);
            TotVol += dummydVol;
            AvCon += conTrap*dummydVol;
        }
    }

    return AvCon/TotVol;
}

void Tri3TH::CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){

    Matd2x2 DMat; 
    Matd2x3 DB;
    ColVecd3 dummyCon; // For element nodal concentration.
    int iNode=0;

    // Precompute thermal constants to avoid division inside any loop
    const double invRT = 1.0 / (R * T); 

    // --- OPTIMIZATION: Extract variant references ONCE outside the loops ---
    auto& ref_intPtFlux = std::get<std::vector<ColVecd3>>(intPtFlux);
    auto& ref_nodFlux   = std::get<std::vector<ColVecd2>>(nodFlux);

    if (Trapping=="GBTrapping"){       // GB
        
        TrapGB* gbTrap = dynamic_cast<TrapGB*>(mat);
        Matd2x2 TMat;
        Matd2x3 TB; 

        ColVecd3 dummyElNod_gPhi; // vector of element nodal gPhi.

        // Loop through all elements
        for(int iElem=0; iElem<nElements; iElem++){

            // Get element nodal concentration and gPhi from the solution vector. 
            for(int iDof=0; iDof<nElConDofs; iDof++){
                int gNod = accessVec(elemConDof, iElem, iDof);
                dummyCon[iDof]        = globalBuffer[gNod];
                dummyElNod_gPhi[iDof] = nod_gPhi.at(gNod);
            }

            // Gauss points
            for(int iGaus=0; iGaus<nElGauss; iGaus++){

                const double& gPhi = accessVec(el_gPhi, iElem, iGaus);
                const double& IntPtCon = accessVec(shapeFunc, iGaus)*dummyCon;
                const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  
                const Matd2x3& dummyBMat = accessVec(BMat, iElem, iGaus);
                const RowVecd3& N_i = accessVec(shapeFunc, iGaus);

                DMat = std::get<Matd2x2>(gbTrap->CalcDMatx(gPhi, T));
                TMat = std::get<Matd2x2>(gbTrap->CalcTMatx(gPhi, T));
                DB = DMat*dummyBMat;
                TB = TMat*dummyBMat;

                // Int pt flux
                auto& flux = accessVec(elFlux, iElem, iGaus);
                flux = -DB * dummyCon + TB * (IntPtCon * dummyElNod_gPhi);

                int outIdx = accessVec(elemIDs, iElem) * nElGauss + iGaus;
                ref_intPtFlux[outIdx][0] = flux[0];
                ref_intPtFlux[outIdx][1] = flux[1];
                ref_intPtFlux[outIdx][2] = 0.0;

                const auto& NodeConn = elemNodeConn[iElem];

                // Nodal values
                for(int iNode=0; iNode<nElNodes; iNode++){
                    double weight = N_i[iNode] * dummydVol;
                    int gNod = NodeConn[iNode];
                    
                    ref_nodFlux[gNod] += flux * weight;
                    nodCount.at(gNod) += weight;
                }
            }
        }

    } else if (Trapping=="HLGBTrapping"){ 

        TrapHLGB* hlgbTrap = dynamic_cast<TrapHLGB*>(mat);
        ColVecd3 dummyElNod_gPhi_HAGB, dummyElNod_gPhi_LAGB;

        // Pre-scale material properties
        const double c_zeta_HAGB = hlgbTrap->get_zeta_HAGB() * invRT;
        const double c_zeta_LAGB = hlgbTrap->get_zeta_LAGB() * invRT;

        // Loop through all elements
        for(int iElem=0; iElem<nElements; iElem++){

            // Get element nodal concentration and gPhi from the solution vector. 
            for(int iDof=0; iDof<nElConDofs; iDof++){
                int gNod = accessVec(elemConDof, iElem, iDof);
                dummyCon[iDof]             = globalBuffer[gNod];
                dummyElNod_gPhi_HAGB[iDof] = nod_gPhi_HAGB[gNod];
                dummyElNod_gPhi_LAGB[iDof] = nod_gPhi_LAGB[gNod];
            }

            // Gauss points
            for(int iGaus=0; iGaus<nElGauss; iGaus++){

                const Matd2x3& dummyBMat = accessVec(BMat, iElem, iGaus);
                const RowVecd3& N_i = accessVec(shapeFunc, iGaus);
                const double& IntPtCon = accessVec(shapeFunc, iGaus)*dummyCon;
                const double& dummydVol = accessVec(intPtVol, iElem, iGaus); 

                DMat = std::get<Matd2x2>(hlgbTrap->CalcDMatx(1.0, T));
                DB = DMat*dummyBMat;

                // Int pt flux
                auto& flux = accessVec(elFlux, iElem, iGaus);
                
                flux = -DB * (
                    -(dummyCon) 
                    + (c_zeta_HAGB * IntPtCon * dummyElNod_gPhi_HAGB)
                    + (c_zeta_LAGB * IntPtCon * dummyElNod_gPhi_LAGB)
                );

                const auto& NodeConn = elemNodeConn[iElem];

                int outIdx = accessVec(elemIDs, iElem) * nElGauss + iGaus;
                ref_intPtFlux[outIdx][0] = flux[0];
                ref_intPtFlux[outIdx][1] = flux[1];
                ref_intPtFlux[outIdx][2] = 0.0;

                for(int iNode=0; iNode<nElNodes; iNode++){
                    double weight = N_i[iNode] * dummydVol;
                    int gNod = NodeConn[iNode];
                    
                    ref_nodFlux[gNod] += flux * weight;
                    nodCount.at(gNod) += weight;
                }
            }
        }

    } else if (Trapping=="2PhaseTrapping"){        // 2Phase

        TrapPhase* phaseTrap = dynamic_cast<TrapPhase*>(mat);
        ColVecd3 dummyElNod_phi_j, dummyElNod_gPhi_ij, dummyElNod_gPhi_ii,
        dummyElNod_gPhi_jj;

        const double c_zeta_j  = phaseTrap->get_zeta_j() * invRT;
        const double c_zeta_ij = phaseTrap->get_zeta_ij() * invRT;
        const double c_zeta_jj = phaseTrap->get_zeta_jj() * invRT;
        const double c_zeta_ii = phaseTrap->get_zeta_ii() * invRT;

        // Loop through all elements
        for(int iElem=0; iElem<nElements; iElem++){

            // Loop through element nodes to get nodal values.
            for(int iDof=0; iDof<nElConDofs; iDof++){
                int gNod = accessVec(elemConDof, iElem, iDof);
                dummyCon[iDof]           = globalBuffer[gNod];
                dummyElNod_phi_j[iDof]   = nod_phi_j.at(gNod);
                dummyElNod_gPhi_ii[iDof] = nod_gPhi_ii.at(gNod);
                dummyElNod_gPhi_ij[iDof] = nod_gPhi_ij.at(gNod);
                dummyElNod_gPhi_jj[iDof] = nod_gPhi_jj.at(gNod);
            }

            // Gauss points
            for(int iGaus=0; iGaus<nElGauss; iGaus++){

                const double& IntPtCon = accessVec(shapeFunc, iGaus)*dummyCon;

                const double& phi_j = accessVec(el_phi_j, iElem, iGaus);  // values of the current int-pt
                const double& gPhi_ii = accessVec(el_gPhi_ii, iElem, iGaus);
                const double& gPhi_ij = accessVec(el_gPhi_ij, iElem, iGaus);
                const double& gPhi_jj = accessVec(el_gPhi_jj, iElem, iGaus);

                DMat = std::get<Matd2x2>(phaseTrap->CalcDMatx(phi_j, T));

                const Matd2x3& dummyBMat = accessVec(BMat, iElem, iGaus); // derivative matrix for the given gauss point.
                const RowVecd3& N_i = accessVec(shapeFunc, iGaus);
                const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  // Volume of the current int-pt 

                DB = DMat*dummyBMat;

                // Int pt flux
                auto& flux = accessVec(elFlux, iElem, iGaus);
                flux = DB * (
                    -(dummyCon) 
                    + (c_zeta_j  * IntPtCon * dummyElNod_phi_j)
                    + (c_zeta_ii * IntPtCon * dummyElNod_gPhi_ii)
                    + (c_zeta_ij * IntPtCon * dummyElNod_gPhi_ij)
                    + (c_zeta_jj * IntPtCon * dummyElNod_gPhi_jj)
                );

                const auto& NodeConn = elemNodeConn.at(iElem);

                int outIdx = accessVec(elemIDs, iElem) * nElGauss + iGaus;
                ref_intPtFlux[outIdx][0] = flux[0];
                ref_intPtFlux[outIdx][1] = flux[1];
                ref_intPtFlux[outIdx][2] = 0.0;

                for(int iNode=0; iNode<nElNodes; iNode++){
                    double weight = N_i[iNode] * dummydVol;
                    int gNod = NodeConn[iNode]; // Vector lookups using index are much faster
                    
                    ref_nodFlux[gNod] += flux * weight;
                    nodCount.at(gNod) += weight;
                }
            }
        }
    }

    // // TODO: For debug!
    // cout << elFlux.at(0).at(0) << "\n\n";
}

void Tri3TH::CalcFsrc(const double conB, BaseTrapping* mat, double* FsrcBuffer, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){
    
}

void Tri3TH::CalcElCon(const double* globalBuffer){

    ColVecd3 nodCon;

    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElConDofs; iDof++){
            nodCon(iDof) = globalBuffer[accessVec(elemConDof, iElem, iDof)];
        }

        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            const RowVecd3& N_i = accessVec(shapeFunc, iGaus);
            accessVec(elCon, iElem, iGaus) = N_i*nodCon;

        }
    }
}

void Tri3TH::CalcPrecipitaitonCon(BaseTrapping* mat, const double T){

    TrapPhase* phaseDriftMat = dynamic_cast<TrapPhase*>(mat);

    double Xmax = phaseDriftMat->get_Xmax();
    double kf0 = phaseDriftMat->get_kf0();
    double Qf = phaseDriftMat->get_Qf();
    double kb0 = phaseDriftMat->get_kb0();
    double Qb = phaseDriftMat->get_Qb();
    double mask;
    double conTrap; double conTrapOld;
    double conL; double deltsConTrap;

    double kf = kf0 * std::exp(-Qf / (R * T));
    double kb = kb0 * std::exp(-Qb / (R * T));

    const double epsilon = 0.1; 

    for(int iElem=0; iElem<nElements; iElem++){
        
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            double phi = accessVec(el_phi_j, iElem, iGaus);
            
            // Mask is 1.0 ONLY if phi is practically 0.0
            double mask = (phi < epsilon) ? 1.0 : 0.0;

            conTrapOld = accessVec(elConTrap, iElem, iGaus);
            conL = accessVec(elCon, iElem, iGaus);

            conTrap = (conTrapOld + dt * mask * kf * conL * Xmax) / (1.0 + dt * mask * (kf * conL + kb));

            accessVec(elConTrap, iElem, iGaus) = conTrap;
            accessVec(elDeltaTrap, iElem, iGaus) = conTrap - conTrapOld;

        }
    }
}

void Tri3TH::CalcTrapRHS(double* FtrapBuffer){

    ColVecd3 dummyFtrap;

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            const RowVecd3& N_i = accessVec(shapeFunc, iGaus);
            const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  

            dummyFtrap = N_i.transpose()*accessVec(elDeltaTrap, iElem, iGaus)*dummydVol;

            for(int iNod2 = 0; iNod2 < nElConDofs; ++iNod2) {
                int gNod = elemNodeConn[iElem][iNod2];
                FtrapBuffer[gNod] += dummyFtrap[iNod2];
            }
        }
    }
}

void Tri3TH::CalcNodTrapCon(vector<double>& nodTrapCon, vector<double>& nodCount){

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            const double& conTrap = accessVec(elConTrap, iElem, iGaus);
            const double& dummydVol = accessVec(intPtVol, iElem, iGaus);  
            const RowVecd3& N_i = accessVec(shapeFunc, iGaus);

            // Nodal values
            int iNode = 0;
            for(auto iNod2=accessVec(elemNodeConn, iElem).begin(); iNod2!=accessVec(elemNodeConn, iElem).end(); iNod2++){

                double weight = N_i[iNode] * dummydVol;
                accessVec(nodTrapCon, *iNod2) += conTrap * weight;;
                nodCount.at(*iNod2) +=  weight;;
                iNode += 1;

            }
        }
    }
}
