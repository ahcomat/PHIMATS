#include<iostream>
#include<algorithm>

#include "FiniteElements/Trapping/Tri6TH.h"
#include "Materials/Trapping/TrapGB.h"
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

Tri6TH::Tri6TH(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, Logger& logger)
    : BaseElemTrap(2, 6, 6, 3, logger){ // nElDim, nElNodes, nElConDofs, nElGauss 

    InitShapeFunc();
    ReadElementsData(H5File_in, H5File_in, iSet);
    InitializeElements(Nodes, H5File_in);
}   

Tri6TH::~Tri6TH(){

    // Exit message
    cout << "Tri6TH elements exited correctly" << "\n";
}

void Tri6TH::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    gaussPts.push_back({1.0/6.0, 1.0/6.0});
    gaussPts.push_back({2.0/3.0, 1.0/6.0});
    gaussPts.push_back({1.0/6.0, 2.0/3.0});

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nElGauss);
    shapeFuncDeriv.resize(nElGauss);

    for(int i=0; i<nElGauss; i++){
        shapeFunc.at(i) =  CalcShapeFunc(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
        shapeFuncDeriv.at(i) = CalcShapeFuncDeriv(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
    }
}

RowVecd6 Tri6TH:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd6 shape;

    shape(0) = (1.0 - xi - eta) * (2.0 * (1.0 - xi - eta) - 1.0); 
    shape(1) = xi * (2.0 * xi - 1.0);                             
    shape(2) = eta * (2.0 * eta - 1.0);                           
    shape(3) = 4.0 * xi * (1.0 - xi - eta);                       
    shape(4) = 4.0 * xi * eta;                                    
    shape(5) = 4.0 * eta * (1.0 - xi - eta);                      

    return shape;
}

Matd2x6 Tri6TH::CalcShapeFuncDeriv(double xi, double eta){

    // dN_ji
    Matd2x6 shapeDeriv;

    shapeDeriv(0, 0) = -3.0 + 4.0 * xi + 4.0 * eta;   
    shapeDeriv(0, 1) = 4.0 * xi - 1.0;               
    shapeDeriv(0, 2) = 0.0;                          
    shapeDeriv(0, 3) = 4.0 - 8.0 * xi - 4.0 * eta;   
    shapeDeriv(0, 4) = 4.0 * eta;                    
    shapeDeriv(0, 5) = -4.0 * eta;                   

    shapeDeriv(1, 0) = -3.0 + 4.0 * xi + 4.0 * eta;   
    shapeDeriv(1, 1) = 0.0;                          
    shapeDeriv(1, 2) = 4.0 * eta - 1.0;              
    shapeDeriv(1, 3) = -4.0 * xi;                    
    shapeDeriv(1, 4) = 4.0 * xi;                     
    shapeDeriv(1, 5) = 4.0 - 4.0 * xi - 8.0 * eta;    
             
    return shapeDeriv;
}

void Tri6TH::InitializeElements(Nodes &Nodes, H5IO &H5File_in){

    // Initialize the storages
    elFlux.resize(nElements);
    elStiffMatx.resize(nElements);
    elCapMatx.resize(nElements);

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd6x2 dummyElNodCoord;   // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    if (Trapping=="GBTrapping"){       // GB

        ColVecd6 dummyElNod_gPhi;  // For element nodal values of phi.

        // Initialize the storages for int-pt gPhi
        el_gPhi.resize(nElements);
        // Read nodal values of gPhi
        nod_gPhi.resize(nNodes);
        string dsetName = "gPhi";
        H5File_in.ReadField1D(dsetName, nod_gPhi);

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

    } else if (Trapping=="2PhaseTrapping") {       // Phase 

        ColVecd6 dummyElNod_phi_j; ColVecd6 dummyElNod_gPhi_jj;
        ColVecd6 dummyElNod_gPhi_ii; ColVecd6 dummyElNod_gPhi_ij;

        // Initialize the storages for int-pt traps
        el_phi_j.resize(nElements); el_gPhi_ii.resize(nElements);
        el_gPhi_jj.resize(nElements); el_gPhi_ij.resize(nElements);

        // Read nodal values of traps
        nod_phi_j.resize(nNodes);
        string dsetName = "phi_j";
        H5File_in.ReadField1D(dsetName, nod_phi_j);

        nod_gPhi_ii.resize(nNodes);
        dsetName = "gPhi_ff";
        H5File_in.ReadField1D(dsetName, nod_gPhi_ii);

        nod_gPhi_ij.resize(nNodes);
        dsetName = "gPhi_fM";
        H5File_in.ReadField1D(dsetName, nod_gPhi_ij);

        nod_gPhi_jj.resize(nNodes);
        dsetName = "gPhi_MM";
        H5File_in.ReadField1D(dsetName, nod_gPhi_jj);

        // Loop through elements.
        for(int iElem=0; iElem<nElements; iElem++){

            elFlux.at(iElem).resize(nElGauss);

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

    }  else {

        cerr << "Error trapping flag: " << Trapping << "\n";
        cerr << "Terminating!\n\n";
        exit(10);
    }
}

RowVecd2 Tri6TH::getGaussCart(RowVecd6& sFunc, Matd6x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Tri6TH::CalcCartDeriv(Matd6x2& elNodCoord, Matd2x6& sFuncDeriv, const double& wt, double& intVol, Matd2x6& cartDeriv){

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

void Tri6TH::getInPtCoords(T_nodStres& glIntPtCoords){

    for(int iElem=0; iElem<nElements; iElem++){

        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            std::get<std::vector<ColVecd3>>(glIntPtCoords).at(elemIDs.at(iElem)*nElGauss+iGaus)[0] =  gaussPtCart.at(iElem).at(iGaus)[0];
            std::get<std::vector<ColVecd3>>(glIntPtCoords).at(elemIDs.at(iElem)*nElGauss+iGaus)[1] =  gaussPtCart.at(iElem).at(iGaus)[1];
            std::get<std::vector<ColVecd3>>(glIntPtCoords).at(elemIDs.at(iElem)*nElGauss+iGaus)[2] =  0;

        }
    }
}

// void Tri6TH::CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi){


// }

void Tri6TH::CalcElemStiffMatx(BaseTrapping* mat, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){

    Matd2x2 DMat; 

    vector<Matd6x6> elKDMatx(nElements); 
    vector<Matd6x6> elKTMatx(nElements); 

    double dummydVol;       // dummy for int-pt volume.

    if (Trapping=="GBTrapping"){       // GB

        Matd2x2 TMat; 

        ColVecd6 dummyElNod_gPhi;  // For element nodal values of phi.
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

                const Matd2x6& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
                const RowVecd6& dummyShFunc = shapeFunc.at(iGauss);
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
    } else if (Trapping=="2PhaseTrapping"){        // Phase

        ColVecd6 dummyElNod_phi_j, dummyElNod_gPhi_ij, dummyElNod_gPhi_ii,
        dummyElNod_gPhi_jj;
        double phi_j, gPhi_ii, gPhi_ij, gPhi_jj;              

        vector<int> NodeConn;

        double zeta_j = dynamic_cast<TrapPhase*>(mat)->get_zeta_j();
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

                phi_j = el_phi_j.at(iElem).at(iGauss);  // values of the current int-pt
                gPhi_ii = el_gPhi_ii.at(iElem).at(iGauss);
                gPhi_ij = el_gPhi_ij.at(iElem).at(iGauss);
                gPhi_jj = el_gPhi_jj.at(iElem).at(iGauss);

                DMat = std::get<Matd2x2>(mat->CalcDMatx(phi_j, T));

                const Matd2x6& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
                const RowVecd6& dummyShFunc = shapeFunc.at(iGauss);
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
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elKDMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
    elCapMatxVariant = &elCapMatx;
}

double Tri6TH::CalcAvCon(const double* globalBuffer){

    double IntPtCon, AvCon = 0, TotVol = 0;
    ColVecd6 dummyCon; 

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

void Tri6TH::CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){

    Matd2x2 DMat; 
    double IntPtCon;   // Integration point concnetration
    ColVecd6 dummyCon; // for element nodal concentration.
    int iNode=0;

    if (Trapping=="GBTrapping"){       // GB
        
        Matd2x2 TMat; 
        double gPhi;
        ColVecd6 dummyElNod_gPhi; // for element nodal concentration.

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

                DMat = std::get<Matd2x2>(mat->CalcDMatx(gPhi, T));
                TMat = std::get<Matd2x2>(dynamic_cast<TrapGB*>(mat)->CalcTMatx(gPhi, T));

                // Int pt flux
                elFlux.at(iElem).at(iGaus) = - DMat*BMat.at(iElem).at(iGaus)*dummyCon 
                                            + TMat*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_gPhi;

                // elFlux.at(iElem).at(iGaus) = - DMat*BMat.at(iElem).at(iGaus)*dummyCon;
                // elFlux.at(iElem).at(iGaus) = TMat*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_gPhi;

                std::get<std::vector<ColVecd3>>(intPtFlux).at(elemIDs.at(iElem)*nElGauss+iGaus)[0] =  elFlux.at(iElem).at(iGaus)[0];
                std::get<std::vector<ColVecd3>>(intPtFlux).at(elemIDs.at(iElem)*nElGauss+iGaus)[1] =  elFlux.at(iElem).at(iGaus)[1];
                std::get<std::vector<ColVecd3>>(intPtFlux).at(elemIDs.at(iElem)*nElGauss+iGaus)[2] =  0;

                // // Nodal values
                // for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){
                //     std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += elFlux.at(iElem).at(iGaus);
                //     nodCount.at(*iNod2) += 1;
                // }

                // Nodal values
                iNode = 0;
                for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                    std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += elFlux.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                    nodCount.at(*iNod2) += shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                    iNode += 1;
                }
            }
        }

    } else if (Trapping=="2PhaseTrapping"){        // Phase

        ColVecd6 dummyElNod_phi_j, dummyElNod_gPhi_ij, dummyElNod_gPhi_ii,
        dummyElNod_gPhi_jj;
        double phi_j, gPhi_ii, gPhi_ij, gPhi_jj;              

        double zeta_j = dynamic_cast<TrapPhase*>(mat)->get_zeta_j();
        double zeta_ij = dynamic_cast<TrapPhase*>(mat)->get_zeta_ij();
        double zeta_jj = dynamic_cast<TrapPhase*>(mat)->get_zeta_jj();
        double zeta_ii = dynamic_cast<TrapPhase*>(mat)->get_zeta_ii();

        // Integration point values.
        for(int iElem=0; iElem<nElements; iElem++){

            // Loop through element nodes to get nodal values.
            for(int iDof=0; iDof<nElConDofs; iDof++){
                dummyCon[iDof] = globalBuffer[elemConDof.at(iElem).at(iDof)];
                dummyElNod_phi_j[iDof] = nod_phi_j.at(elemConDof.at(iElem).at(iDof));
                dummyElNod_gPhi_ii[iDof] = nod_gPhi_ii.at(elemConDof.at(iElem).at(iDof));
                dummyElNod_gPhi_ij[iDof] = nod_gPhi_ij.at(elemConDof.at(iElem).at(iDof));
                dummyElNod_gPhi_jj[iDof] = nod_gPhi_jj.at(elemConDof.at(iElem).at(iDof));
            }

            // Gauss points
            for(int iGaus=0; iGaus<nElGauss; iGaus++){

                IntPtCon = shapeFunc.at(iGaus)*dummyCon;

                phi_j = el_phi_j.at(iElem).at(iGaus);  // values of the current int-pt
                gPhi_ii = el_gPhi_ii.at(iElem).at(iGaus);
                gPhi_ij = el_gPhi_ij.at(iElem).at(iGaus);
                gPhi_jj = el_gPhi_jj.at(iElem).at(iGaus);

                DMat = std::get<Matd2x2>(mat->CalcDMatx(phi_j, T));

                // Int pt flux
                elFlux.at(iElem).at(iGaus) = - DMat*BMat.at(iElem).at(iGaus)*dummyCon 
                                             + DMat*zeta_j/(R*T)*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_phi_j
                                             + DMat*zeta_ii/(R*T)*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_gPhi_ii
                                             + DMat*zeta_ij/(R*T)*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_gPhi_ij
                                             + DMat*zeta_jj/(R*T)*IntPtCon*BMat.at(iElem).at(iGaus)*dummyElNod_gPhi_jj;

                std::get<std::vector<ColVecd3>>(intPtFlux).at(elemIDs.at(iElem)*nElGauss+iGaus)[0] =  elFlux.at(iElem).at(iGaus)[0];
                std::get<std::vector<ColVecd3>>(intPtFlux).at(elemIDs.at(iElem)*nElGauss+iGaus)[1] =  elFlux.at(iElem).at(iGaus)[1];
                std::get<std::vector<ColVecd3>>(intPtFlux).at(elemIDs.at(iElem)*nElGauss+iGaus)[2] =  0;

                // Nodal values
                iNode = 0;
                for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                    std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += elFlux.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                    nodCount.at(*iNod2) += shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                    iNode += 1;
                }
            }
        }
    }

    // // TODO: For debug!
    // cout << elFlux.at(0).at(0) << "\n\n";
}

void Tri6TH::CalcFsrc(const double conB, BaseTrapping* mat, double* FsrcBuffer, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr){
    
}

void Tri6TH::CalcElCon(const double* globalBuffer){

}