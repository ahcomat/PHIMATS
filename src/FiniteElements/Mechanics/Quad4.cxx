#include<iostream>
#include<algorithm>

#include"FiniteElements/Mechanics/Quad4.h"
#include"Materials/Mechanics/IsoHard.h"

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

Quad4::Quad4(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, string matModel, Logger& logger)
    : BaseElemMech(2, 4, 4, 2, 3, 8, matModel, logger){ // nElDim, nElNodes, nElGauss, dispDofs, nElStres, nElDispDofs

    if (materialModel != "Elastic" && materialModel != "ElastoPlastic") {
        
        throw std::invalid_argument("Invalid material model: < " + materialModel + " >\nAllowed models are < Elastic, ElastoPlastic >\n");
    }

    InitShapeFunc();
    ReadElementsData(H5File_in, H5File_mesh, iSet);
    InitializeElements(Nodes);
}   

Quad4::~Quad4(){}

void Quad4::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    vector<double> ip = {-0.57735027, 0.57735027};
    vector<double> dummy(nElDim);
    for(std::size_t iGaus=0; iGaus<ip.size(); iGaus++){
        for(std::size_t jGauss=0; jGauss<ip.size(); jGauss++){
            dummy.at(0) = ip.at(iGaus);
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

RowVecd4 Quad4:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd4 shape;

    shape(0) = (1.0-eta-xi+xi*eta)*0.25;
    shape(1) = (1.0-eta+xi-xi*eta)*0.25;
    shape(2) = (1.0+eta+xi+xi*eta)*0.25;
    shape(3) = (1.0+eta-xi-xi*eta)*0.25;

    return shape;
}

Matd2x4 Quad4::CalcShapeFuncDeriv(double xi, double eta){

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

void Quad4::InitializeElements(Nodes &Nodes){

    // Initialize the vector containing each element stiffness matrix.
    elStiffMatx.resize(nElements); 

    // Move after allocation for `elStiffMatx`.
    elStiffMatxVariant = &elStiffMatx;

    // Initialize the storages for int-pt stresses/strains
    elStres.resize(nElements); elStran.resize(nElements); elDStran.resize(nElements);

    // For elastitcy total strain is elastic strain!
    elStrain_e_Variant = &elStran;

    if (materialModel=="ElastoPlastic"){
        elStran_e.resize(nElements); elStran_e_old.resize(nElements);   // Elastic strain tensors
        elStran_p.resize(nElements); elStran_p_old.resize(nElements);   // Plastic strain tensors
        elStran_eq.resize(nElements); elStran_eq_old.resize(nElements); // Equivalent platsic strain
        elStres_eq.resize(nElements); // Equivalent (von Mises) stress
        elStres_h.resize(nElements);  // Hydrostatic stress
        elRho.resize(nElements);      // Norm dislocation density
        el_wp.resize(nElements); el_wp_old.resize(nElements);     // Plastic work density
        

        elStrain_e_Variant = &elStran_e;
        el_wp_ptr = &el_wp;
    }     

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd4x2 dummyElNodCoord; // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.
    BuMat.resize(nElements);  // Initialize the size of BuMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    // Loop through elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStres.at(iElem).resize(nElGauss);
        elStran.at(iElem).resize(nElGauss);
        elDStran.at(iElem).resize(nElGauss);

        // Initilize to zeros.
        for (int iGaus=0; iGaus<nElGauss; iGaus++){
            elStres.at(iElem).at(iGaus).setZero();
            elStran.at(iElem).at(iGaus).setZero();
            elDStran.at(iElem).at(iGaus).setZero();
        }

        if (materialModel=="ElastoPlastic"){

            elStran_e.at(iElem).resize(nElGauss);
            elStran_p.at(iElem).resize(nElGauss); 
            elStres_eq.at(iElem).resize(nElGauss);
            elStran_eq.at(iElem).resize(nElGauss);
            elStres_h.at(iElem).resize(nElGauss);
            elRho.at(iElem).resize(nElGauss);
            accessVec(el_wp, iElem).resize(nElGauss);

            elStran_e_old.at(iElem).resize(nElGauss);
            elStran_p_old.at(iElem).resize(nElGauss); 
            elStran_eq_old.at(iElem).resize(nElGauss);
            accessVec(el_wp_old, iElem).resize(nElGauss);

            // Initilize to zeros.
            for (int iGaus=0; iGaus<nElGauss; iGaus++){

                elStran_e.at(iElem).at(iGaus).setZero();
                elStran_p.at(iElem).at(iGaus).setZero();
                elStran_eq.at(iElem).at(iGaus) = 0;
                elStres_eq.at(iElem).at(iGaus) = 0;
                elStres_h.at(iElem).at(iGaus) = 0;
                elRho.at(iElem).at(iGaus) = 0;
                accessVec(el_wp, iElem, iGaus) = 0;

                elStran_e_old.at(iElem).at(iGaus).setZero();
                elStran_p_old.at(iElem).at(iGaus).setZero();
                elStran_eq_old.at(iElem).at(iGaus) = 0;
            }
        }

        BMat.at(iElem).resize(nElGauss);
        BuMat.at(iElem).resize(nElGauss);

        // Loop through nodes to get coordinates.
        for(int iNod=0; iNod<nElNodes; iNod++){

            dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(0);
            dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(1);
        }
        elemNodCoord.at(iElem) = dummyElNodCoord;

        // Loop through integration points.
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            // For safety
            BMat.at(iElem).at(iGaus).setZero();
            BuMat.at(iElem).at(iGaus).setZero();
        
            // Cart coord of iGaus point.
            dummyElemGauss.at(iGaus) = getGaussCart(shapeFunc.at(iGaus), dummyElNodCoord);
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGaus), wts.at(iGaus), dummyIntVol.at(iGaus), BMat.at(iElem).at(iGaus), BuMat.at(iElem).at(iGaus));
        }
        gaussPtCart.at(iElem) = dummyElemGauss;
        intPtVol.at(iElem) = dummyIntVol;
    }
}

RowVecd2 Quad4::getGaussCart(RowVecd4& sFunc, Matd4x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Quad4::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv, Matd3x8& strainMat){

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

    // Strain matrix of the current integration point du_kl.
    for(int iNod=0; iNod<nElNodes; iNod++){

        strainMat(0,2*iNod)   = cartDeriv(0,iNod);
        strainMat(0,2*iNod+1) = 0;
        strainMat(1,2*iNod)   = 0;
        strainMat(1,2*iNod+1) = cartDeriv(1,iNod);
        strainMat(2,2*iNod)   = cartDeriv(1,iNod);
        strainMat(2,2*iNod+1) = cartDeriv(0,iNod);
    }
}

void Quad4::CalcElemStiffMatx(T_DMatx CMatx){

    Matd3x3 DMat = std::get<Matd3x3>(CMatx);

    double dummydVol;   // dummy for int-pt volume.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

        // Integration over all Gauss points.
        for (int iGaus=0; iGaus<nElGauss; iGaus++){

            const Matd3x8& dummyBu = BuMat.at(iElem).at(iGaus); // Strain matrix for the given gauss point.
            dummydVol = intPtVol.at(iElem).at(iGaus);  // Volume of the current integration point 

            // [B_kl]^T D_kk B_kl
            elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*DMat*dummyBu*dummydVol;
        }  
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elStiffMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n";
}

void Quad4::CalcStres(T_DMatx CMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<double>& nodCount){

    ColVecd8 dummyDisp; // for element nodal displacement.
    ColVecd8 dummyForc; // for element nodal internal force.
    int iNode;  // counter for the number of nodes.

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElDispDofs; iDof++){
            dummyDisp(iDof) = globalBuffer[elemDispDof.at(iElem).at(iDof)];
        }

        // Gauss points
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            // Int pt values
            elStran.at(iElem).at(iGaus) = BuMat.at(iElem).at(iGaus)*dummyDisp;
            elStres.at(iElem).at(iGaus) = std::get<Matd3x3>(CMatx)*elStran.at(iElem).at(iGaus);

            dummyForc = BuMat.at(iElem).at(iGaus).transpose()*elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);

            // Sometimes it throws segmentation fault, have no idea why (⊙_⊙)？
            for(int pom=0; pom<nElDispDofs; pom++){
                Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
            }

            // Nodal values
            iNode = 0;
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd3>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd3>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += 1;

                // std::get<std::vector<ColVecd3>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                // std::get<std::vector<ColVecd3>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                // nodCount.at(*iNod2) += 1;
            }
        }
    }

    // // TODO: For debug!
    // cout << elStran.at(0).at(0) << "\n\n";
}

void Quad4::CalcElDStran(const double* globalBuffer){

    ColVecd8 dummyDisp; // for element nodal displacement.

    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElDispDofs; iDof++){
            dummyDisp(iDof) = globalBuffer[elemDispDof.at(iElem).at(iDof)];
        }

        for(int iGaus=0; iGaus<nElGauss; iGaus++){
            elDStran.at(iElem).at(iGaus) = BuMat.at(iElem).at(iGaus)*dummyDisp;
        }
    }
}

void Quad4::CalcElStran(const double* globalBuffer){
    
    ColVecd8 dummyDisp; // for element nodal displacement.

    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElDispDofs; iDof++){
            dummyDisp(iDof) = globalBuffer[elemDispDof.at(iElem).at(iDof)];
        }

        for(int iGaus=0; iGaus<nElGauss; iGaus++){
            elStran.at(iElem).at(iGaus) = BuMat.at(iElem).at(iGaus)*dummyDisp;
        }
    }
}

void Quad4::CalcNodVals( T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<double>& nodStres_h, vector<double>& nodRho, vector<double>& nodCount){

    try {
        if (elStran_e.data() == nullptr){
            throw runtime_error("Plasicity strain container vectors were not allocated.\n\n       Please add the keyword argument < Elastoplastic > in the element constructor.\n");
        }
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
    }

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){
        // Gauss points
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            // Nodal values
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd3>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd3>>(nodStran_e).at(*iNod2) += elStran_e.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd3>>(nodStran_p).at(*iNod2) += elStran_p.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd3>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                nodStran_eq.at(*iNod2) += elStran_eq.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                nodStres_eq.at(*iNod2) += elStres_eq.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                nodStres_h.at(*iNod2) += elStres_h.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                nodRho.at(*iNod2) += elRho.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += intPtVol.at(iElem).at(iGaus);
            }
        }
    }
}

void Quad4::CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep){

    IsoHard* plasticMat = dynamic_cast<IsoHard*>(mat);

    try {
        if (elStran_e.data() == nullptr){
            throw runtime_error("Plasicity strain container vectors were not allocated.\n\n       Please add the keyword argument < Elastoplastic > in the element constructor.\n");
        }
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
    }

    if (!updateStiffMat){

        for(int iElem=0; iElem<nElements; iElem++){
            for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                plasticMat->ReturnMapping2D(elDStran.at(iElem).at(iGaus),
                                            elStres.at(iElem).at(iGaus),
                                            elStran_e.at(iElem).at(iGaus),
                                            elStran_p.at(iElem).at(iGaus),
                                            elStran_eq.at(iElem).at(iGaus), 
                                            elStres_eq.at(iElem).at(iGaus),
                                            elStres_h.at(iElem).at(iGaus),
                                            elRho.at(iElem).at(iGaus),
                                            elStran_e_old.at(iElem).at(iGaus),
                                            elStran_p_old.at(iElem).at(iGaus),
                                            elStran_eq_old.at(iElem).at(iGaus), iStep);

            }
        }

    } else {  // Update the element stiffness matrix

        double dummydVol;   // dummy for int-pt volume.

        for(int iElem=0; iElem<nElements; iElem++){

            elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

            for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                plasticMat->ReturnMapping2D(elDStran.at(iElem).at(iGaus),
                                            elStres.at(iElem).at(iGaus),
                                            elStran_e.at(iElem).at(iGaus),
                                            elStran_p.at(iElem).at(iGaus),
                                            elStran_eq.at(iElem).at(iGaus), 
                                            elStres_eq.at(iElem).at(iGaus),
                                            elStres_h.at(iElem).at(iGaus),
                                            elRho.at(iElem).at(iGaus),
                                            elStran_e_old.at(iElem).at(iGaus),
                                            elStran_p_old.at(iElem).at(iGaus),
                                            elStran_eq_old.at(iElem).at(iGaus), iStep);

                const Matd3x8& dummyBu = BuMat.at(iElem).at(iGaus); // Strain matrix for the given gauss point.
                dummydVol = intPtVol.at(iElem).at(iGaus);  // Volume of the current integration point 

                // [B_kl]^T D_kk B_kl
                elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*std::get<Matd3x3>(plasticMat->getCMatx())*dummyBu*dummydVol;
            }
        }
    }
}

void Quad4::CalcRetrunMapping_PFF(BaseMechanics* mat, const bool& updateStiffMat, int iStep, const std::vector<std::vector<double>>* gPhi_d_ptr){

    IsoHard* plasticMat = dynamic_cast<IsoHard*>(mat);

    if (!updateStiffMat){

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGaus=0; iGaus<nElGauss; iGaus++){
            
            plasticMat->ReturnMapping2D_PFF(accessVec(elDStran, iElem, iGaus),
                                        accessVec(elStres, iElem, iGaus),
                                        accessVec(elStran_e, iElem, iGaus),
                                        accessVec(elStran_p, iElem, iGaus),
                                        accessVec(elStran_eq, iElem, iGaus),
                                        accessVec(elStres_eq, iElem, iGaus),
                                        accessVec(elStres_h, iElem, iGaus),
                                        accessVec(elRho, iElem, iGaus),
                                        accessVec(elStran_e_old, iElem, iGaus),
                                        accessVec(elStran_p_old, iElem, iGaus),
                                        accessVec(elStran_eq_old, iElem, iGaus), iStep,
                                        accessVec(*gPhi_d_ptr, iElem, iGaus),
                                        accessVec(el_wp_old, iElem, iGaus),
                                        accessVec(el_wp, iElem, iGaus));

        }
    }

    } else {  // Update the element stiffness matrix

        double dummydVol;   // dummy for int-pt volume.

        for(int iElem=0; iElem<nElements; iElem++){

            elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

            for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                plasticMat->ReturnMapping2D_PFF(accessVec(elDStran, iElem, iGaus),
                                            accessVec(elStres, iElem, iGaus),
                                            accessVec(elStran_e, iElem, iGaus),
                                            accessVec(elStran_p, iElem, iGaus),
                                            accessVec(elStran_eq, iElem, iGaus),
                                            accessVec(elStres_eq, iElem, iGaus),
                                            accessVec(elStres_h, iElem, iGaus),
                                            accessVec(elRho, iElem, iGaus),
                                            accessVec(elStran_e_old, iElem, iGaus),
                                            accessVec(elStran_p_old, iElem, iGaus),
                                            accessVec(elStran_eq_old, iElem, iGaus), iStep,
                                            accessVec(*gPhi_d_ptr, iElem, iGaus),
                                            accessVec(el_wp_old, iElem, iGaus),
                                            accessVec(el_wp, iElem, iGaus));

                const Matd3x8& dummyBu = BuMat.at(iElem).at(iGaus); // Strain matrix for the given gauss point.
                dummydVol = intPtVol.at(iElem).at(iGaus);  // Volume of the current integration point 

                // [B_kl]^T D_kk B_kl
                elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*std::get<Matd3x3>(plasticMat->getCMatx())*dummyBu*dummydVol;
            }
        }
    }

}

void Quad4::CalcFint(double* Fint){

    ColVecd8 dummyForc; // for element nodal internal force.

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            dummyForc = BuMat.at(iElem).at(iGaus).transpose()*elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);

            for(int pom=0; pom<nElDispDofs; pom++){
                Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
            }
        }
    }
}

void Quad4::getNew(){

    elStran_e_old = elStran_e;
    elStran_p_old = elStran_p;
    elStran_eq_old = elStran_eq;
    el_wp_old = el_wp;

}

