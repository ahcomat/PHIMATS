#include<iostream>
#include<algorithm>

#include"FiniteElements/Mechanics/Quad4Axi.h"
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

Quad4Axi::Quad4Axi(H5IO &H5File_in, H5IO &H5File_mesh, Nodes &Nodes, int iSet, string matModel, Logger& logger)
    : BaseElemMech(2, 4, 4, 2, 4, 8, matModel, logger){ // nElDim, nElNodes, nElGauss, dispDofs, nElStres, nElDispDofs

    if (materialModel != "Elastic" && materialModel != "ElastoPlastic") {
        
        throw std::invalid_argument("Invalid material model: < " + materialModel + " >\nAllowed models are < Elastic, ElastoPlastic >\n");
    }

    InitShapeFunc();
    ReadElementsData(H5File_in, H5File_mesh, iSet);
    InitializeElements(Nodes);
}   

Quad4Axi::~Quad4Axi(){}

void Quad4Axi::InitShapeFunc(){

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
        accessVec(shapeFunc, i) =  CalcShapeFunc(accessVec(gaussPts, i, 0), accessVec(gaussPts, i, 1));
        accessVec(shapeFuncDeriv, i) =  CalcShapeFuncDeriv(accessVec(gaussPts, i, 0), accessVec(gaussPts, i, 1));
    }
}

RowVecd4 Quad4Axi:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd4 shape;

    shape(0) = (1.0-eta-xi+xi*eta)*0.25;
    shape(1) = (1.0-eta+xi-xi*eta)*0.25;
    shape(2) = (1.0+eta+xi+xi*eta)*0.25;
    shape(3) = (1.0+eta-xi-xi*eta)*0.25;

    return shape;
}

Matd2x4 Quad4Axi::CalcShapeFuncDeriv(double xi, double eta){

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

void Quad4Axi::InitializeElements(Nodes &Nodes){

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
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGaus), shapeFunc.at(iGaus), wts.at(iGaus), dummyIntVol.at(iGaus), BMat.at(iElem).at(iGaus), BuMat.at(iElem).at(iGaus));
        }
        gaussPtCart.at(iElem) = dummyElemGauss;
        intPtVol.at(iElem) = dummyIntVol;
    }
}

RowVecd2 Quad4Axi::getGaussCart(RowVecd4& sFunc, Matd4x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Quad4Axi::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const RowVecd4& shFunc, 
                             const double& wt, double& intVol, Matd2x4& cartDeriv, Matd4x8& strainMat) {
    
    // Jacobian J = dN/dxi * nodes
    Matd2x2 jacMat = sFuncDeriv * elNodCoord;
    
    // Radius at the current Gauss point: r = sum(N_i * r_i)
    // Assuming elNodCoord column 0 is 'r' and column 1 is 'z'
    double radius = shFunc.dot(elNodCoord.col(0));

    // Volume for axisymmetry: dV = 2 * PI * r * det(J) * w
    intVol = 2.0 * M_PI * radius * jacMat.determinant() * wt;

    // Cartesian derivatives dN/dx = J^-1 * dN/dxi
    cartDeriv = jacMat.inverse() * sFuncDeriv;

    // Build BuMat (4 x 8)
    strainMat.setZero();
    for(int iNod=0; iNod<nElNodes; iNod++) {
        // 0: Radial strain (eps_rr = du_r/dr)
        strainMat(0, 2*iNod)   = cartDeriv(0, iNod);
        
        // 1: Axial strain (eps_zz = du_z/dz)
        strainMat(1, 2*iNod+1) = cartDeriv(1, iNod);
        
        // 2: Hoop strain (eps_theta = u_r/r)
        strainMat(2, 2*iNod)   = shFunc(iNod) / radius;
        
        // 3: Shear strain (gamma_rz = du_r/dz + du_z/dr)
        strainMat(3, 2*iNod)   = cartDeriv(1, iNod);
        strainMat(3, 2*iNod+1) = cartDeriv(0, iNod);
    }
}

void Quad4Axi::CalcElemStiffMatx(T_DMatx CMatx) {

    // Lift the Material Matrix once at the start
    // Standardizing on Matd4x4 for Axisymmetry
    const Matd4x4& DMat = std::get<Matd4x4>(CMatx);

    for(int iElem = 0; iElem < nElements; iElem++) {

        // Lift the reference to the element stiffness matrix to avoid repeated vector lookups
        Matd8x8& Ke = accessVec(elStiffMatx, iElem);
        Ke.setZero(); 

        // Integration over Gauss points
        for (int iGaus = 0; iGaus < nElGauss; iGaus++) {

            // Lift references to B-matrix and Volume
            const Matd4x8& B = accessVec(BuMat, iElem, iGaus);
            const double dV  = accessVec(intPtVol, iElem, iGaus);

            /* OPTIMIZATION: 
               Instead of (B^T * D * B), we calculate (D * B) first.
               D(4x4) * B(4x8) is much cheaper than B^T(8x4) * D(4x4).
               Eigen's expression templates will optimize the 4x4 * 4x8 product 
               extremely well, often using SIMD instructions.
            */
            const Matd4x8 DB = DMat * B;

            Ke.noalias() += B.transpose() * DB * dV;
        }  
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elStiffMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n";
}

void Quad4Axi::CalcStres(T_DMatx CMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<double>& nodCount){

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
            elStres.at(iElem).at(iGaus) = std::get<Matd4x4>(CMatx)*elStran.at(iElem).at(iGaus);

            dummyForc = BuMat.at(iElem).at(iGaus).transpose()*elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);

            for(int pom=0; pom<nElDispDofs; pom++){
                Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
            }

            // Nodal values
            iNode = 0;
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd4>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd4>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += 1;

                // std::get<std::vector<ColVecd4>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                // std::get<std::vector<ColVecd4>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                // nodCount.at(*iNod2) += 1;
            }
        }
    }

    // // TODO: For debug!
    // cout << elStran.at(0).at(0) << "\n\n";
}

void Quad4Axi::CalcElDStran(const double* globalBuffer){

    ColVecd8 dummyDisp; // for element nodal displacement.

    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElDispDofs; iDof++){
            dummyDisp(iDof) = globalBuffer[accessVec(elemDispDof, iElem, iDof)];
        }

        for(int iGaus=0; iGaus<nElGauss; iGaus++){
            accessVec(elDStran, iElem, iGaus) = accessVec(BuMat, iElem, iGaus)*dummyDisp;
        }
    }
}

void Quad4Axi::CalcElStran(const double* globalBuffer){
    
    ColVecd8 dummyDisp; // for element nodal displacement.

    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElDispDofs; iDof++){
            dummyDisp(iDof) = globalBuffer[accessVec(elemDispDof, iElem, iDof)];
        }

        for(int iGaus=0; iGaus<nElGauss; iGaus++){
            accessVec(elStran, iElem, iGaus) = accessVec(BuMat, iElem, iGaus)*dummyDisp;
        }
    }
}

void Quad4Axi::CalcNodVals( T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<double>& nodStres_h, vector<double>& nodRho, vector<double>& nodCount){

    try {
        if (elStran_e.data() == nullptr){
            throw runtime_error("Plasicity strain container vectors were not allocated.\n    Please add the keyword argument < Elastoplastic > in the element constructor.\n");
        }
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        logger.log("Exception caught in CalcNodVals::CalcNodVals\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("    Critical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);

    }

    // LIFT: Extract references once before the hot loops
    // This removes runtime type-checking from the nested loops
    vector<ColVecd4>& v_nodStran   = std::get<std::vector<ColVecd4>>(nodStran);
    vector<ColVecd4>& v_nodStran_e = std::get<std::vector<ColVecd4>>(nodStran_e);
    vector<ColVecd4>& v_nodStran_p = std::get<std::vector<ColVecd4>>(nodStran_p);
    vector<ColVecd4>& v_nodStres   = std::get<std::vector<ColVecd4>>(nodStres);

    // Optimized Hot Loops
    for (int iElem = 0; iElem < nElements; iElem++) {
        for (int iGaus = 0; iGaus < nElGauss; iGaus++) {
            
            // Pre-fetch Gauss point values to avoid redundant accessVec calls inside the node loop
            const double vol = accessVec(intPtVol, iElem, iGaus);
            const auto& es   = accessVec(elStres, iElem, iGaus);
            const auto& en   = accessVec(elStran, iElem, iGaus);
            const auto& ene  = accessVec(elStran_e, iElem, iGaus);
            const auto& enp  = accessVec(elStran_p, iElem, iGaus);
            
            const double s_eq = accessVec(elStres_eq, iElem, iGaus);
            const double n_eq = accessVec(elStran_eq, iElem, iGaus);
            const double s_h  = accessVec(elStres_h, iElem, iGaus);
            const double rho  = accessVec(elRho, iElem, iGaus);

            for (auto nodeIdx : elemNodeConn.at(iElem)) {
                // Now using direct vector indexing and pre-fetched values
                v_nodStran[nodeIdx]   += en * vol;
                v_nodStran_e[nodeIdx] += ene * vol;
                v_nodStran_p[nodeIdx] += enp * vol;
                v_nodStres[nodeIdx]   += es * vol;

                nodStran_eq[nodeIdx] += n_eq * vol;
                nodStres_eq[nodeIdx] += s_eq * vol;
                nodStres_h[nodeIdx]  += s_h * vol;
                nodRho[nodeIdx]      += rho * vol;
                nodCount[nodeIdx]    += vol;
            }
        }
    }
}

void Quad4Axi::CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep){

    IsoHard* plasticMat = dynamic_cast<IsoHard*>(mat);

    // LIFT: Resolve the variant types once for the whole element set
    Matd4x4& Cep = std::get<Matd4x4>(plasticMat->getCMatx_ep());   // Elastoplastic stiffness matrix
    const Matd4x4& Ce = std::get<Matd4x4>(plasticMat->getCMatx()); // Elastic stiffness
    
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
                
                plasticMat->ReturnMappingAxi(accessVec(elDStran, iElem, iGaus),
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
                                        Ce, Cep);

            }
        }

    } else {  // Update the element stiffness matrix

        double dummydVol;   // dummy for int-pt volume.

        for(int iElem=0; iElem<nElements; iElem++){

            elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

            for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                plasticMat->ReturnMappingAxi(accessVec(elDStran, iElem, iGaus),
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
                                        Ce, Cep);

                const Matd4x8& dummyBu = accessVec(BuMat, iElem, iGaus); // Strain matrix for the given gauss point.
                dummydVol = accessVec(intPtVol, iElem, iGaus);  // Volume of the current integration point 

                // [B_kl]^T D_kk B_kl
                accessVec(elStiffMatx, iElem).noalias() +=  dummyBu.transpose()*Cep*dummyBu*dummydVol;
            }
        }
    }
}

void Quad4Axi::CalcRetrunMapping_PFF(BaseMechanics* mat, const bool& updateStiffMat, int iStep, const std::vector<std::vector<double>>* gPhi_d_ptr){

    IsoHard* plasticMat = dynamic_cast<IsoHard*>(mat);

    // LIFT: Resolve the variant types once for the whole element set
    Matd4x4& Cep = std::get<Matd4x4>(plasticMat->getCMatx_ep());   // Elastoplastic stiffness matrix
    const Matd4x4& Ce = std::get<Matd4x4>(plasticMat->getCMatx()); // Elastic stiffness

    if (!updateStiffMat){

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGaus=0; iGaus<nElGauss; iGaus++){
            
            plasticMat->ReturnMappingAxi_PFF(accessVec(elDStran, iElem, iGaus),
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
                                        accessVec(el_wp, iElem, iGaus),
                                        Ce, Cep);

        }
    }

    } else {  // Update the element stiffness matrix

        double dummydVol;   // dummy for int-pt volume.

        for(int iElem=0; iElem<nElements; iElem++){

            elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

            for(int iGaus=0; iGaus<nElGauss; iGaus++){
                
                plasticMat->ReturnMappingAxi_PFF(accessVec(elDStran, iElem, iGaus),
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
                                            accessVec(el_wp, iElem, iGaus),
                                            Ce, Cep);

                const Matd4x8& dummyBu = accessVec(BuMat, iElem, iGaus); // Strain matrix for the given gauss point.
                dummydVol = accessVec(intPtVol, iElem, iGaus);  // Volume of the current integration point 

                // [B_kl]^T D_kk B_kl
                elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*Cep*dummyBu*dummydVol;
            }
        }
    }

}

void Quad4Axi::CalcFint(double* Fint){

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

void Quad4Axi::getNew(){

    elStran_e_old = elStran_e;
    elStran_p_old = elStran_p;
    elStran_eq_old = elStran_eq;
    el_wp_old = el_wp;

}

