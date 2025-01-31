#include<iostream>

#include"FiniteElements/Mechanics/Hex8.h"
#include"Materials/Mechanics/IsoHard.h"
#include <iomanip> // For std::setprecision


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

Hex8::Hex8(H5IO &H5File_in, Nodes &Nodes, int iSet, string matModel, Logger& logger)
    : BaseElemMech(3, 8, 8, 3, 6, 24,  matModel, logger){ // nElDim, nElNodes, nElGauss, dispDofs, nElStres, nElDispDofs

    if (materialModel != "Elastic" && materialModel != "ElastoPlastic") {
        
        throw std::invalid_argument("Invalid material model: < " + materialModel + " >\nAllowed models are < Elastic, ElastoPlastic >\n");
    }

    InitShapeFunc();
    ReadElementsData(H5File_in, iSet);
    InitializeElements(Nodes);
}

Hex8::~Hex8(){

    // Exit message
    cout << "Elements exited correctly" << "\n";
}

void Hex8::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    vector<double> ip = {-0.57735027, 0.57735027};
    vector<double> dummy(nElDim);

    for(int iGaus=0; iGaus<2; iGaus++){
        for(int jGaus=0; jGaus<2; jGaus++){
            for(int kGaus=0; kGaus<2; kGaus++){

                dummy.at(0) = ip.at(iGaus);
                dummy.at(1) = ip.at(jGaus);
                dummy.at(2) = ip.at(kGaus);

                gaussPts.push_back(dummy);
            }
        }
    }

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nElGauss);
    shapeFuncDeriv.resize(nElGauss);

    for(int i=0; i<nElGauss; i++){
        shapeFunc.at(i) =  CalcShapeFunc(gaussPts.at(i).at(0), gaussPts.at(i).at(1), gaussPts.at(i).at(2));
        shapeFuncDeriv.at(i) = CalcShapeFuncDeriv(gaussPts.at(i).at(0), gaussPts.at(i).at(1), gaussPts.at(i).at(2));
    }
}

RowVecd8 Hex8::CalcShapeFunc(double xi, double eta, double zeta){

    // N_i
    RowVecd8 shape;

    shape(0) = 0.125*(1.0-xi)*(1.0-eta)*(1.0-zeta);
    shape(1) = 0.125*(1.0+xi)*(1.0-eta)*(1.0-zeta);
    shape(2) = 0.125*(1.0+xi)*(1.0+eta)*(1.0-zeta);
    shape(3) = 0.125*(1.0-xi)*(1.0+eta)*(1.0-zeta);
    shape(4) = 0.125*(1.0-xi)*(1.0-eta)*(1.0+zeta);
    shape(5) = 0.125*(1.0+xi)*(1.0-eta)*(1.0+zeta);
    shape(6) = 0.125*(1.0+xi)*(1.0+eta)*(1.0+zeta);
    shape(7) = 0.125*(1.0-xi)*(1.0+eta)*(1.0+zeta);

    return shape;
}

Matd3x8 Hex8::CalcShapeFuncDeriv(double xi, double eta, double zeta){

    // dN_ji
    Matd3x8 shapeDeriv;

    shapeDeriv(0,0) = -0.125*(1.0-eta)*(1.0-zeta);
    shapeDeriv(0,1) =  0.125*(1.0-eta)*(1.0-zeta);
    shapeDeriv(0,2) =  0.125*(1.0+eta)*(1.0-zeta);
    shapeDeriv(0,3) = -0.125*(1.0+eta)*(1.0-zeta);
    shapeDeriv(0,4) = -0.125*(1.0-eta)*(1.0+zeta);
    shapeDeriv(0,5) =  0.125*(1.0-eta)*(1.0+zeta);
    shapeDeriv(0,6) =  0.125*(1.0+eta)*(1.0+zeta);
    shapeDeriv(0,7) = -0.125*(1.0+eta)*(1.0+zeta);
    
    shapeDeriv(1,0) = -0.125*(1.0-xi)*(1.0-zeta);
    shapeDeriv(1,1) = -0.125*(1.0+xi)*(1.0-zeta);
    shapeDeriv(1,2) =  0.125*(1.0+xi)*(1.0-zeta);
    shapeDeriv(1,3) =  0.125*(1.0-xi)*(1.0-zeta);
    shapeDeriv(1,4) = -0.125*(1.0-xi)*(1.0+zeta);
    shapeDeriv(1,5) = -0.125*(1.0+xi)*(1.0+zeta);
    shapeDeriv(1,6) =  0.125*(1.0+xi)*(1.0+zeta);
    shapeDeriv(1,7) =  0.125*(1.0-xi)*(1.0+zeta);

    shapeDeriv(2,0) = -0.125*(1.0-xi)*(1.0-eta);
    shapeDeriv(2,1) = -0.125*(1.0+xi)*(1.0-eta);
    shapeDeriv(2,2) = -0.125*(1.0+xi)*(1.0+eta);
    shapeDeriv(2,3) = -0.125*(1.0-xi)*(1.0+eta);
    shapeDeriv(2,4) =  0.125*(1.0-xi)*(1.0-eta);
    shapeDeriv(2,5) =  0.125*(1.0+xi)*(1.0-eta);
    shapeDeriv(2,6) =  0.125*(1.0+xi)*(1.0+eta);
    shapeDeriv(2,7) =  0.125*(1.0-xi)*(1.0+eta);
             
    return shapeDeriv;
}

void Hex8::InitializeElements(Nodes &Nodes){

    // Initialize the vector containing each element stiffness matrix.
    elStiffMatx.resize(nElements); 

    // Move after allocation for `elStiffMatx`.
    elStiffMatxVariant = &elStiffMatx;

    // Initialize the storages for int-pt stresses/strains
    elStres.resize(nElements); elStran.resize(nElements); elDStran.resize(nElements);

    if (materialModel=="ElastoPlastic"){
        elStran_e.resize(nElements); elStran_e_old.resize(nElements);   // Elastic strain tensors
        elStran_p.resize(nElements); elStran_p_old.resize(nElements);   // Plastic strain tensors
        elStran_eq.resize(nElements); elStran_eq_old.resize(nElements); // Equivalent platic strain
        elStres_eq.resize(nElements); // Equivalent (von Mises) stress
    }

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd8x3 dummyElNodCoord; // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd3> dummyElemGauss(nElGauss); // For element Gauss points.

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

            elStran_e_old.at(iElem).resize(nElGauss);
            elStran_p_old.at(iElem).resize(nElGauss); 
            elStran_eq_old.at(iElem).resize(nElGauss);

            // Initilize to zeros.
            for (int iGaus=0; iGaus<nElGauss; iGaus++){
                
                elStran_e.at(iElem).at(iGaus).setZero();
                elStran_p.at(iElem).at(iGaus).setZero();
                elStres_eq.at(iElem).at(iGaus) = 0;
                elStran_eq.at(iElem).at(iGaus) = 0;

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
            dummyElNodCoord(iNod, 2) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(2);

        }
        elemNodCoord.at(iElem) = dummyElNodCoord;

        // Loop through integration points.
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            /*
            In 3D the values are directly set, so we have to initialize it to zeros or we will get garbage
            */ 
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

RowVecd3 Hex8::getGaussCart(RowVecd8& sFunc, Matd8x3& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Hex8::CalcCartDeriv(Matd8x3& elNodCoord, Matd3x8& sFuncDeriv, const double& wt, double& intVol, Matd3x8& cartDeriv, Matd6x24& strainMat){

    // Calculates the jacobian matrix J_jj = dN_ji x_ij
    Matd3x3 jacMat = sFuncDeriv*elNodCoord;

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

        strainMat(0,3*iNod)   = cartDeriv(0,iNod);
        strainMat(1,3*iNod+1) = cartDeriv(1,iNod);
        strainMat(2,3*iNod+2) = cartDeriv(2,iNod);
        
        strainMat(3,3*iNod+1) = cartDeriv(2,iNod);
        strainMat(3,3*iNod+2) = cartDeriv(1,iNod);
        
        strainMat(4,3*iNod)   = cartDeriv(2,iNod);
        strainMat(4,3*iNod+2) = cartDeriv(0,iNod);
        
        strainMat(5,3*iNod)   = cartDeriv(1,iNod);
        strainMat(5,3*iNod+1) = cartDeriv(0,iNod);
    }
}

void Hex8::CalcElemStiffMatx(T_DMatx DMatx){

    // Matd6x24 dummyBu;   // dummy for strain matrix.
    double dummydVol;   // dummy for int-pt volume.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

        // Integration over all Gauss points.
        for (int iGaus=0; iGaus<nElGauss; iGaus++){

            const Matd6x24& dummyBu = BuMat.at(iElem).at(iGaus); // Strain matrix for the given gauss point.
            dummydVol = intPtVol.at(iElem).at(iGaus);  // Volume of the current integration point 

            // [B_kl]^T D_kk B_kl
            elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*std::get<Matd6x6>(DMatx)*dummyBu*dummydVol;
        }
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elStiffMatx)
    //     cout << iStifMat << "\n\n";
    // cout << std::setprecision(15) << elStiffMatx.at(44)(16,21) << "\n";
    
    // for (int i = 0; i < elStiffMatx.at(45).rows(); ++i) {
    //     for (int j = 0; j < elStiffMatx.at(45).cols(); ++j) {
    //         std::cout << elStiffMatx.at(45)(i, j);
    //         if (j < elStiffMatx.at(45).cols() - 1) {
    //             std::cout << ", "; // Use a comma to separate columns
    //         }
    //     }
    //     std::cout << "\n"; // New line after each row
    // }
    // std::cout << "\n"; // New line after each row

}

void Hex8::CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount){

    ColVecd24 dummyDisp; // for element nodal displacement.
    ColVecd24 dummyForc; // for element nodal internal force.

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
            elStres.at(iElem).at(iGaus) = std::get<Matd6x6>(DMatx)*elStran.at(iElem).at(iGaus);

            dummyForc = BuMat.at(iElem).at(iGaus).transpose()*elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);

            // Sometimes it throws segmentation fault, have no idea why (⊙_⊙)？
            for(int pom=0; pom<nElDispDofs; pom++){
                Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
            }

            // Nodal values
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd6>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd6>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += 1;
            }
        }
    }

    // // TODO: For debug!
    // cout << elStran.at(0).at(0) << "\n\n";
}

void Hex8::CalcNodVals( T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<int>& nodCount){

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

                std::get<std::vector<ColVecd6>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd6>>(nodStran_e).at(*iNod2) += elStran_e.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd6>>(nodStran_p).at(*iNod2) += elStran_p.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd6>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus);
                nodStran_eq.at(*iNod2) += elStran_eq.at(iElem).at(iGaus);
                nodStres_eq.at(*iNod2) += elStres_eq.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += 1;
            }
        }
    }
}

void Hex8::CalcElDStran(const double* globalBuffer){

    ColVecd24 dummyDisp; // for element nodal displacement.

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

void Hex8::CalcElStran(const double* globalBuffer){

    ColVecd24 dummyDisp; // for element nodal displacement.

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

void Hex8::CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep){

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
                
                plasticMat->ReturnMapping3D(elDStran.at(iElem).at(iGaus),
                                            elStres.at(iElem).at(iGaus),
                                            elStran_e.at(iElem).at(iGaus),
                                            elStran_p.at(iElem).at(iGaus),
                                            elStran_eq.at(iElem).at(iGaus), 
                                            elStres_eq.at(iElem).at(iGaus),
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
                
                plasticMat->ReturnMapping3D(elDStran.at(iElem).at(iGaus),
                                            elStres.at(iElem).at(iGaus),
                                            elStran_e.at(iElem).at(iGaus),
                                            elStran_p.at(iElem).at(iGaus),
                                            elStran_eq.at(iElem).at(iGaus), 
                                            elStres_eq.at(iElem).at(iGaus),
                                            elStran_e_old.at(iElem).at(iGaus),
                                            elStran_p_old.at(iElem).at(iGaus),
                                            elStran_eq_old.at(iElem).at(iGaus), iStep);

                const Matd6x24& dummyBu = BuMat.at(iElem).at(iGaus); // Strain matrix for the given gauss point.
                dummydVol = intPtVol.at(iElem).at(iGaus);  // Volume of the current integration point 

                // [B_kl]^T D_kk B_kl
                elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*std::get<Matd6x6>(plasticMat->getDMatx())*dummyBu*dummydVol;
            }
        }
    }
}

void Hex8::CalcFint(double* Fint){

    ColVecd24 dummyForc; // for element nodal internal force.

    for(int iElem=0; iElem<nElements; iElem++){
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            dummyForc = BuMat.at(iElem).at(iGaus).transpose()*elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);

            for(int pom=0; pom<nElDispDofs; pom++){
                Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
            }
        }
    }
}

void Hex8::getNew(){

    elStran_e_old = elStran_e;
    elStran_p_old = elStran_p;
    elStran_eq_old = elStran_eq;

}