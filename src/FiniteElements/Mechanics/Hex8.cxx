#include<iostream>

#include"FiniteElements/Mechanics/Hex8.h"

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

Hex8::Hex8(H5IO &H5File_in, Nodes &Nodes, int iSet)
    : BaseElemMech(3, 8, 3, 6, 24, 8){ // nElDim, nElNodes, dispDofs, nElStres, nElDispDofs, nElGauss

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

    for(int iGauss=0; iGauss<2; iGauss++){
        for(int jGauss=0; jGauss<2; jGauss++){
            for(int kGauss=0; kGauss<2; kGauss++){

                dummy.at(0) = ip.at(iGauss);
                dummy.at(1) = ip.at(jGauss);
                dummy.at(2) = ip.at(kGauss);

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

    // Initialize the storages for int-pt stresses/strains
    elStres.resize(nElements); elStran.resize(nElements);      

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
        for(int iGauss=0; iGauss<nElGauss; iGauss++){
        
            // Cart coord of iGauss point.
            dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss), BuMat.at(iElem).at(iGauss));
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

    Matd6x6 DMat = std::get<Matd6x6>(DMatx);

    elStiffMatx.resize(nElements); // Initialize the vector containing each element stiffness matrix.

    // Matd6x24 dummyBu;   // dummy for strain matrix.
    double dummydVol;   // dummy for int-pt volume.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

        // Integration over all Gauss points.
        for (int iGauss=0; iGauss<nElGauss; iGauss++){

            const Matd6x24& dummyBu = BuMat.at(iElem).at(iGauss); // Strain matrix for the given gauss point.
            dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current integration point 

            // [B_kl]^T D_kk B_kl
            elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*DMat*dummyBu*dummydVol;
        }  
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elStiffMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
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
