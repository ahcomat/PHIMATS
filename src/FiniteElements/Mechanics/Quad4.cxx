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

Quad4::Quad4(H5IO &H5File_in, Nodes &Nodes, int iSet)
    : BaseElemMech(2, 4, 2, 3, 8, 4){ // nElDim, nElNodes, dispDofs, nElStres, nElDispDofs, nElGauss

    InitShapeFunc();
    ReadElementsData(H5File_in, iSet);
    InitializeElements(Nodes);
}   

Quad4::~Quad4(){

    // Exit message
    cout << "Quad4 elements exited correctly" << "\n";
}

void Quad4::InitShapeFunc(){

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
    elStres.resize(nElements); elStran.resize(nElements); 

    if (materialModel=="ElastoPlastic"){
        elStran_e.resize(nElements); elStran_p.resize(nElements); // Elastic and plastic strain tensors
        elStres_eq.resize(nElements); elStran_eq.resize(nElements); // Equivalent platic strain/stress
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

        if (materialModel=="ElastoPlastic"){

            elStran_e.at(iElem).resize(nElGauss);
            elStran_p.at(iElem).resize(nElGauss); 
            elStres_eq.at(iElem).resize(nElGauss);
            elStran_eq.at(iElem).resize(nElGauss);

            // Initilize to zeros.
            for (int iGaus=0; iGaus<nElGauss; iGaus++){

                elStran_e.at(iElem).at(iGaus).setZero();
                elStran_p.at(iElem).at(iGaus).setZero();
                elStran_eq.at(iElem).at(iGaus) = 0;
                elStres_eq.at(iElem).at(iGaus) = 0;
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
        for(int iGauss=0; iGauss<nElGauss; iGauss++){
        
            // Cart coord of iGauss point.
            dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss), BuMat.at(iElem).at(iGauss));
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

void Quad4::CalcElemStiffMatx(T_DMatx DMatx){

    Matd3x3 DMat = std::get<Matd3x3>(DMatx);

    double dummydVol;   // dummy for int-pt volume.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

        // Integration over all Gauss points.
        for (int iGauss=0; iGauss<nElGauss; iGauss++){

            const Matd3x8& dummyBu = BuMat.at(iElem).at(iGauss); // Strain matrix for the given gauss point.
            dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current integration point 

            // [B_kl]^T D_kk B_kl
            elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*DMat*dummyBu*dummydVol;
        }  
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elStiffMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n";
}

void Quad4::CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount){

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
            elStres.at(iElem).at(iGaus) = std::get<Matd3x3>(DMatx)*elStran.at(iElem).at(iGaus);

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

void Quad4::CalcNodVals( T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<int>& nodCount){


}

void Quad4::CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep){

}

void Quad4::CalcFint(double* Fint){

}

