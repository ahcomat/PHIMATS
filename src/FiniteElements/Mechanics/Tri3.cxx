#include<iostream>
#include<algorithm>

#include"FiniteElements/Mechanics/Tri3.h"

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

Tri3::Tri3(H5IO &H5File_in, Nodes &Nodes, int iSet, string matModel, Logger& logger)
    : BaseElemMech(2, 3, 1, 2, 3, 6,  matModel, logger){ // nElDim, nElNodes, nElGauss, dispDofs, nElStres, nElDispDofs

    if (materialModel != "Elastic" && materialModel != "ElastoPlastic") {
        
        throw std::invalid_argument("Invalid material model: < " + materialModel + " >\nAllowed models are < Elastic, ElastoPlastic >\n");
    }

    InitShapeFunc();
    ReadElementsData(H5File_in, iSet);
    InitializeElements(Nodes);
}

Tri3::~Tri3(){

    // Exit message
    cout << "Tri3 elements exited correctly" << "\n";
}

void Tri3::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    vector<double> ip = {1.0/3.0, 1.0/3.0};
    gaussPts.push_back(ip);

    // for (auto& s : gaussPts[0])
    //     cout << s << "\n";

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nElNodes);
    shapeFuncDeriv.resize(nElNodes);

    for(int i=0; i<nElGauss; i++){
        shapeFunc.at(i) =  CalcShapeFunc(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
        shapeFuncDeriv.at(i) = CalcShapeFuncDeriv(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
    }
}

RowVecd3 Tri3:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd3 shape;

    shape(0) = xi;
    shape(1) = eta;
    shape(2) = 1 - xi - eta;

    return shape;
}

Matd2x3 Tri3::CalcShapeFuncDeriv(double xi, double eta){

    // dN_ji
    Matd2x3 shapeDeriv;

    shapeDeriv(0,0) = 1;
    shapeDeriv(0,1) = 0;
    shapeDeriv(0,2) = -1;

    shapeDeriv(1,0) = 0;
    shapeDeriv(1,1) = 1;
    shapeDeriv(1,2) = -1;
             
    return shapeDeriv;
}

void Tri3::InitializeElements(Nodes &Nodes){

    nDof = dispDofs*nNodes;      // Calc total number of dips DOFs for element set.

    // Initialize the storages for int-pt stresses/strains
    elStres.resize(nElements); elStran.resize(nElements);      

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd3x2 dummyElNodCoord; // For node coordinates.

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

RowVecd2 Tri3::getGaussCart(RowVecd3& sFunc, Matd3x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Tri3::CalcCartDeriv(Matd3x2& elNodCoord, Matd2x3& sFuncDeriv, const double& wt, double& intVol, Matd2x3& cartDeriv, Matd3x6& strainMat){

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

void Tri3::CalcElemStiffMatx(T_DMatx CMatx){

    Matd3x3 DMat = std::get<Matd3x3>(CMatx);

    elStiffMatx.resize(nElements); // Initialize the vector containing each element stiffness matrix.

    double dummydVol;   // dummy for int-pt volume.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.         

        // Integration over all Gauss points.
        for (int iGauss=0; iGauss<nElGauss; iGauss++){

            const Matd3x6& dummyBu = BuMat.at(iElem).at(iGauss); // Strain matrix for the given gauss point.
            dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current integration point 

            // [B_kl]^T D_kk B_kl
            elStiffMatx.at(iElem).noalias() += dummyBu.transpose()*DMat*dummyBu*dummydVol;
        }  
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elStiffMatx)
    //     cout << iStifMat << "\n\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
}

void Tri3::CalcStres(T_DMatx CMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<double>& nodCount){

    ColVecd6 dummyDisp; // for element nodal displacement.
    ColVecd6 dummyForc; // for element nodal internal force.

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

            for(int pom=0; pom<nElDispDofs; pom++){
                Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
            }

            // Nodal values
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd3>>(nodStran).at(*iNod2) += elStran.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                std::get<std::vector<ColVecd3>>(nodStres).at(*iNod2) += elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += intPtVol.at(iElem).at(iGaus);
            }
        }
    }

    // // TODO: For debug!
    // cout << elStran.at(0).at(0) << "\n\n";
}

void Tri3::CalcElStran(const double* globalBuffer){
    
}

void Tri3::CalcNodVals( T_nodStres& nodStres, T_nodStres& nodStran, T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<double>& nodStres_h, vector<double>& nodRho, vector<double>& nodCount){
}

void Tri3::CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep){

}

