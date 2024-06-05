#include<iostream>
#include<algorithm>

#include"FiniteElements/Mechanics/Quad4.h"

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

Quad4::Quad4(H5IO &H5File_in, Nodes &Nodes)
    : BaseElemMech(2, 4, 2, 3, 8, 4){ // nDim, nElNodes, dispDofs, nStres, nElDispDofs, nGauss

    InitShapeFunc();
    ReadElementsData(H5File_in);
    InitializeElements(Nodes);
}   

Quad4::~Quad4(){


    // Exit message
    cout << "Quad4 elements exited correctly" << "\n";
}

void Quad4::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    vector<double> ip = {-0.57735027, 0.57735027};
    vector<double> dummy(nDim);
    for(std::size_t iGauss=0; iGauss<ip.size(); iGauss++){
        for(std::size_t jGauss=0; jGauss<ip.size(); jGauss++){
            dummy.at(0) = ip.at(iGauss);
            dummy.at(1) = ip.at(jGauss);
            gaussPts.push_back(dummy);
        }
    }

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nGauss);
    shapeFuncDeriv.resize(nGauss);

    for(int i=0; i<nGauss; i++){
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

void Quad4::ReadElementsData(H5IO &H5File_in){

    string dsetName;
    dsetName = "SimulationParameters/nElements";
    nElements = H5File_in.ReadScalar(dsetName);
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nNodes";
    nNodes = H5File_in.ReadScalar(dsetName);

    // Initialize the size.
    elemNodeConn.resize(nElements);  
    elemDispDof.resize(nElements);

    dsetName = "SimulationParameters/nElementSets";
    nElementSets = H5File_in.ReadScalar(dsetName);

    // Read node connectivity.
    vector<int> dummy(nElNodes);
    for (int iElem=0; iElem<nElements; iElem++){
        dsetName = "NodeConnectivity/Element_"+to_string(iElem);
        H5File_in.ReadFieldInt1D(dsetName, dummy);

        elemNodeConn.at(iElem) = dummy;
        elemDispDof.at(iElem) = CalcElemDispDof(iElem);
    }

    // for (auto& s : elemNodeConn[0])
    //     cout << s << "\n"; 
}

vector<int> Quad4::CalcElemDispDof(int iElem){

    vector<int> dispDof(nElDispDofs);
    for(int iNod=0; iNod<nElNodes; iNod++){

        dispDof.at(nDim*iNod) = nDim*elemNodeConn.at(iElem).at(iNod);
        dispDof.at(nDim*iNod+1) = nDim*elemNodeConn.at(iElem).at(iNod)+1;
    }

    return dispDof;
}

// vector<int> Quad4::getNodeConnect(int iElem){

//     return elemNodeConn.at(iElem);
// }

void Quad4::InitializeElements(Nodes &Nodes){

    nDof = dispDofs*nNodes;      // Calc total number of dips DOFs for element set.

    // Initialize the storages for int-pt stresses/strains
    elStres.resize(nElements); elStran.resize(nElements);

    // Initialize the storages for nodal stresses/strains and counts
    nodStres.resize(nNodes); nodStran.resize(nNodes); nodCount.resize(nNodes);      

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd4x2 dummyElNodCoord; // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.
    BuMat.resize(nElements);  // Initialize the size of BuMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nGauss);  // For integration point volume.

    // Loop through elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStres.at(iElem).resize(nGauss);
        elStran.at(iElem).resize(nGauss);
        BMat.at(iElem).resize(nGauss);
        BuMat.at(iElem).resize(nGauss);

        // Loop through nodes to get coordinates.
        for(int iNod=0; iNod<nElNodes; iNod++){

            dummyElNodCoord(iNod, 0) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(0);
            dummyElNodCoord(iNod, 1) = Nodes.getNodCoord(elemNodeConn.at(iElem).at(iNod)).at(1);
        }
        elemNodCoord.at(iElem) = dummyElNodCoord;

        // Loop through integration points.
        for(int iGauss=0; iGauss<nGauss; iGauss++){
        
            // Cart coord of iGauss point.
            dummyElemGauss.at(iGauss) = getGaussCart(shapeFunc.at(iGauss), dummyElNodCoord);
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss), BuMat.at(iElem).at(iGauss));
        }
        gaussPtCart.at(iElem) = dummyElemGauss;
        intPtVol.at(iElem) = dummyIntVol;
    }
}

// Matd4x2 Quad4::getElemNodCoord(int iElem){
    
//     return elemNodCoord.at(iElem);
// }

RowVecd2 Quad4::getGaussCart(RowVecd4& sFunc, Matd4x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Quad4::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv, Matd3x8& strainMat){

    // Calculates the jacobian matrix J_jj = dN_ji x_ij
    Matd2x2 jacMat = sFuncDeriv*elNodCoord;

    // Jacobian determinant. Since all the weights for quad4 are "1", it is also the volume of the integration point. 
    intVol = jacMat.determinant();
        
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

    elStiffMatx.resize(nElements); // Initialize the vector containing each element stiffness matrix.

    Matd3x8 dummyBu;    // dummy for strain matrix.
    double dummydVol;   // dummy for int-pt volume.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero(); // Must be populated with zeros.

        // Integration over all Gauss points.
        for (int iGauss=0; iGauss<nGauss; iGauss++){

            dummyBu = BuMat.at(iElem).at(iGauss); // Strain matrix for the given gauss point.
            dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current integration point 

            // [B_kl]^T D_kk B_kl
            elStiffMatx.at(iElem) += dummyBu.transpose()*std::get<Matd3x3>(DMatx)*dummyBu*dummydVol;
        }        
    }

    // cout << elStiffMatx.at(0) << "\n\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
}

void Quad4::CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint){

    ColVecd8 dummyDisp; // for element nodal displacement.
    ColVecd8 dummyForc; // for element nodal internal force.

    for(int iNod=0; iNod<nNodes; iNod++){
        nodStran.at(iNod).setZero();
        nodStres.at(iNod).setZero();
        nodCount.at(iNod) = 0;
    }

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal displacements from the solution vector. 
        for(int iDof=0; iDof<nElDispDofs; iDof++){
            dummyDisp(iDof) = globalBuffer[elemDispDof.at(iElem).at(iDof)];
        }

        // Gauss points
        for(int iGaus=0; iGaus<nGauss; iGaus++){

            // Int pt values
            elStran.at(iElem).at(iGaus) = BuMat.at(iElem).at(iGaus)*dummyDisp;
            elStres.at(iElem).at(iGaus) = std::get<Matd3x3>(DMatx)*elStran.at(iElem).at(iGaus);

            dummyForc = BuMat.at(iElem).at(iGaus).transpose()*elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);

            // Sometimes it throws segmentation fault, have no idea why (⊙_⊙)？
            for(int pom=0; pom<nElDispDofs; pom++){
                Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
            }

            // Nodal values
            for(auto iNod=elemNodeConn.at(iElem).begin(); iNod!=elemNodeConn.at(iElem).end(); iNod++){

                nodStran.at(*iNod) += elStran.at(iElem).at(iGaus);
                nodStres.at(*iNod) += elStres.at(iElem).at(iGaus);
                nodCount.at(*iNod) += 1;
            }
            
        }
    }

    // Number averaging the nodal values
    for(int iNod=0; iNod<nNodes; iNod++){
        
        nodStran.at(iNod) = nodStran.at(iNod)/nodCount.at(iNod);
        nodStres.at(iNod) = nodStres.at(iNod)/nodCount.at(iNod);
    }
}

void Quad4::WriteOut(H5IO &H5File_out){

    // Stresses and strains
    H5File_out.WriteStres3("Strain", nNodes, nStres, nodStran);
    H5File_out.WriteStres3("Stress", nNodes, nStres, nodStres);
}
