/**
 * @file Quad4.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing quad4 elements. 
 * @date 2024-05-22
 * 
 * @copyright Copyright (c) 2024
 * 
 * Updates (when, what and who)
 * 
 */

#include<iostream>
#include<algorithm>

#include"FiniteElements/Quad4.h"

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
    : Elements(2, 4, 2, 3, 8, 4){ // nDim, nElNodes, dispDofs, nStres, nElDispDofs, nGauss

    InitShapeFunc();
    ReadElementsData(H5File_in);
    InitializeElements(Nodes);
    // InitPETSC();
}

Quad4::~Quad4(){

    // Exit message
    cout << "Quad4 elements exited correctly" << "\n";
}

void Quad4::InitShapeFunc(){

    // Initialize the Gauss points vectors.
    vector<double> ip = {-0.57735027, 0.57735027};
    vector<double> dummy(nDim);
    for(int iGauss=0; iGauss<ip.size(); iGauss++){
        for(int jGauss=0; jGauss<ip.size(); jGauss++){
            dummy.at(0) = ip.at(iGauss);
            dummy.at(1) = ip.at(jGauss);
            gaussPts.push_back(dummy);
        }
    }

    //Initialize shape functions and derivatives in natural coordinates.
    shapeFunc.resize(nGauss);
    shapeFuncDeriv.resize(nGauss);

    for(int i=0; i<nGauss; i++){
        shapeFunc.at(i) = getShapeFunc(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
        shapeFuncDeriv.at(i) = getShapeFuncDeriv(gaussPts.at(i).at(0), gaussPts.at(i).at(1));
    }
}

RowVecd4 Quad4::getShapeFunc(double xi, double eta){

    // N_i
    RowVecd4 shape;

    shape(0) = (1.0-eta-xi+xi*eta)*0.25;
    shape(1) = (1.0-eta+xi-xi*eta)*0.25;
    shape(2) = (1.0+eta+xi+xi*eta)*0.25;
    shape(3) = (1.0+eta-xi-xi*eta)*0.25;

    return shape;
}

Matd2x4 Quad4::getShapeFuncDeriv(double xi, double eta){

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
        elemDispDof.at(iElem) = getElemDispDof(iElem);
    }

    // for (auto& s : elemNodeConn[0])
    //     cout << s << "\n"; 

    // /**
    //  * Read Dirichlet BCs
    //  */
    // dsetName = "SimulationParameters/nPresDofs";
    // nPresDofs = H5File_in.ReadScalar(dsetName);
    // PetscMalloc1(nPresDofs, &presDofs);
    // PetscMalloc1(nPresDofs, &presVals); 
    // PetscMalloc1(nTotDof, &Fint); 

    // vector<double> dummy2(3);
    // for (int iPresDof=0; iPresDof<nPresDofs; iPresDof++){
    //     // Read values
    //     dsetName = "PrescribedDOFs/Prescribed_"+to_string(iPresDof);
    //     H5File_in.ReadFieldDoub1D(dsetName, dummy2);
    //     // Assign values
    //     presDofs[iPresDof] = nDim*dummy2.at(0)+dummy2.at(1); // nDim*iNode+dof
    //     presVals[iPresDof] = dummy2.at(2);
    // }
}

vector<int> Quad4::getElemDispDof(int iElem){

    vector<int> dispDof(nElDispDofs);
    for(int iNod=0; iNod<nElNodes; iNod++){

        dispDof.at(2*iNod) = 2*elemNodeConn.at(iElem).at(iNod);
        dispDof.at(2*iNod+1) = 2*elemNodeConn.at(iElem).at(iNod)+1;
    }

    return dispDof;
}

// vector<int> Quad4::getNodeConnect(int iElem){

//     return elemNodeConn.at(iElem);
// }

void Quad4::InitializeElements(Nodes &Nodes){

    nNodes = Nodes.getNNodes();   // Get total number of nodes. 
    nTotDof = dispDofs*nNodes;    // Calc total number of DOFs.

    // Initialize the storages for int-pt stresses/strains
    elStres.resize(nElements); elStran.resize(nElements);

    // Initialize the storages for nodal stresses/strains and counts
    nodStres.resize(nNodes); nodStran.resize(nNodes); nodCount.resize(nNodes);      

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    // Eigen::Matrix<double, 4, 2> dummyElNodCoord; // For node coordinates.
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
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss), BuMat.at(iElem).at(iGauss));
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

void Quad4::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, double& intVol, Matd2x4& cartDeriv, Matd3x8& strainMat){

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
    double dummydVol;   // dummy for int point volume.

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

    // cout << elStiffMatx.at(10) << "\n";
    elStiffMatxVariant = &elStiffMatx;
}

// void Quad4::setDirichBC(){

//     // MatZeroRowsColumns(A, nPresDofs, presDofs, 1.0, NULL, NULL);
//     MatZeroRows(A, nPresDofs, presDofs, 1.0, NULL, NULL);

//     MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

//     VecSetValues(b, nPresDofs, presDofs, presVals, ADD_VALUES); // Check for incremental loads "ADD_VALUES"
//     VecAssemblyBegin(b); VecAssemblyEnd(b);
// }

// Vec& Quad4::getB(){

//     return b;
// }

// Vec& Quad4::getX(){

//     return x;
// }

// Mat& Quad4::getA(){

//     return A;
// }

// void Quad4::CalcStres(Matd3x3 DMatx, bool nodStrFlag){

//     VecGetArrayRead(x, &globalBuffer);

//     ColVecd8 dummyDisp; // for element nodal displacement.
//     ColVecd8 dummyForc; // for element nodal internal force.

//     // Set zeros, otherwise will get garbage \o/ memory values.
//     // Internal force
//     for(int iDof=0; iDof<nTotDof; iDof++){
        
//         Fint[iDof] = 0;
//     }
//     // Nodal stress
//     if(nodStrFlag){
//         for(int iNod=0; iNod<nNodes; iNod++){
//             nodStran.at(iNod).setZero();
//             nodStres.at(iNod).setZero();
//             nodCount.at(iNod) = 0;

//             Fint[iNod] = 0;
//         }
//     }


//     /**
//      * Integration point values. 
//      * 
//      */
//     for(int iElem=0; iElem<nElements; iElem++){

//         // Get element nodal displacements from the solution vector. 
//         for(int iDof=0; iDof<nElDispDofs; iDof++){

//             dummyDisp(iDof) = globalBuffer[elemDispDof.at(iElem).at(iDof)];
//         }

//         // Gauss points
//         for(int iGaus=0; iGaus<nGauss; iGaus++){

//             // Int pt values
//             elStran.at(iElem).at(iGaus) = BuMat.at(iElem).at(iGaus)*dummyDisp;
//             elStres.at(iElem).at(iGaus) = DMatx*elStran.at(iElem).at(iGaus);

//             dummyForc = BuMat.at(iElem).at(iGaus).transpose()*elStres.at(iElem).at(iGaus)*intPtVol.at(iElem).at(iGaus);

//             // Sometimes it throws segmentation fault, have no idea why (⊙_⊙)？
//             for(int pom=0; pom<nElDispDofs; pom++){

//                 Fint[elemDispDof.at(iElem).at(pom)] += dummyForc(pom);
//             }

//             // Nodal values
//             if(nodStrFlag){
//                 for(auto iNod=elemNodeConn.at(iElem).begin(); iNod!=elemNodeConn.at(iElem).end(); iNod++){

//                     nodStran.at(*iNod) += elStran.at(iElem).at(iGaus);
//                     nodStres.at(*iNod) += elStres.at(iElem).at(iGaus);
//                     nodCount.at(*iNod) += 1;
//                 }
//             }
//         }
//     }

//     VecRestoreArrayRead(x, &globalBuffer);

//     /**
//      * Nodal value by averaging
//      * 
//      */

//     if(nodStrFlag){
//         for(int iNod=0; iNod<nNodes; iNod++){

//             nodStran.at(iNod) = nodStran.at(iNod)/nodCount.at(iNod);
//             nodStres.at(iNod) = nodStres.at(iNod)/nodCount.at(iNod);
//         }
//     }
// }

// void Quad4::WriteOut(H5IO &H5File_out){

//     // Displacements
//     VecGetArrayRead(x, &globalBuffer);
//     H5File_out.WriteArray_1D("Disp", nTotDof, globalBuffer);
//     VecRestoreArrayRead(x, &globalBuffer);

//     H5File_out.WriteArray_1D("Force", nTotDof, Fint);

//     // Stresses and strains
//     H5File_out.WriteStres3("Strain", nNodes, nStres, nodStran);
//     H5File_out.WriteStres3("Stress", nNodes, nStres, nodStres);
// }
