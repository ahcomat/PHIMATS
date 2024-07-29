#include<iostream>
#include<algorithm>

#include"FiniteElements/Transport/Tri3T.h"

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

Tri3T::Tri3T(H5IO &H5File_in, Nodes &Nodes, int iSet)
    : BaseElemTransport(2, 3, 3, 3){ // nElDim, nElNodes, nElConDofs, nElGauss 

    InitShapeFunc();
    ReadElementsData(H5File_in, iSet);
    InitializeElements(Nodes);
}   

Tri3T::~Tri3T(){

    // Exit message
    cout << "Tri3T elements exited correctly" << "\n";
}

void Tri3T::InitShapeFunc(){

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

RowVecd3 Tri3T:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd3 shape;

    shape(0) = xi;
    shape(1) = eta;
    shape(2) = 1 - xi - eta;

    return shape;
}

Matd2x3 Tri3T::CalcShapeFuncDeriv(double xi, double eta){

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

void Tri3T::InitializeElements(Nodes &Nodes){

    // Initialize the storages for int-pt flux 
    elFlux.resize(nElements);      

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd3x2 dummyElNodCoord; // For node coordinates.

    gaussPtCart.resize(nElements);  // Initialize the size of the Cart Gauss points.
    vector<RowVecd2> dummyElemGauss(nElGauss); // For element Gauss points.

    BMat.resize(nElements);   // Initialize the size of BMatrix.

    intPtVol.resize(nElements);   
    vector<double> dummyIntVol(nElGauss);  // For integration point volume.

    // Loop through elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elFlux.at(iElem).resize(nElGauss);
        gaussPtCart.at(iElem).resize(nElGauss);
        BMat.at(iElem).resize(nElGauss);

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
            CalcCartDeriv(dummyElNodCoord, shapeFuncDeriv.at(iGauss), wts.at(iGauss), dummyIntVol.at(iGauss), BMat.at(iElem).at(iGauss));
        }
        gaussPtCart.at(iElem) = dummyElemGauss;
        intPtVol.at(iElem) = dummyIntVol;
    }
}

RowVecd2 Tri3T::getGaussCart(RowVecd3& sFunc, Matd3x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Tri3T::CalcCartDeriv(Matd3x2& elNodCoord, Matd2x3& sFuncDeriv, const double& wt, double& intVol, Matd2x3& cartDeriv){

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

void Tri3T::CalcElemStiffMatx(T_DMatx KMatx, double s){

    Matd2x2 KMat = std::get<Matd2x2>(KMatx);

    elStiffMatx.resize(nElements);
    elCapMatx.resize(nElements);
    vector<Matd3x3> elKdMatx(nElements);

    double dummydVol;   // dummy for int-pt volume.

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero();
        elKdMatx.at(iElem).setZero();  // Must be populated with zeros.
        elCapMatx.at(iElem).setZero(); // Must be populated with zeros.                  

        // Integration over all Gauss points.
        for (int iGauss=0; iGauss<nElGauss; iGauss++){

            const Matd2x3& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
            const RowVecd3& dummyShFunc = shapeFunc.at(iGauss);
            dummydVol = intPtVol.at(iElem).at(iGauss);  // Volume of the current integration point 

            // [B_ji]^T k_jj B_ji
            elKdMatx.at(iElem).noalias() += dummyBMat.transpose()*KMat*dummyBMat*dummydVol;
            // [N_i]^T s N_i
            elCapMatx.at(iElem).noalias() += s*(dummyShFunc.transpose()*dummyShFunc)*dummydVol;
        }

        elStiffMatx.at(iElem) = dt*elKdMatx.at(iElem) + elCapMatx.at(iElem);
    }

    // // TODO: For debug!
    // for (auto& iStifMat : elKdMatx)
    //     cout << iStifMat << "\n\n";
    // cout << elStiffMatx.at(0) << "\n\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
    elCapMatxVariant = &elCapMatx;
}

void Tri3T::CalcFlux(T_DMatx KMatx, const double* globalBuffer, T_nodStres& nodFlux, vector<double>& nodCount){

    ColVecd3 dummyCon; // for element nodal concentration.

    // Integration point values.
    for(int iElem=0; iElem<nElements; iElem++){

        // Get element nodal concentration from the solution vector. 
        for(int iDof=0; iDof<nElConDofs; iDof++){
            dummyCon(iDof) = globalBuffer[elemConDof.at(iElem).at(iDof)];
        }

        // Gauss points
        for(int iGaus=0; iGaus<nElGauss; iGaus++){

            // Int pt flux
            elFlux.at(iElem).at(iGaus) = -(std::get<Matd2x2>(KMatx)*BMat.at(iElem).at(iGaus))*dummyCon;

            // Nodal values
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                // TODO: There are some issues here
                std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += elFlux.at(iElem).at(iGaus);
                nodCount.at(*iNod2) += 1;
            }
        }
    }

    // // TODO: For debug!
    // cout << elFlux.at(0).at(0) << "\n\n";
}

double Tri3T::CalcAvCon(const double* globalBuffer){

    double IntPtCon, AvCon = 0, TotVol = 0;
    ColVecd3 dummyCon; 

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