#include<iostream>
#include<algorithm>
#include <omp.h>

#include"FiniteElements/Transport/Quad4T.h"

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

Quad4T::Quad4T(H5IO &H5File_in, Nodes &Nodes, int iSet)
    : BaseElemTrans(2, 4, 4, 4){ // nElDim, nElNodes, nElConDofs, nElGauss 

    InitShapeFunc();
    ReadElementsData(H5File_in, iSet);
    InitializeElements(Nodes);
}   

Quad4T::~Quad4T(){

    // Exit message
    cout << "Quad4T elements exited correctly" << "\n";
}

void Quad4T::InitShapeFunc(){

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

RowVecd4 Quad4T:: CalcShapeFunc(double xi, double eta){

    // N_i
    RowVecd4 shape;

    shape(0) = (1.0-eta-xi+xi*eta)*0.25;
    shape(1) = (1.0-eta+xi-xi*eta)*0.25;
    shape(2) = (1.0+eta+xi+xi*eta)*0.25;
    shape(3) = (1.0+eta-xi-xi*eta)*0.25;

    return shape;
}

Matd2x4 Quad4T::CalcShapeFuncDeriv(double xi, double eta){

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

void Quad4T::InitializeElements(Nodes &Nodes){

    // Initialize the storages for int-pt flux 
    elFlux.resize(nElements);      

    elemNodCoord.resize(nElements); // Initialize the size of node coordinates.
    Matd4x2 dummyElNodCoord; // For node coordinates.

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

RowVecd2 Quad4T::getGaussCart(RowVecd4& sFunc, Matd4x2& elNodCoord){

    return sFunc*elNodCoord;  // N_i x_ij
}

void Quad4T::CalcCartDeriv(Matd4x2& elNodCoord, Matd2x4& sFuncDeriv, const double& wt, double& intVol, Matd2x4& cartDeriv){

    // Calculates the jacobian matrix J_jj = dN_ji x_ij
    Matd2x2 jacMat = sFuncDeriv*elNodCoord;

    // Jacobian determinant.
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
}

void Quad4T::CalcElemStiffMatx(T_DMatx KMatx, double s){

    Matd2x2 KMat = std::get<Matd2x2>(KMatx);

    elStiffMatx.resize(nElements);
    elCapMatx.resize(nElements);
    vector<Matd4x4> elKdMatx(nElements);

    double dummydVol;   // dummy for int-pt volume.

    // // Set the number of threads
    // omp_set_num_threads(4); // Set to the desired number of threads
    // // Parallelize the outer loop
    // #pragma omp parallel for

    // Loop through all elements.
    for(int iElem=0; iElem<nElements; iElem++){

        elStiffMatx.at(iElem).setZero();
        elKdMatx.at(iElem).setZero();  // Must be populated with zeros.
        elCapMatx.at(iElem).setZero(); // Must be populated with zeros.                  

        // Integration over all Gauss points.
        for (int iGauss=0; iGauss<nElGauss; iGauss++){

            const Matd2x4& dummyBMat = BMat.at(iElem).at(iGauss); // derivative matrix for the given gauss point.
            const RowVecd4& dummyShFunc = shapeFunc.at(iGauss);
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
    // cout << elStiffMatx.at(0) << "\n";

    // Pointer to the vector, not the vector itself.
    elStiffMatxVariant = &elStiffMatx;
    elCapMatxVariant = &elCapMatx;
}

void Quad4T::CalcFlux(T_DMatx KMatx, const double* globalBuffer, T_nodStres& nodFlux, vector<double>& nodCount){

    ColVecd4 dummyCon; // for element nodal concentration.
    int iNode;  // counter for the number of nodes.

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
            iNode = 0;
            for(auto iNod2=elemNodeConn.at(iElem).begin(); iNod2!=elemNodeConn.at(iElem).end(); iNod2++){

                std::get<std::vector<ColVecd2>>(nodFlux).at(*iNod2) += elFlux.at(iElem).at(iGaus)*shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                nodCount.at(*iNod2) += shapeFunc.at(iGaus)[iNode]*wts.at(iGaus);
                iNode += 1;
            }
        }
    }

    // // TODO: For debug!
    // cout << elFlux.at(0).at(0) << "\n\n";
}
