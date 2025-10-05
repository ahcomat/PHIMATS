#include "FiniteElements/Trapping/BaseElemTrap.h"
#include "Materials/Trapping/TrapGB.h"
#include "Materials/Trapping/TrapPhase.h"
#include "Materials/Trapping/MechTrap.h"

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

void BaseElemTrap::ReadElementsData(H5IO &H5File_in, int iSet){

    string dsetName;
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nElements";
    nElements = H5File_in.ReadScalar(dsetName);
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/nNodes";
    nNodes = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/dt";
    dt = H5File_in.ReadScalar(dsetName);
    dsetName = "SimulationParameters/Trapping";
    Trapping = H5File_in.ReadString(dsetName);
    
    // Initialize the size.
    elemNodeConn.resize(nElements);  
    elemConDof.resize(nElements);
    elemIDs.resize(nElements);

    // Read global element IDs
    dsetName = "Elements/ElementSet_"+std::to_string(iSet)+"/ElementSetIDs";
    H5File_in.ReadField1D(dsetName, elemIDs);

    // Read node connectivity
    dsetName = "SimulationParameters/nTotElements";
    int totElements = H5File_in.ReadScalar(dsetName);  // Total number of elements
    vector<vector<int>> totElNodConnectivity(totElements);  // Node connectivity for all elements
    dsetName = "NodeConnectivity";
    H5File_in.ReadField2D(dsetName, totElements, nElNodes, totElNodConnectivity);

    // Node connectivity for the element set
    for (int iElem=0; iElem<nElements; iElem++){

        elemNodeConn.at(iElem) = totElNodConnectivity.at(elemIDs.at(iElem));
        elemConDof.at(iElem) = totElNodConnectivity.at(elemIDs.at(iElem));
    }

    // TODO: For debug!
    // for (auto& s : elemNodeConn[0])
    //     cout << s << "\n"; 
}

void BaseElemTrap::ReadNodalStress(H5IO &H5File_stress, int iStep){

    H5File_stress.ReadField1D("Stress_h/Step_"+std::to_string(iStep), nod_sigma_h);
    H5File_stress.ReadField1D("Rho/Step_"+std::to_string(iStep), nod_rho);

}

void BaseElemTrap::CalcEquilibriumBC(BaseTrapping* mat, double* presVals, int* presDofs, const int nPresDofs, const double conB, const double T){

    try{

        if (Trapping=="MechTrapping" || Trapping=="MechTrappingPFF"){                   // Stresses and dislocations

            double Vh = dynamic_cast<MechTrap*>(mat)->get_Vh();
            double zeta_rho = dynamic_cast<MechTrap*>(mat)->get_zeta_rho();

            for (int iPresDof=0; iPresDof<nPresDofs; iPresDof++){
                presVals[iPresDof] = conB*exp((nod_sigma_h[presDofs[iPresDof]]*Vh +
                                               nod_rho[presDofs[iPresDof]]*zeta_rho)/(R*T));
            }

        } else {
            throw std::runtime_error("unsupported trapping type < " + Trapping + " >");
        }

    } catch (const std::runtime_error& e) {

            logger.log("\nException caught in BaseElemTrap::setEquilibriumBC:\n", "", false);
            logger.log("    " + std::string(e.what()), "", false);
            logger.log("\nCritical error encountered. Terminating!\n", "", false);
            exit(EXIT_FAILURE);
    
        }
}

