/**
 * @file BaseElemTrap.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for trapping elements. 
 * @date 2024-06-21
 * 
 * @copyright Copyright (C) 2025 Abdelrahman Hussein
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  
 */

#ifndef BASEELEMTRAP_H
#define BASEELEMTRAP_H

#include "H5IO.h"
#include "FiniteElements/BaseElements.h"
#include "Materials/Trapping/BaseTrapping.h"

class BaseElemTrap: public BaseElements{

public:

BaseElemTrap(int nElDim, int nElNodes, int nElGauss, int nElConDofs, Logger& logger)
    : BaseElements(nElDim, nElNodes, nElGauss, logger), nElConDofs(nElConDofs) {};

/**
 * @brief Get dimensions of the element. 
 * 
 * @return int 
 */
int get_nDim() const { return nElDim; };

/**
 * @brief Get the number of element Gauss points. 
 * 
 * @return int 
 */
int get_nGauss() const { return nElGauss; };

/**
 * @brief Get the number of element con DOFs. Same as number of nodes. 
 * 
 * @return int 
 */
int get_nElConDofs() const { return nElConDofs; };

/**
 * @brief Set new time increment `dt`.
 * 
 * @param dtNew 
 */
void set_dt(double dtNew) { dt = dtNew; };

/**
 * @brief Return a const reference to the `elemConDof`. 
 * 
 * @return const vector<vector<int>>& 
 */
const vector<vector<int>>& get_elemConDof() const { return elemConDof; };

/**
 * @brief Reads the data `nElements`, `nElementSets` and `elemNodeConn` from hdf5 file.
 * 
 * @param H5File_in  
 */
void ReadElementsData(H5IO &H5File_in, int iSet);

/**
 * @brief Function to read the nodal stresses and normalized dislocation density.
 * 
 * @param H5File_stress Output HDF5 file from mechanical simulation.
 * @param iStep Time step to read. 
 */
void ReadNodalStress(H5IO &H5File_stress, int iStep);

/**
 * @brief Calculates the equilibrium concentration at boundary nodes. 
 * 
 * @param mat Trapping material pointer.
 * @param presVals double* prescribed nodal values.
 * @param presDofs int* container for the degrees of freedom of the boundary nodes. 
 * @param nPresDofs int total number of degrees of freedom.
 * @param conB double constant boundary concentration. 
 * @param T double temperature. 
 */
void CalcEquilibriumBC(BaseTrapping* mat, double* presVals, int* presDofs, const int nPresDofs, const double conB, const double T);


/**
 * @brief Return const reference to the vector of element capacitance matrix c_ii.
 * 
 * @return const vector<T_ElStiffMatx>& 
 */
const T_ElStiffMatx& getElCapMatx() const { return elCapMatxVariant; }

/**
 * @brief Get a constant reference to `elCon`.
 * 
 * @return const std::vector<std::vector<double>>& 
 */
const std::vector<std::vector<double>>& getElCon() const ;

/**
 * @brief Get the int-pt coordinates
 * 
 * @param glIntPtCoords 
 */
virtual void getInPtCoords(T_nodStres& glIntPtCoords) = 0;

/**
 * @brief Calculates the element stiffness and capacitances matrix.
 * 
 * @param mat 
 * @param T 
 */
virtual void CalcElemStiffMatx(BaseTrapping* mat, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr = nullptr) = 0;

/**
 * @brief Calculates the volume averaged concentration
 * 
 * @param globalBuffer Solution vector
 * @return double Volume averaged concentration
 */
virtual double CalcAvCon(const double* globalBuffer) = 0;

// /**
//  * @brief Evaluates the gradients of scalar field at the int-points and maps them to the nodes. 
//  * 
//  * @param nodGrad
//  * @param nodCount 
//  */
// virtual void CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapSigmaH) = 0;

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
virtual void CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr = nullptr) = 0;

/**
 * @brief 
 * 
 * @param conB 
 * @param mat 
 * @param FsrcBuffer 
 * @param T 
 * @param elPhi_d_ptr 
 */
virtual void CalcFsrc(const double conB, BaseTrapping* mat, double* FsrcBuffer, const double T, const std::vector<std::vector<double>>* elPhi_d_ptr) = 0;

/**
 * @brief Calculated the integration point concentration. 
 * 
 * @param globalBuffer 
 */
virtual void CalcElCon(const double* globalBuffer) = 0;

protected:      

/// @brief Number of element concentration (temperature) dofs.
const int nElConDofs;            

/// @brief Time increment.   
double dt;          

/// @brief Trapping type.
string Trapping;      

/// @brief Int-pt phi [nElGauss]. 
vector<vector<double>> el_gPhi;   

/// @brief Int-pt phi [nElGauss]. 
vector<vector<double>> el_gPhi_HAGB; 

/// @brief Int-pt phi [nElGauss]. 
vector<vector<double>> el_gPhi_LAGB; 

/// @brief Int-pt phi_j [nElGauss]. 
vector<vector<double>> el_phi_j;    

/// @brief Int-pt gPhi_jj [nElGauss]. 
vector<vector<double>> el_gPhi_jj;        

/// @brief Int-pt gPhi_ij [nElGauss]. 
vector<vector<double>> el_gPhi_ij;  

/// @brief Int-pt gPhi_ii [nElGauss].
vector<vector<double>> el_gPhi_ii;       

/// @brief Int-pt concentration [nElGauss].
vector<vector<double>> elCon; 

/// @brief Pointer of int-pt concentraiton [nElGauss]. 
vector<vector<double>>* elCon_ptr;

/// @brief nodal values of phi [nTotNodes]
vector<double> nod_gPhi;

/// @brief nodal values of phi [nTotNodes]
vector<double> nod_gPhi_HAGB;

/// @brief nodal values of phi [nTotNodes]
vector<double> nod_gPhi_LAGB;

/// @brief nodal values of phi_j [nTotNodes]
vector<double> nod_phi_j;  

/// @brief nodal values of gPhi_jj [nTotNodes]
vector<double> nod_gPhi_jj;      

/// @brief nodal values of gPhi_ij [nTotNodes]
vector<double> nod_gPhi_ij; 

/// @brief nodal values of gPhi_ii [nTotNodes]
vector<double> nod_gPhi_ii;    

/// @brief nodal values of hydrostatic stress [nTotNodes]
vector<double> nod_sigma_h; 

/// @brief nodal values of normalized dislocation density [nTotNodes]
vector<double> nod_rho; 

/// @brief Element concentration (temperature) dofs. In this case, it is identical to `elemNodeConn`.
vector<vector<int>> elemConDof;    

/// @brief Variant for returning elCapMatx. 
T_ElStiffMatx elCapMatxVariant;   

};
#endif