/**
 * @file BaseElemTrap.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for trapping elements. 
 * @date 2024-06-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASEELEMTRAP_H
#define BASEELEMTRAP_H

#include "FiniteElements/BaseElements.h"
#include "Materials/Trapping/BaseTrapping.h"

class BaseElemTrap: public BaseElements{

public:

BaseElemTrap(int nElDim, int nElNodes, int nElConDofs, int nElGauss)
    : nElDim(nElDim), nElNodes(nElNodes), nElConDofs(nElConDofs), nElGauss(nElGauss) {};

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
virtual void CalcElemStiffMatx(BaseTrapping* mat, const double T) = 0;

/**
 * @brief Calculates the volume averaged concentration
 * 
 * @param globalBuffer Solution vector
 * @return double Volume averaged concentration
 */
virtual double CalcAvCon(const double* globalBuffer) = 0;

/**
 * @brief Evaluates the gradients of scalar field at the int-points and maps them to the nodes. 
 * 
 * @param nodGrad
 * @param nodCount 
 */
virtual void CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapSigmaH) = 0;

/**
 * @brief Return const reference to the vector of element capacitance matrix c_ii.
 * 
 * @return const vector<T_ElStiffMatx>& 
 */
const T_ElStiffMatx& getElCapMatx() const { return elCapMatxVariant; }

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
virtual void CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, T_nodStres& intPtFlux, vector<double>& nodCount, const double T) = 0;

protected:

const int nElDim;              /// @brief Spatial dimensions of the element.
const int nElNodes;            /// @brief Number of nodes per element.
const int nElConDofs;          /// @brief Number of element concentration (temperature) dofs.
const int nElGauss;            /// @brief Number of gauss points.

double dt;                     /// @brief Time increment.     
int Trapping;                  /// @brief Flag for trapping type.

vector<vector<double>> intPtVol;    /// @brief Int-pt volume.       

vector<vector<int>> elemConDof;    /// @brief Element concentration (temperature) dofs. In this case, it is identical to `elemNodeConn`.

T_ElStiffMatx elCapMatxVariant;   /// @brief Variant for returning elCapMatx. 

};
#endif