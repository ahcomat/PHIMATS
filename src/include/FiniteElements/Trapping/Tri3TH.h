/**
 * @file Tri3TH.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing Tri3TH elements for grain boundary trapping. 
 * @date 2024-07-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef TRI3TH_H
#define TRI3TH_H

#include "BaseElemTrap.h"
#include "Nodes.h"

class Tri3TH: public BaseElemTrap{

public:

Tri3TH(H5IO &H5File_in, Nodes &Nodes, int iSet);   

~Tri3TH() override ;

/**
 * @brief Initializes the data for the shape functions.
 * 
 */
void InitShapeFunc();

/**
 * @brief Returns the shape functions of int-pts in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return RowVecd3
 */
RowVecd3 CalcShapeFunc(double xi, double eta);

/**
 * @brief Returns the shape function derivatives of int-pts in natural coordinates.
 * 
 * @param xi 
 * @param eta 
 * @return Matd2x3 
 */
Matd2x3 CalcShapeFuncDeriv(double xi, double eta);

/**
 * @brief Initializes the data for all elements: `elemNodCoord`, `gaussPtCart`, `BMat`, `BuMat`
 *  and `intPtVol`.
 * 
 * @param Nodes 
 */
void InitializeElements(Nodes& Nodes, H5IO &H5File_in);

/**
 * @brief Get the cartesian coordinates of gauss points `N_i x_ij`.
 * 
 * @param elCoord 
 * @param sFunc 
 * @return RowVecd2 
 */
RowVecd2 getGaussCart(RowVecd3& sFunc, Matd3x2& elCoord);

/**
 * @brief Calculates `intPtVol`, `BMat` and `BuMat`.
 * 
 * @param elNodCoord 
 * @param sFuncDeriv 
 * @param intVol 
 * @param cartDeriv 
 * @param strainMat 
 */
void CalcCartDeriv(Matd3x2& elNodCoord, Matd2x3& sFuncDeriv, const double& wt, double& intVol, Matd2x3& cartDeriv);

/**
 * @brief Get the int-pt coordinates
 * 
 * @param glIntPtCoords 
 */
void getInPtCoords(T_nodStres& glIntPtCoords) override;

/**
 * @brief Evaluates the gradients of phi at the int-points and maps them to the nodes. 
 * 
 * @param nodGrad
 * @param nodCount 
 * 
 * @todo Move to BaseElemTrap.
 */
void CalcGrad(T_nodStres& nodGrad, vector<double>& nodCount, double* nodLapPhi) override;

/**
 * @brief Calculates the element stiffness and capacitance matrix for all elements.
 * 
 * @param mat Material
 * @param T Current temperature
 */
void CalcElemStiffMatx(BaseTrapping* mat, const double T) override;

/**
 * @brief Updates the element stiffness matrix for new temperature.
 * 
 * @param mat 
 * @param T 
 */
void UpdateElemStiffMatx(BaseTrapping* mat, const double T) override;

/**
 * @brief Calculates the volume averaged concentration
 * 
 * @param globalBuffer Solution vector
 * @return Volume averaged concentration
 */
double CalcAvCon(const double* globalBuffer) override;

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
void CalcFlux(BaseTrapping* mat, const double* globalBuffer, T_nodStres& nodFlux, vector<double>& nodCount, const double T) override;

private:

const double R = 8.31446261815324; /// @brief Universal gas constant [J/mol.K]

const vector<double> wts{1.0/6.0, 1.0/6.0, 1.0/6};  /// @brief Weights of the gauss points [nElGauss].

vector<RowVecd3> shapeFunc; /// @brief Values of the shape functions at integration points in natural coordinates [nElNodes].
vector<Matd2x3> shapeFuncDeriv; /// @brief Values of the shape function derivatives at integration points in natural coordinates [nElDim, nElNodes]. 

vector<Matd3x2> elemNodCoord;   /// @brief Node Coordinates [nElDim, nElNodes]. 

vector<vector<RowVecd2>> gaussPtCart;  /// @brief Cartesian coordinates of Gauss points for all elements [nElDim]. 

vector<vector<ColVecd2>> elFlux;         /// @brief Int-pt flux [nElStres]

vector<vector<double>> el_gPhi;          /// @brief Int-pt phi [nElStres]. 

vector<vector<double>> el_martensite;    /// @brief Int-pt martensite [nElStres]. 
vector<vector<double>> el_gPhiMM;        /// @brief Int-pt gPhiMM [nElStres]. 
vector<vector<double>> el_gPhifM;        /// @brief Int-pt gPhifM [nElStres]. 
vector<vector<double>> el_gPhiff;        /// @brief Int-pt gPhiff [nElStres]. 

vector<double> nod_gPhi;        /// @brief nodal values of phi [nTotNodes]

vector<double> nod_martensite;  /// @brief nodal values of martensite [nTotNodes]
vector<double> nod_gPhiMM;      /// @brief nodal values of gPhiMM [nTotNodes]
vector<double> nod_gPhifM;      /// @brief nodal values of gPhifM [nTotNodes]
vector<double> nod_gPhiff;      /// @brief nodal values of gPhiff [nTotNodes]

vector<vector<Matd2x3>> BMat;   /// @brief Derivatives (scalar) matrix [nElDim, nElNodes].

vector<Matd3x3> elStiffMatx;    /// @brief Element stiffness matrix [nElDispDofs, nElDispDofs].    
vector<Matd3x3> elCapMatx;      /// @brief Element capacitance matrix [nElDispDofs, nElDispDofs].   

};
#endif