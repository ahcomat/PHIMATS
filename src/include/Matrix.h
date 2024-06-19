/**
 * @file Matrix.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief This file includes custom type aliases based on the library `Eigen`. 
 *        It also contains `std::variant` return types based on these aliases 
 *        to be used for derived classes.          
 * 
 * @date 2024-05-18
 * 
 * @copyright Copyright (c) 2024
 *  
 */

#include"Eigen/Dense"
#include<variant>
#include<vector>

using namespace std;

// Type aliases for `Eigen` vectors and matrices

typedef Eigen::Vector<double, 2> ColVecd2;    
typedef Eigen::Vector<double, 3> ColVecd3;
typedef Eigen::Vector<double, 4> ColVecd4;    
typedef Eigen::Vector<double, 6> ColVecd6;                 
typedef Eigen::Vector<double, 3> ColVecd3;                  
typedef Eigen::Vector<double, 8> ColVecd8;        
typedef Eigen::Vector<double, 24> ColVecd24;        

typedef Eigen::RowVector<double, 2> RowVecd2; 
typedef Eigen::RowVector<double, 3> RowVecd3;             
typedef Eigen::RowVector<double, 4> RowVecd4;
typedef Eigen::RowVector<double, 8> RowVecd8; 

typedef Eigen::Matrix<double, 2, 2> Matd2x2;
typedef Eigen::Matrix<double, 3, 3> Matd3x3; 
typedef Eigen::Matrix<double, 4, 4> Matd4x4; 
typedef Eigen::Matrix<double, 6, 6> Matd6x6;
typedef Eigen::Matrix<double, 8, 8> Matd8x8;
typedef Eigen::Matrix<double, 24, 24> Matd24x24; 

typedef Eigen::Matrix<double, 4, 2> Matd4x2;           
typedef Eigen::Matrix<double, 2, 4> Matd2x4;

typedef Eigen::Matrix<double, 2, 3> Matd2x3;            
typedef Eigen::Matrix<double, 3, 2> Matd3x2;   

typedef Eigen::Matrix<double, 3, 6> Matd3x6;     
typedef Eigen::Matrix<double, 3, 8> Matd3x8;
typedef Eigen::Matrix<double, 8, 3> Matd8x3; 

typedef Eigen::Matrix<double, 6, 24> Matd6x24; 


// Aliases based on `std::variant`.
// NOTE: This is type alias, not variable definition.


// For some reason this doesn't work O.o and makes extremely unexpected behavior 
// /**
//  * @brief Variants `KMatx`.
//  * 
//  */
// using T_KMatx = std::variant<Matd2x2, Matd3x3>;

/**
 * @brief Variants `DMatx`.
 * 
 */
using T_DMatx = std::variant<Matd2x2, Matd3x3, Matd6x6>;

/**
 * @brief Variants `elStiffMatx`. It is a pointer because we need this to avoid copy in `Elements::getElStiffMatx`.
 * 
 */
using T_ElStiffMatx = std::variant<vector<Matd4x4>*, vector<Matd6x6>*, vector<Matd8x8>*, vector<Matd24x24>*>;

/**
 * @brief Variants for vector (tensor in Voigt notation) quantities.
 * 
 */
using T_nodStres = std::variant<vector<ColVecd2>, vector<ColVecd3>, vector<ColVecd6>>;