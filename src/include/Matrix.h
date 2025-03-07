/**
 * @file Matrix.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief This file includes custom type aliases based on the library `Eigen`. 
 *        It also contains `std::variant` return types based on these aliases 
 *        to be used for derived classes.          
 * 
 * @date 2024-05-18
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
typedef Eigen::RowVector<double, 6> RowVecd6;
typedef Eigen::RowVector<double, 8> RowVecd8; 

typedef Eigen::Matrix<double, 2, 2, Eigen::RowMajor> Matd2x2;
typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Matd3x3;  // `RowMajor` for non-symmetric element stiffness matrix.
typedef Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Matd4x4; 
typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> Matd6x6;
typedef Eigen::Matrix<double, 8, 8, Eigen::RowMajor> Matd8x8;
typedef Eigen::Matrix<double, 24, 24, Eigen::RowMajor> Matd24x24; 

typedef Eigen::Matrix<double, 4, 2> Matd4x2;           
typedef Eigen::Matrix<double, 2, 4> Matd2x4;

typedef Eigen::Matrix<double, 2, 3> Matd2x3;            
typedef Eigen::Matrix<double, 3, 2> Matd3x2;   

typedef Eigen::Matrix<double, 2, 6> Matd2x6;            
typedef Eigen::Matrix<double, 6, 2> Matd6x2;

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
using T_ElStiffMatx = std::variant<vector<Matd3x3>*, vector<Matd4x4>*, vector<Matd6x6>*, vector<Matd8x8>*, vector<Matd24x24>*>;

/**
 * @brief Variants for vector (tensor in Voigt notation) quantities.
 * 
 */
using T_nodStres = std::variant<vector<ColVecd2>, vector<ColVecd3>, vector<ColVecd4>, vector<ColVecd6>>;