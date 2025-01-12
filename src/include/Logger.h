/**
 * @file Logger.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for handling messages.          
 * 
 * @date 2025-01-12
 * 
 * @copyright Copyright (C) 2024 Abdelrahman Hussein
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

#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream>

#include <petscsys.h>  // PETSc header for parallel utilities

class Logger {

public:

Logger(const std::string& fileName = "", MPI_Comm comm = PETSC_COMM_WORLD);
~Logger();

void log(const std::string& message, const std::string& level = "INFO");
void showIntroMessage();

private:

std::ofstream logFile;
bool logToFile;
int rank;    // MPI rank
int size;    // MPI size

// Get the current timestamp
std::string getCurrentTime() const;

};
#endif 
