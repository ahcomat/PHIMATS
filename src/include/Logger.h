/**
 * @file Logger.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief Class for managing application logging.
 * 
 * Provides functionality for logging messages with various severity levels (INFO, WARNING, ERROR).      
 * 
 * @date 2025-01-12
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

#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream>
#include <ctime>

#include <petscsys.h>  // PETSc header for parallel utilities

using namespace std;


class Logger {

public:

Logger(const std::string& fileName = "", MPI_Comm comm = PETSC_COMM_WORLD);
~Logger();

string applyColor(const std::string& level);

string stripAnsiCodes(const std::string& input);

void log(const std::string& message, const std::string& level = "INFO", bool includeTimestamp=true);

/**
 * @brief Show PhiMATS intro message. 
 * 
 */
void IntroMessage();

/**
 * @brief Show entering loop message
 * 
 */
void LoopMessage();

/**
 * @brief Show message for step increment. 
 * 
 */
void StepIncrement(const int& iStep);

/**
 * @brief Show message for step increment. 
 * 
 */
void FieldOutput(const int& iStep);

void StartTimer();

void ExitMessage();

private:

std::time_t startTime;
std::time_t endTime;

std::ofstream logFile;
bool logToFile;
int rank;    // MPI rank
int size;    // MPI size

// Get the current timestamp
std::string getCurrentTime() const;

};
#endif 
