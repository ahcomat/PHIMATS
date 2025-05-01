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

#include <mpi.h>  

using namespace std;


class Logger {

public:

Logger(MPI_Comm comm, const std::string& fileName = "");
~Logger();

string applyColor(const std::string& level);

string stripAnsiCodes(const std::string& input);

void log(const std::string& message, const std::string& level = "INFO", bool includeTimestamp=true);

/**
 * @brief Print PHIMATS intro message to terminal. 
 * 
 */
void IntroMessage();

/**
 * @brief Print entering loop message to terminal.
 * 
 */
void LoopMessage();

/**
 * @brief Print message for step increment to terminal.
 * 
 * @param iStep Step number. ks
 */
void StepIncrement(const int& iStep);

/**
 * @brief Print message for full field increment to terminal. 
 * 
 */
void FieldOutput(const int& iStep);

/**
 * @brief Start timer for simulations. 
 * 
 */
void StartTimer();

/**
 * @brief Print simulation completion message and total execution time to terminal. 
 * 
 */
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
