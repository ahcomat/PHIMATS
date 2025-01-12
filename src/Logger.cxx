#include "Logger.h"
#include <iostream>
#include <ctime>
#include <stdexcept>

std::string Logger::getCurrentTime() const {
    std::time_t now = std::time(nullptr);
    char buf[80];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return std::string(buf);
}

Logger::Logger(const std::string& fileName, MPI_Comm comm)
    : logToFile(!fileName.empty()) {
    // Get rank and size from PETSc communicator
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Open log file if specified and rank is 0
    if (logToFile && rank == 0) {
        logFile.open(fileName, std::ios::out | std::ios::app);
        if (!logFile.is_open()) {
            throw std::runtime_error("Failed to open log file: " + fileName);
        }
    }
}

Logger::~Logger() {
    if (logToFile && logFile.is_open() && rank == 0) {
        logFile.close();
    }
}

void Logger::log(const std::string& message, const std::string& level) {
    std::string fullMessage = "[" + getCurrentTime() + "] [Rank " + std::to_string(rank) +
                              "] [" + level + "] " + message;

    // Only rank 0 logs to console and file
    if (rank == 0) {
        std::cout << fullMessage << std::endl;
        if (logToFile) {
            logFile << fullMessage << std::endl;
        }
    }
}

void Logger::showIntroMessage() {
    if (rank == 0) {
        log("**********************************************");
        log("*       Welcome to My Parallel Simulation    *");
        log("*       Version: 1.0                         *");
        log("*       Author: Your Name                    *");
        log("*       Description: Finite Element Solver   *");
        log("**********************************************");
    }
}
