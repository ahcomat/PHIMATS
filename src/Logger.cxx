#include "Logger.h"
#include "Version.h"
#include <iostream>
#include <stdexcept>
#include <regex>

std::string Logger::getCurrentTime() const {
    
    std::time_t now = std::time(nullptr);
    char buf[80];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return std::string(buf);
}

Logger::Logger(MPI_Comm comm, const std::string& fileName)
    : logToFile(!fileName.empty()) {

    // Get rank and size from PETSc communicator
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Open log file if specified and rank is 0
    if (logToFile && rank == 0) {
        logFile.open(fileName, std::ios::out | std::ios::trunc);
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

string Logger::applyColor(const std::string& level) {

    if (level == "INFO") {
        return "\033[32m" + level  + "\033[0m"; // Green
    } else if (level == "WARNING") {
        return "\033[33m" + level  + "\033[0m"; // Yellow
    } else if (level == "ERROR" || level == "CRITICAL") {
        return "\033[31m" + level  + "\033[0m"; // Red
    }
    return level ; // Default: No color
    
}


std::string Logger::stripAnsiCodes(const std::string& input) {
    // Regex to match ANSI escape sequences
    static const std::regex ansiRegex(R"(\x1B\[[0-9;]*m)");
    return std::regex_replace(input, ansiRegex, "");
}


void Logger::log(const std::string& message, const std::string& level, bool includeTimestamp) {
    
    std::string fullMessage;

    if (includeTimestamp) {
        fullMessage = "[" + getCurrentTime() + "] [" + applyColor(level) + "] " + message;
    } else {
        fullMessage = message;
    }   

    // Only rank 0 logs to console and file
    if (rank == 0) {
        // Console output: Include ANSI codes for styling
        std::string consoleMessage = fullMessage;
        std::cout << consoleMessage << std::endl;

        // File output: Plain text without ANSI escape codes
        if (logToFile) {
            std::string plainMessage = stripAnsiCodes(fullMessage); // Strip ANSI codes
            logFile << plainMessage << std::endl;
        }
    }
}


void Logger::IntroMessage() {
    if (rank == 0) {
        const std::string redPhi = "\033[31mφ\033[0m";  // Red phi (phi) for terminal output
        const std::string plainPhi = "φ";              // Plain phi for log file output

        log("", "", false);
        log("****************************************************************************", "", false);
        log("*                                                                          *", "", false);
        log("*       Phase-field Multiphysics Materials Simulator (" + plainPhi + "\033[31mMATS\033[0m)               *", "", false);
        log("*       Version: " + std::string(VERSION_STRING) + "                                                    *", "", false);
        log("*       Release date: " + std::string(VERSION_DATE)+ "                                           *", "", false);
        log("*       Website: " + std::string(PROJECT_WEBSITE)+ "                              *", "", false);
        log("*       For citation, please use:                                          *", "", false);
        log("*       Abdelrahman Hussein, Finite Element Theory for PHIMATS, 2025       *","", false);
        log("*       DOI: 10.48550/ARXIV.2502.16283                                     *","", false);
        log("*                                                                          *", "", false);
        log("****************************************************************************", "", false);
        log("", "", false);
    }
}

void Logger::LoopMessage(){

    if (rank == 0) {
        log("   >>> Entering the solver loop <<<", "INFO");
        log("", "", false);
    }
}

void Logger::StepIncrement(const int& iStep){

        if (rank == 0) {
        log("   Increment " + std::to_string(iStep), "INFO");
        log("", "", false);
    }
}

void Logger::FieldOutput(const int& iStep){

        if (rank == 0) {
        log("   Field output for Step_" + std::to_string(iStep), "INFO");
        log("", "", false);
    }
}

void Logger::StartTimer() {
    startTime = std::time(nullptr); // Record the start time
}

void Logger::ExitMessage(){

    endTime = std::time(nullptr); // Record the current (or end) time

    if (rank == 0) {
        log("   Simulation completed successfully", "INFO");
        log("   Total run time is: "+to_string(std::difftime(endTime, startTime)/60.0)+" minutes." , "INFO");
        log("", "", false);
    }
}

