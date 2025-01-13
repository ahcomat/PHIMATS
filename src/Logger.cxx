#include "Logger.h"
#include "Version.h"
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

string Logger::LevelToString(LogLevel level) {

    switch (level) {
        case INFO: return "INFO";
        case WARNING: return "WARNING";
        case ERROR: return "ERROR";
        default: return "UNKNOWN";
    }
    
}

void Logger::log(const std::string& message, const std::string& level, bool includeTimestamp) {
    std::string fullMessage;
    if (includeTimestamp) {
        fullMessage = "[" + getCurrentTime() + "] [" + level + "] " + message;
    } else {
        fullMessage = message;
    }

    // Only rank 0 logs to console and file
    if (rank == 0) {
        // Console output: Include ANSI codes for styling
        std::string consoleMessage = fullMessage;
        std::cout << consoleMessage << std::endl;

        // File output: Exclude ANSI escape sequences (clean plain text)
        if (logToFile) {
            logFile << fullMessage << std::endl;
        }
    }
}

void Logger::IntroMessage() {
    if (rank == 0) {
        const std::string redPhi = "\033[31mφ\033[0m";  // Red phi (phi) for terminal output
        const std::string plainPhi = "φ";              // Plain phi for log file output

        log("", "", false);
        log("*       Phase-field Multiphysics Materials Simulator (" + redPhi + "MATS) ", "", false);
        log("*       Version: " + std::string(VERSION_STRING), "", false);
        log("*       Release date: " + std::string(VERSION_DATE), "", false);
        // log("*       Website: " + std::string(PROJECT_WEBSITE), "", false);
        // log("*       For citation, please use: " + std::string(PROJECT_CITATION), "", false);
        log("", "", false);
    }
}

void Logger::LoopMessage(){

    if (rank == 0) {
        log("", "", false);
        log("   >>> Entering the solver loop <<<", LevelToString(INFO));
        log("", "", false);
    }
}

void Logger::StepIncrement(const int& iStep){

        if (rank == 0) {
        log("", "", false);
        log("   Increment " + std::to_string(iStep), LevelToString(INFO));
        log("", "", false);
    }
}

void Logger::FieldOutput(const int& iStep){

        if (rank == 0) {
        log("", "", false);
        log("   Field output for Step_" + std::to_string(iStep), LevelToString(INFO));
        log("", "", false);
    }
}

