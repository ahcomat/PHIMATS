#include <iostream>
#include <stdexcept>
#include <cstdlib>  
#include "Materials/PFF/PFF.h"

PFF::PFF(string dimensions, H5IO& H5File, int iSet, Logger& logger)
    : BasePFF(dimensions, logger) {

    try {

        wc = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/PFF/wc");
        const_ell = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/PFF/const_ell");

    } catch (const std::exception& e){

        logger.log("Exception caught in PFF material constructor:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("    Critical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
        
    }
}

