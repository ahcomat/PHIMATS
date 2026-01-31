#include <iostream>
#include <stdexcept>
#include <cstdlib>  
#include "Materials/PFF/ChemoMech.h"

ChemoMech::ChemoMech(string dimensions, H5IO& H5File, int iSet, Logger& logger)
    : BasePFF(dimensions, logger) {

    try {

        wc = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/ChemoMech/wc");
        // const_ell = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/ChemoMech/const_ell");

        wc_min = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/ChemoMech/wc_min");
        beta = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/ChemoMech/beta");
        c_DTB = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/ChemoMech/c_DTB");
        c_crit = H5File.ReadScalar("Materials/Material_" + to_string(iSet)+"/ChemoMech/c_crit");

    } catch (const std::exception& e){

        logger.log("Exception caught in ChemoMech material constructor:\n", "ERROR", true);
        logger.log("    " + std::string(e.what()), "", false);
        logger.log("    Critical error encountered. Terminating!\n", "", false);
        exit(EXIT_FAILURE);
        
    }
}

