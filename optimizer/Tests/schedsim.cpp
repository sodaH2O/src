#include "Simulation/OpalSimulation.h"
#include "Util/Types.h"

#include <string>
#include <iostream>
#include <map>

void test_callback(reqVars_t res);

int main(int argc, char** argv) {

    Param_t params;
    params.insert(std::pair<std::string, double>("QBUNCH", 1e-10));
    params.insert(std::pair<std::string, double>("NSLICE", 100));

    reqVars_t requestedVars;
    reqVarInfo_t v1 = {EVALUATE, 0.0, 1};
    reqVarInfo_t v2 = {EVALUATE, 0.0, 1};
    requestedVars.insert(std::pair< std::string, reqVarInfo_t >("E", v1));
    requestedVars.insert(std::pair< std::string, reqVarInfo_t >("dE", v2));

    SimulationScheduler *scheduler = SimulationScheduler::Instance();
    scheduler_m->dispatchSimulation("test-sim", params, requestedVars, test_callback);

    reqVars_t::iterator it;
    for(it = res.begin(); it != res.end(); it++)
        std::cout << it->first << " = " << it->second.value << std::endl;

    return 0;
}
