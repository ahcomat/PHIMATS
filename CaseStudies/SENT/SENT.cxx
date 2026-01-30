#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "FiniteElements/Mechanics/Quad4.h"
#include "FiniteElements/PFF/Quad4PFF.h"
#include "Materials/Mechanics/IsoHard.h"
#include "Materials/PFF/PFF.h"
#include "Models/MechModel.h"
#include "Models/PFFModel.h"
#include "Solvers/LinearSolver.h"

struct PhysicsIO {
    H5IO* in;
    H5IO* out;
    string modelName;

    PhysicsIO(const string& name, Logger& logger) : modelName(name) {
        in  = new H5IO(name + ".in.hdf5", logger);
        out = new H5IO(name + ".out.hdf5", logger);
    }
    ~PhysicsIO() { delete in; delete out; }
};

using namespace std;

int main(int argc, char **argv){

    // Initialize PETSc (always here!)
    PetscInitialize(&argc, &argv, NULL, NULL);

    // Read inputs -----------

    string SimulName = "SENT_EPTH";
    
    // Logger object for handling terminal user interface
    Logger logger(PETSC_COMM_WORLD, SimulName+".log");
    logger.StartTimer();
    logger.IntroMessage();

  // Initialize I/O files  -------------------------

    const string meshFileName = SimulName + ".mesh.hdf5";
    H5IO meshH5File(meshFileName, logger);

    PhysicsIO mech(SimulName + ".mech", logger);
    PhysicsIO pff(SimulName + ".pff", logger);

    // Nodes  -------------------------

    Nodes Nodes;
    Nodes.ReadNodes(*mech.in, meshH5File);

    // Material  -------------------------

    vector<BaseMechanics*> mechMatVec;
    mechMatVec.push_back(new IsoHard("2D", *mech.in, 1, logger));

    vector<BasePFF*> pffMatVec;
    pffMatVec.push_back(new PFF("2D", *pff.in, 1, logger));
    
    // Elements  -------------------------

    vector<BaseElemMech*> mechElemVec;
    mechElemVec.push_back(new Quad4(*mech.in, meshH5File, Nodes, 1, "ElastoPlastic", logger));

    vector<BaseElemPFF*> pffElemVec;
    pffElemVec.push_back(new Quad4PFF(*pff.in, meshH5File, Nodes, 1, pffMatVec[0]->get_const_ell(), logger));

    // ----------- Initialize the system -----------

    // Mechanics model -----------
    MechModel modelMechanics(mechElemVec, *mech.in, logger, "DIRECT", 1);
    
    // Dirichlet BCs
    modelMechanics.InitializeDirichBC(*mech.in);

    // Calculate element stiffness matrix
    modelMechanics.CalcElemStiffMatx(mechElemVec, mechMatVec);

    // Assemble global stiffness matrix
    modelMechanics.Assemble(mechElemVec);

    // Write initial conditions
    modelMechanics.CalcNodVals(mechElemVec);
    modelMechanics.WriteOut(mechElemVec, *mech.out, to_string(0));

    // PFF model -----------
    PFFModel modelPFF(pffElemVec, *pff.in, logger);

    // Sets a constant critical work energy density
    modelPFF.set_const_wc(pffElemVec, pffMatVec);

    // Calculate gPhi_d
    modelPFF.Calc_gPhi_d(pffElemVec);

    // Write initial conditions
    modelPFF.WriteOut(*pff.out, to_string(0));

    // Linear solver for PFF
    LinearSolver pffLinearSolver(modelPFF.getK(), logger);

    // Solver loop -----------

        logger.LoopMessage();

    // Solve the system
    int nSteps = modelMechanics.get_nSteps();

    for (int iStep=1; iStep<nSteps+1; iStep++){

        logger.StepIncrement(iStep);

        // Mechanics -----------

        // Apply BCs
        modelMechanics.setDirichBC();

        // Solve mechanics
        modelMechanics.SolveSNES(mechElemVec, mechMatVec, iStep, pffElemVec);
        
        // Map stresses and strains to nodes
        modelMechanics.CalcNodVals(mechElemVec);

        // Write field output
        modelMechanics.WriteOut(mechElemVec, *mech.out, to_string(iStep));

        // PFF -----------

        // Spectral decomposition of elastic strain energy density
        modelPFF.CalcPsiSpectral(pffElemVec, mechElemVec, mechMatVec);

        // *******************

        // Chose a crack driving force by uncommenting ONLY one of the options below:

        // // Brittle crack driving force
        // modelPFF.CalcDrivForcB(pffElemVec);   

        // // Elastoplastic crack driving force
        // modelPFF.CalcDrivForcEP(pffElemVec, mechElemVec);

        // Elastoplastic crack driving force with a threshold
        modelPFF.CalcDrivForcEP_TH(pffElemVec, mechElemVec, 0.25);       

        // *******************
        
        // Build element stiffness matrix
        modelPFF.CalcElemStiffMatx(pffElemVec);
        // Assemble global stiffness matrix
        modelPFF.Assemble(pffElemVec);
        // Calculate RHS
        modelPFF.CalcFH(pffElemVec);

        // Solve PFF (linear)
        pffLinearSolver.Solve(modelPFF.getX(), modelPFF.getF());
        // Calculate gPhi_d
        modelPFF.Calc_gPhi_d(pffElemVec);
        // Write field output
        modelPFF.WriteOut(*pff.out, to_string(iStep));

        //     -----------

        logger.FieldOutput(iStep);

    }

    logger.ExitMessage();

    // Deallocate

    // Elements -----------
    for (auto* elem : mechElemVec) { // Elements
        delete elem;
    }

        for (auto* elem : pffElemVec) { // Elements
        delete elem;
    }

    // Materials -----------
    for (auto* mat : mechMatVec) { // Mech material
        delete mat;
    }

    for (auto* mat : pffMatVec) { // Pff material
        delete mat;
    }

    return 0;
}