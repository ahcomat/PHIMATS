#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "FiniteElements/Mechanics/Quad4.h"
#include "FiniteElements/PFF/Quad4PFF.h"
#include "FiniteElements/Trapping/Quad4TH.h"
#include "Materials/Mechanics/IsoHard.h"
#include "Models/MechModel.h"
#include "Materials/PFF/ChemoMech.h"
#include "Materials/Trapping/MechTrap.h"
#include "Models/TrappingModel.h"
#include "Models/PFFModel.h"
#include "Solvers/LinearSolver.h"
#include <memory>

struct PhysicsIO {
    std::unique_ptr<H5IO> in;
    std::unique_ptr<H5IO> out;
    string modelName;

    PhysicsIO(const string& name, Logger& logger) : modelName(name) {
        in  = std::make_unique<H5IO>(name + ".in.hdf5", logger);
        out = std::make_unique<H5IO>(name + ".out.hdf5", logger);
    }
    
    ~PhysicsIO() = default;

};

using namespace std;

int main(int argc, char **argv){

    // Initialize PETSc (always here!)
    PetscInitialize(&argc, &argv, NULL, NULL);

    {

    // Read inputs -----------

    string SimulName = "CT_115MPa";

    // Logger object for handling terminal user interface
    Logger logger(PETSC_COMM_WORLD, SimulName+".log");
    logger.StartTimer();
    logger.IntroMessage();

    // Initialize I/O files  -------------------------

    const string meshFileName = "../CT.mesh.hdf5";
    H5IO meshH5File(meshFileName, logger);

    PhysicsIO mech(SimulName + ".mech", logger);
    PhysicsIO pff(SimulName + ".pff", logger);
    PhysicsIO diff(SimulName + ".diff", logger);

    // Nodes  -------------------------

    Nodes Nodes;
    Nodes.ReadNodes(*mech.in, meshH5File);

    // Material  -------------------------

    vector<BaseMechanics*> mechMatVec;
    mechMatVec.push_back(new IsoHard("2D", *mech.in, 1, logger));

    vector<BaseTrapping*> diffMatVec;
    diffMatVec.push_back(new MechTrap("2D", *diff.in, 1, logger));

    vector<BasePFF*> pffMatVec;
    pffMatVec.push_back(new ChemoMech("2D", *pff.in, 1, logger));
    
    // Elements  -------------------------

    vector<BaseElemMech*> mechElemVec;
    mechElemVec.push_back(new Quad4(*mech.in, meshH5File, Nodes, 1, "ElastoPlastic", logger));

    vector<BaseElemTrap*> diffElemVec;
    diffElemVec.push_back(new Quad4TH(*diff.in, meshH5File, Nodes, 1, logger));

    vector<BaseElemPFF*> pffElemVec;
    pffElemVec.push_back(new Quad4PFF(*pff.in, meshH5File, Nodes, 1, pffMatVec[0]->get_const_ell(), logger));

    // ----------- Initialize the system -----------

    // Mechanics model -----------

    MechModel modelMechanics(mechElemVec, *mech.in, logger, "DIRECT", 3);
    
    // Dirichlet BCs
    modelMechanics.InitializeDirichBC(*mech.in);

    // Calculate element stiffness matrix
    modelMechanics.CalcElemStiffMatx(mechElemVec, mechMatVec);

    // Assemble global stiffness matrix
    modelMechanics.Assemble(mechElemVec);

    // Write initial conditions
    modelMechanics.CalcNodVals(mechElemVec);
    modelMechanics.WriteOut(mechElemVec, *mech.out, to_string(0));

    // Diffusion model -----------
        
    TrappingModel modelDiff(diffElemVec, *diff.in, logger);

    // Sets initial uniform concentration
    modelDiff.setUniformCon(0);
    
    // Calculates integration point concentration
    modelDiff.CalcElCon(diffElemVec);

    // Read initial nodal stresses from the input stress file.
    modelDiff.ReadNodalStress(diffElemVec, *mech.out, 0);

    // Calculate the stiffness matrix of all elements
    modelDiff.CalcElemStiffMatx(diffElemVec, diffMatVec, pffElemVec);

    // Dirichlet BCs (has to use this function even if only zero flux BCs are imposed)
    modelDiff.InitializeBC(*diff.in);

    // Assemble global stiffness matrix
    modelDiff.Assemble(diffElemVec);

    // Output Initial values
    modelDiff.WriteAvCon(diffElemVec, *diff.out, 0);
    modelDiff.WriteOut(*diff.out, to_string(0));
    modelDiff.CalcFlux(diffElemVec, diffMatVec, pffElemVec);
    modelDiff.WriteFlux(*diff.out, to_string(0));

    // Solver: Linear with direct solver
    LinearSolver diffLinearSolver(modelDiff.getK(), logger);

    // PFF model -----------

    PFFModel modelPFF(pffElemVec, *pff.in, logger);

    // Sets a constant critical work energy density
    modelPFF.set_decay_wc(pffElemVec, pffMatVec, diffElemVec);

    // Calculate gPhi_d
    modelPFF.Calc_gPhi_d(pffElemVec);

    // Map int-pt values to nodes
    modelPFF.CalcNodVals(pffElemVec, mechElemVec);

    // Write initial conditions
    modelPFF.WriteOut(*pff.out, to_string(0));

    // Linear solver for PFF
    LinearSolver pffLinearSolver(modelPFF.getK(), logger);

    // Solver loop -----------

    logger.LoopMessage();

    int nSteps = modelMechanics.get_nSteps();

    for (int iStep=1; iStep<nSteps+1; iStep++){

        logger.StepIncrement(iStep);

        //-------------------------------------------------------//

        // Mechanics -----------

        // Apply BCs
        modelMechanics.setDirichBC();

        // Non-linear solver
        modelMechanics.SolveSNES(mechElemVec, mechMatVec, iStep, pffElemVec);
        
        // Map stresses and strains to nodes
        modelMechanics.CalcNodVals(mechElemVec);

        // Write field output
        modelMechanics.WriteOut(mechElemVec, *mech.out, to_string(iStep));

        //-------------------------------------------------------//

        // Diffusion -----------

        // Read stress for the step
        modelDiff.ReadNodalStress(diffElemVec, *mech.out, iStep);
        
        // The stiffness matrix and KSP have to be updated for the new stress
        modelDiff.CalcElemStiffMatx(diffElemVec, diffMatVec, pffElemVec);
        modelDiff.Assemble(diffElemVec);
        diffLinearSolver.UpdateKSP(modelDiff.getK());

        // Set equilibrium concentration boundary conditions
        modelDiff.CalcEquilibriumBC(diffElemVec, diffMatVec);

        // Sets BC and updates the RHS. Required even for zero flux BCs
        modelDiff.setBC();

        // Make the crack a hydrogen source. Comment out if sink. 
        modelDiff.CalcFsrc(diffElemVec, diffMatVec, pffElemVec);

        // Solve diffuison
        diffLinearSolver.Solve(modelDiff.getX(), modelDiff.getF());

        // Calculates integration point concentration
        modelDiff.CalcElCon(diffElemVec);

        // Write field output 
        modelDiff.WriteOut(*diff.out, to_string(iStep));
        modelDiff.WriteAvCon(diffElemVec, *diff.out, iStep);
        modelDiff.CalcFlux(diffElemVec, diffMatVec, pffElemVec);
        modelDiff.WriteFlux(*diff.out, to_string(iStep));

        //-------------------------------------------------------//

        // PFF -----------

        modelPFF.set_decay_wc(pffElemVec, pffMatVec, diffElemVec);

        // Spectral decomposition of elastic strain energy density
        modelPFF.CalcPsiSpectral(pffElemVec, mechElemVec, mechMatVec);

        // *******************     

        // Ductile crack driving force with a threshold
        modelPFF.CalcDrivForcHybridDuctile_TH(pffElemVec, mechElemVec, 1.0, 0.6); 

        // *******************
        
        // Build element stiffness matrix
        modelPFF.CalcElemStiffMatx(pffElemVec);
        // Assemble global stiffness matrix
        modelPFF.Assemble(pffElemVec);
        pffLinearSolver.UpdateKSP(modelPFF.getK());
        // Calculate RHS
        modelPFF.CalcFH(pffElemVec);

        // Solve PFF (linear)
        pffLinearSolver.Solve(modelPFF.getX(), modelPFF.getF());
        // Calculate gPhi_d
        modelPFF.Calc_gPhi_d(pffElemVec);
        // Map int-pt values to nodes
        modelPFF.CalcNodVals(pffElemVec, mechElemVec);
        // Write field output
        modelPFF.WriteOut(*pff.out, to_string(iStep));

        //     -----------

        logger.FieldOutput(iStep);

    }

    logger.ExitMessage();

    // Deallocate

    // Elements -----------
    for (auto* elem : mechElemVec) delete elem;
    for (auto* elem : diffElemVec) delete elem;
    for (auto* elem : pffElemVec) delete elem;

    // Materials -----------
    for (auto* mat : mechMatVec) delete mat;
    for (auto* mat : diffMatVec) delete mat;
    for (auto* mat : pffMatVec) delete mat;
    
    }

    // Finalize PETSc -----------
    
    PetscFinalize();

  return 0;
}