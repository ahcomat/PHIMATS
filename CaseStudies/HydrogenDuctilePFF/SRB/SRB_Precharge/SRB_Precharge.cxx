#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "Solvers/LinearSolver.h"
#include "FiniteElements/Trapping/Quad4THAxi.h"
#include "Materials/Trapping/MechTrap.h"
#include "Models/TrappingModel.h"

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
    {
        // Read inputs -----------

        string SimulName = "SRB_Precharge";

        // Logger object for handling terminal user interface
        Logger logger(PETSC_COMM_WORLD, SimulName+".log");
        logger.StartTimer();
        logger.IntroMessage();

        // Initialize I/O files  -------------------------

        const string meshFileName = "../SmoothAxi.mesh.hdf5";
        H5IO meshH5File(meshFileName, logger);

        PhysicsIO diff(SimulName + ".diff", logger);

        // Nodes  -------------------------

        Nodes Nodes;
        Nodes.ReadNodes(*diff.in, meshH5File);

        // Material  -------------------------

        vector<BaseTrapping*> diffMatVec;
        diffMatVec.push_back(new MechTrap("2D", *diff.in, 1, logger));
        
        // Elements  -------------------------

        // Finite elements: Trapping element type
        vector<BaseElemTrap*> diffElemVec;
        diffElemVec.push_back(new Quad4THAxi(*diff.in, meshH5File, Nodes, 1, logger));

        // ----------- Initialize the system -----------

        // Diffusion model -----------
        
        TrappingModel modelDiff(diffElemVec, *diff.in, logger);

        // Sets initial uniform concentration
        modelDiff.setUniformCon(0);
        
        // Calculates integration point concentration
        modelDiff.CalcElCon(diffElemVec);

        // // Read initial nodal stresses from the input stress file.
        // modelDiff.ReadNodalStress(diffElemVec, *mech.out, 0);

        // Calculate the stiffness matrix of all elements
        modelDiff.CalcElemStiffMatx(diffElemVec, diffMatVec);

        // Dirichlet BCs (has to use this function even if only zero flux BCs are imposed)
        modelDiff.InitializeBC(*diff.in);

        // Assemble global stiffness matrix
        modelDiff.Assemble(diffElemVec);

        // Output Initial values
        modelDiff.WriteAvCon(diffElemVec, *diff.out, 0);
        modelDiff.WriteOut(*diff.out, to_string(0));

        // Solver: Linear with direct solver
        LinearSolver diffLinearSolver(modelDiff.getK(), logger);

        // Solver loop -----------

        logger.LoopMessage();

        int nSteps = modelDiff.get_nSteps();

        for (int iStep=1; iStep<nSteps+1; iStep++){

            logger.StepIncrement(iStep);

            //-------------------------------------------------------//

            // Diffusion -----------
            
            // The stiffness matrix and KSP have to be updated for the new stress
            modelDiff.CalcElemStiffMatx(diffElemVec, diffMatVec);
            modelDiff.Assemble(diffElemVec);
            diffLinearSolver.UpdateKSP(modelDiff.getK());

            // Sets BC and updates the RHS. Required even for zero flux BCs
            modelDiff.setBC();

            // Solve diffuison
            diffLinearSolver.Solve(modelDiff.getX(), modelDiff.getF());

            // Calculates integration point concentration
            modelDiff.CalcElCon(diffElemVec);

            // Write field output 
            modelDiff.WriteOut(*diff.out, to_string(iStep));
            modelDiff.WriteAvCon(diffElemVec, *diff.out, iStep);

            //-------------------------------------------------------//

            logger.FieldOutput(iStep);

        }

        logger.ExitMessage();

        // Deallocate

        // Elements -----------

        for (auto* elem : diffElemVec) delete elem;

        // Materials -----------

        for (auto* mat : diffMatVec) delete mat;

    }
  PetscFinalize();

  return 0;
}