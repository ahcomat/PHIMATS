#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "Materials/Trapping/TrapPhase.h"
#include "FiniteElements/Trapping/Tri3TH.h"
#include "Models/TrappingModel.h"
#include "Solvers/LinearSolver.h"
#include "Logger.h"

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

	// Model name, same as `Simul`
	string SimulName = "300_eq_f4_k1_dis";

	// Logger object for handling terminal user interface
	Logger logger(PETSC_COMM_WORLD, SimulName+".log");
	// Start timer
	logger.StartTimer();
	// Print intro message
	logger.IntroMessage();

	// Initialize I/O files

    // Mesh file
    const string meshFileName = "RVE.mesh.hdf5";
    H5IO meshH5File(meshFileName, logger);
    // RVE file
    const string rveFile = "RVE.rve.hdf5";
	H5IO H5File_rve(rveFile, logger);

    PhysicsIO diff(SimulName + ".diff", logger);

	// Materials vector
	vector<BaseTrapping*> matTVec;
	matTVec.push_back(new TrapPhase("2D", *diff.in, 1, logger, true));

	// Nodes
	Nodes Nodes;
	Nodes.ReadNodes(*diff.in, meshH5File);

	// Elements vector
	vector<BaseElemTrap*> diffElemVec;
	diffElemVec.push_back(new Tri3TH(*diff.in, meshH5File, Nodes, 1, logger, &H5File_rve));

	// Initialize the system -----------

	// Model
	TrappingModel model(diffElemVec, *diff.in, logger);

    // Sets initial uniform concentration
  	model.setUniformCon(0.018330);

	// Calculate the stiffness matrix of all elements
	model.CalcElemStiffMatx(diffElemVec, matTVec);
	
	// Initialize boundary conditions
	model.InitializeBC(*diff.in);

	// Assemble global stiffness matrix
	model.Assemble(diffElemVec);

	// Solver: "DIRECT" solver much faster for charging/permeation simulations
	LinearSolver linearSolver(model.getK(), logger, "DIRECT");

	// Write initial conditions
	model.WriteAvCon(diffElemVec, *diff.out, 0);
	model.WriteOut( *diff.out, to_string(0));
    model.WriteTrapCon( *diff.out, to_string(0));
    model.WriteAvTrapCon(diffElemVec, *diff.out, 0);
    // model.CalcFlux(diffElemVec, matTVec);
	// model.WriteFlux(*diff.out, to_string(0));

	// Solver loop -----------

	logger.LoopMessage();

	// Number of steps
	int nSteps = model.get_nSteps();
	// Frequency of full-field output
	int tOut = 1;

	for (int iStep=1; iStep<nSteps+1; iStep++){

		// Print step increment
		logger.StepIncrement(iStep);

        // Calculates integration point concentration
        model.CalcElCon(diffElemVec);

        // Update the carbide-trapped concentration
        model.CalcPrecipitaitonCon(diffElemVec, matTVec);

        // Update the trapped load RHS (must be before setBC!!)
        model.UpdateTrapRHS(diffElemVec);

		// Apply BCs
		model.setBC();
		
		// Solve
		linearSolver.Solve(model.getX(), model.getF());

        // // Calculate the flux 
		// model.CalcFlux(diffElemVec, matTVec);

		// Write averge concetration (every step)
		model.WriteAvCon(diffElemVec, *diff.out, iStep);

        // Write averge trapped concetration
        model.WriteAvTrapCon(diffElemVec, *diff.out, iStep);

        // Map trapped concentration to the nodes
        model.CalcNodTrapCon(diffElemVec);

		// Write field output
		if (!(iStep%tOut)){
			
			// Mobile concentration
			model.WriteOut( *diff.out, to_string(iStep));

            // model.WriteFlux(*diff.out, to_string(iStep));

            // Trapped concentration
            model.WriteTrapCon(*diff.out, to_string(iStep));

			logger.FieldOutput(iStep);
		}

	}

	logger.ExitMessage();

	// Deallocate
	for (auto* elem : diffElemVec) { // Elements
	delete elem;
	}

	for (auto* mat : matTVec) { // Material
	delete mat;
	}
    
    }

    // Finalize PETSc -----------
    
    PetscFinalize();

	return 0;
}