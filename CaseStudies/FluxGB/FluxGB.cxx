#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "Materials/Trapping/TrapGB.h"
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
	string SimulName = "FluxGB";

	// Logger object for handling terminal user interface
	Logger logger(PETSC_COMM_WORLD, SimulName+".log");
	// Start timer
	logger.StartTimer();
	// Print intro message
	logger.IntroMessage();

	// Initialize I/O files

    // Mesh file
    const string meshFileName = SimulName + ".mesh.hdf5";
    H5IO meshH5File(meshFileName, logger);
    // RVE file
    const string rveFile = "GB.rve.hdf5";
	H5IO H5File_rve(rveFile, logger);

    PhysicsIO diff(SimulName + ".diff", logger);

	// Materials vector
	vector<BaseTrapping*> matTVec;
	matTVec.push_back(new TrapGB("2D", *diff.in, 1, logger));

	// Nodes
	Nodes Nodes;
	Nodes.ReadNodes(*diff.in, meshH5File);

	// Elements vector
	vector<BaseElemTrap*> Tri3THElemVec;
	Tri3THElemVec.push_back(new Tri3TH(*diff.in, meshH5File, Nodes, 1, logger, &H5File_rve));

	// Initialize the system -----------

	// Model
	TrappingModel model(Tri3THElemVec, *diff.in, logger);

	// Calculate the stiffness matrix of all elements
	model.CalcElemStiffMatx(Tri3THElemVec, matTVec);
	
	// Initialize boundary conditions
	model.InitializeBC(*diff.in);

	// Assemble global stiffness matrix
	model.Assemble(Tri3THElemVec);

	// Solver: "DIRECT" solver much faster for charging/permeation simulations
	LinearSolver linearSolver(model.getK(), logger, "DIRECT");

	// Write initial conditions
	model.WriteAvCon(Tri3THElemVec, *diff.out, 0);
	model.WriteOut( *diff.out, to_string(0));
	model.CalcFlux(Tri3THElemVec, matTVec);
	model.WriteExitFlux(*diff.out, 0);
	model.WriteFlux(*diff.out, to_string(0));

	// Solver loop -----------

	logger.LoopMessage();

	// Number of steps
	int nSteps = model.get_nSteps();
	// Frequency of full-field output
	int tOut = 1;

	for (int iStep=1; iStep<nSteps+1; iStep++){

		// Print step increment
		logger.StepIncrement(iStep);

		// Apply BCs
		model.setBC();
		
		// Solve
		linearSolver.Solve(model.getX(), model.getF());

		// Calculate the flux 
		model.CalcFlux(Tri3THElemVec, matTVec);

		// Write averge concetration (every step)
		model.WriteAvCon(Tri3THElemVec, *diff.out, iStep);

		// Write exit side flux (every step)
		model.WriteExitFlux(*diff.out, iStep);

		// Write field output
		if (!(iStep%tOut)){
			
			// Flux
			model.WriteFlux(*diff.out, to_string(iStep));
			// Concentration
			model.WriteOut( *diff.out, to_string(iStep));

			logger.FieldOutput(iStep);
		}

	}

	logger.ExitMessage();

	// Deallocate
	for (auto* elem : Tri3THElemVec) { // Elements
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