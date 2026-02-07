#include <iostream>
#include <vector>

#include "petscsys.h"
#include "H5IO.h"
#include "Materials/Trapping/TrapPhase.h"
#include "Nodes.h"
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

	// Initialize PETSc (model object finalizes PETSc in the destructor)
	PetscInitialize(&argc, &argv, NULL, NULL);
    
    {

	// Read inputs -----------

	// Model name, same as `Simul`
	string SimulName = "ChargingHT990";

	// Logger object for handling terminal user interface
	Logger logger(PETSC_COMM_WORLD, SimulName+".log");
	// Start timer
	logger.StartTimer();
	// Print intro message
	logger.IntroMessage();

	// Initialize I/O hdf5 files

	// Mesh file
    const string meshFileName = SimulName + ".mesh.hdf5";
    H5IO meshH5File(meshFileName, logger);
    // RVE file
    const string rveFile = "../HT990.rve.hdf5";
	H5IO H5File_rve(rveFile, logger);

	PhysicsIO diff(SimulName + ".diff", logger);

	// Material vector
	vector<BaseTrapping*> matTVec;
	matTVec.push_back(new TrapPhase("2D", *diff.in, 1, logger));

	// Nodes
	Nodes Nodes;
	Nodes.ReadNodes(*diff.in, meshH5File);

	// Elements vector
	vector<BaseElemTrap*> Tri3THElemVec;
	Tri3THElemVec.push_back(new Tri3TH(*diff.in, meshH5File, Nodes, 1, logger, &H5File_rve));

	// Initialize the system -----------

	// Model
	TrappingModel model(Tri3THElemVec, *diff.in, logger);

	// Calculate the initial stiffness matrix
	model.CalcElemStiffMatx(Tri3THElemVec, matTVec);

	// Initialize boundary conditions
	model.InitializeBC(*diff.in);

	// Assemble global stiffness matrix
	model.Assemble(Tri3THElemVec);

	// Solver: "DIRECT" solver much faster for charging/permeation simulations
	LinearSolver linearSolver(model.getK(), logger, "DIRECT");

	// Write initial conditions
	model.WriteAvCon(Tri3THElemVec, *diff.out, 0);    // Volume average concentration 
	model.WriteOut(*diff.out, to_string(0));		   // Concentration field 

	// Solver loop -----------

	logger.LoopMessage();

	// Sets BC and updates the RHS. Required even for zero flux BCs
	model.setBC();

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

		// Write average concetration (every step)
		model.WriteAvCon(Tri3THElemVec, *diff.out, iStep);

		// Write field output (in this case, every step)
		if (!(iStep%tOut)){
			
			// Full field output
			model.WriteOut(*diff.out, to_string(iStep));
			logger.FieldOutput(iStep);

		}
	}

	logger.ExitMessage();

	// Deallocate -----------

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