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

  	// Initialize PETSc 
	PetscInitialize(&argc, &argv, NULL, NULL);

  // Read inputs -----------

	// Model name, same as `Simul`
	string SimulName = "TDSHT990";

  	// Logger object for handling terminal user interface
	Logger logger(PETSC_COMM_WORLD, SimulName+".log");
	logger.StartTimer();
	logger.IntroMessage();

    // Initialize I/O hdf5 files
	// Mesh file
    const string meshFileName = SimulName + ".mesh.hdf5";
    H5IO meshH5File(meshFileName, logger);
    // RVE file
    const string rveFile = "../HT990.rve.hdf5";
	H5IO H5File_rve(rveFile, logger);

	PhysicsIO diff(SimulName + ".diff", logger);

    // Initial conditions for TDS simulation
    const string infileName2 = "../Charging/ChargingHT990.diff.out.hdf5";
	H5IO H5File_in2(infileName2, logger);

    // Material vector
    vector<BaseTrapping*> matTVec;
    matTVec.push_back(new TrapPhase("2D", *diff.in, 1, logger));

    // Nodes
    Nodes Nodes;
    Nodes.ReadNodes(*diff.in, meshH5File);
    
    // Elements
    vector<BaseElemTrap*> Tri3THElemVec;
    Tri3THElemVec.push_back(new Tri3TH(*diff.in, meshH5File, Nodes, 1, logger, &H5File_rve));

    // Initialize the sytem -----------

    // Model
    TrappingModel model(Tri3THElemVec, *diff.in, logger);
    
    // Calculate the stiffness matrix
    model.CalcElemStiffMatx(Tri3THElemVec, matTVec);

    // Initialize boundary conditions (must be called before `Assemble`!)
    model.InitializeBC(*diff.in);

    // Assemble global stiffness matrix
    model.Assemble(Tri3THElemVec);

    // Solver: "GMRES" solver is more efficient for large system TDS simulations because
    // it handles the frequent updates to the stiffness matrix more effectively
    LinearSolver linearSolver(model.getK(), logger, "GMRES");

    // Read initial data (Last step in the charing simulation)
    model.ReadInitialCon(H5File_in2, 20);

    // Write initial conditions
    model.WriteTemp(*diff.out, 0);
	model.WriteAvCon(Tri3THElemVec, *diff.out, 0);
	model.WriteOut( *diff.out, to_string(0));
  
  // Solver loop -----------

	logger.LoopMessage();

  //Solve the system
	int nSteps = model.get_nSteps();

  // Frequency of full-field output
  int tOut = 1;

  // Heating rate for TDS
	double HR = 900.0/3600.0;

  for (int iStep=1; iStep<nSteps+1; iStep++){

    logger.StepIncrement(iStep);

		// TDS
		model.UpdateTemp(iStep, HR);
		model.WriteTemp(*diff.out, iStep);
		model.CalcElemStiffMatx(Tri3THElemVec, matTVec);
		model.Assemble(Tri3THElemVec);
		linearSolver.UpdateKSP(model.getK());

		// Apply BCs
		model.setBC();
		
		// Solve
		linearSolver.Solve(model.getX(), model.getF());

		// Write total concetration
		model.WriteAvCon(Tri3THElemVec, *diff.out, iStep);

		// Write field output
		if (!(iStep%tOut)){
      // Full field output
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

    return 0;
}