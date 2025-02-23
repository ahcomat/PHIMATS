#include <iostream>
#include <vector>

#include "petscsys.h"
#include "H5IO.h"
#include "Materials/Trapping/TrapPhase.h"
#include "Nodes.h"
#include "FiniteElements/Trapping/Tri3TH.h"
#include "Models/TrappingModel.h"
#include "Solvers/LinearTransport.h"
#include "Logger.h"

using namespace std;

int main(int argc, char **argv){

  	// Initialize PETSc 
	PetscInitialize(&argc, &argv, NULL, NULL);

  // Read inputs -----------

	// Model name, same as `Simul`
	string modelName = "TDSHT990";

  	// Logger object for handling terminal user interface
	Logger logger(PETSC_COMM_WORLD, modelName+".log");
	logger.StartTimer();
	logger.IntroMessage();

  // Initialize I/O hdf5 files
	const string infileName = modelName+"_in.hdf5";
	H5IO H5File_in(infileName, logger);

	const string outfileName = modelName+"_out.hdf5";
	H5IO H5File_out(outfileName, logger);

  // Initial conditions for TDS simulation
  const string infileName2 = "../Charging/ChargingHT990_out.hdf5";
	H5IO H5File_in2(infileName2, logger);

  // Material vector
  vector<BaseTrapping*> matTVec;
  matTVec.push_back(new TrapPhase("2D", H5File_in, 1, logger));

  // Nodes
  Nodes Nodes;
  Nodes.ReadNodes(H5File_in);
 
  // Elements
  vector<BaseElemTrap*> Tri3THElemVec;
  Tri3THElemVec.push_back(new Tri3TH(H5File_in, Nodes, 1, logger));

  // Initialize the sytem -----------

  // Model
  TrappingModel model(Tri3THElemVec, H5File_in, logger);
 
  // Calculate the stiffness matrix
  model.CalcElemStiffMatx(Tri3THElemVec, matTVec);

  // Initialize boundary conditions (must be called before `Assemble`!)
	model.InitializeBC(H5File_in);

  // Assemble global stiffness matrix
  model.Assemble(Tri3THElemVec);

  // Solver: "GMRES" solver is more efficient for large system TDS simulations because
  // it handles the frequent updates to the stiffness matrix more effectively
  LinearTransport linearSolver(model.getK(), logger, "GMRES");

  // Read initial data (Last step in the charing simulation)
	model.ReadInitialCon(H5File_in2, 20);

	// Write initial conditions
  model.WriteTemp(H5File_out, 0);
	model.WriteAvCon(Tri3THElemVec, H5File_out, 0);
	model.WriteOut( H5File_out, to_string(0));
  
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
		model.WriteTemp(H5File_out, iStep);
		model.CalcElemStiffMatx(Tri3THElemVec, matTVec);
		model.Assemble(Tri3THElemVec);
		linearSolver.UpdateKSP(model.getK());

		// Apply BCs
		model.setBC();
		
		// Solve
		linearSolver.Solve(model.getX(), model.getF());

		// Write total concetration
		model.WriteAvCon(Tri3THElemVec, H5File_out, iStep);

		// Write field output
		if (!(iStep%tOut)){
      // Full field output
      model.WriteOut( H5File_out, to_string(iStep));
      logger.FieldOutput(iStep);
		}
	}

  // Deallocate
  for (auto* elem : Tri3THElemVec) { // Elements
    delete elem;
  }

    for (auto* mat : matTVec) { // Material
    delete mat;
  }

  return 0;
}