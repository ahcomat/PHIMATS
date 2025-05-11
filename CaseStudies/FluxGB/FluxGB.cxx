#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "Materials/Trapping/TrapGB.h"
#include "FiniteElements/Trapping/Tri3TH.h"
#include "Models/TrappingModel.h"
#include "Solvers/LinearTransport.h"
#include "Logger.h"

using namespace std;

int main(int argc, char **argv){
  
	// Initialize PETSc (always here!)
	PetscInitialize(&argc, &argv, NULL, NULL);

	// Read inputs -----------

	// Model name, same as `Simul`
	string modelName = "FluxGB";

	// Logger object for handling terminal user interface
	Logger logger(PETSC_COMM_WORLD, modelName+".log");
	// Start timer
	logger.StartTimer();
	// Print intro message
	logger.IntroMessage();

	// Initialize I/O files
	const string infileName = modelName+"_in.hdf5";
	H5IO H5File_in(infileName, logger);

	const string outfileName = modelName+"_out.hdf5";
	H5IO H5File_out(outfileName, logger);

	// Materials vector
	vector<BaseTrapping*> matTVec;
	matTVec.push_back(new TrapGB("2D", H5File_in, 1, logger));

	// Nodes
	Nodes Nodes;
	Nodes.ReadNodes(H5File_in);

	// Elements vector
	vector<BaseElemTrap*> Tri3THElemVec;
	Tri3THElemVec.push_back(new Tri3TH(H5File_in, Nodes, 1, logger));

	// Initialize the system -----------

	// Model
	TrappingModel model(Tri3THElemVec, H5File_in, logger);

	// Calculate the stiffness matrix of all elements
	model.CalcElemStiffMatx(Tri3THElemVec, matTVec);
	
	// Initialize boundary conditions
	model.InitializeBC(H5File_in);

	// Assemble global stiffness matrix
	model.Assemble(Tri3THElemVec);

	// Solver: "DIRECT" solver much faster for charging/permeation simulations
	LinearTransport linearSolver(model.getK(), logger, "DIRECT");

	// Write initial conditions
	model.WriteAvCon(Tri3THElemVec, H5File_out, 0);
	model.WriteOut( H5File_out, to_string(0));
	model.CalcFlux(Tri3THElemVec, matTVec);
	model.WriteExitFlux(H5File_out, 0);
	model.WriteFlux(H5File_out, to_string(0));

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
		model.WriteAvCon(Tri3THElemVec, H5File_out, iStep);

		// Write exit side flux (every step)
		model.WriteExitFlux(H5File_out, iStep);

		// Write field output
		if (!(iStep%tOut)){
			
			// Flux
			model.WriteFlux(H5File_out, to_string(iStep));
			// Concentration
			model.WriteOut( H5File_out, to_string(iStep));

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