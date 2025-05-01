#include <iostream>
#include <vector>

#include "Eigen/Dense"
#include "H5IO.h"
#include "Nodes.h"

#include "FiniteElements/Trapping/Quad4TH.h"
#include "Materials/Trapping/MechTrap.h"
#include "Models/TrappingModel.h"
#include "Solvers/LinearTransport.h"

using namespace std;

int main(int argc, char **argv){

  // Initialize PETSc (always here!)
  PetscInitialize(&argc, &argv, NULL, NULL);

  // Initialize I/O ----------

  string modelName = "CrackTip_Diffusion_130s";

  // Logger object for handling terminal user interface
  Logger logger(PETSC_COMM_WORLD, modelName+".log");
  logger.StartTimer();
  logger.IntroMessage();

  // Input file
  const string inFileName = modelName+"_in.hdf5";
  H5IO H5File_in(inFileName, logger);

  // Output file
  const string outFileName = modelName+"_out.hdf5";
  H5IO H5File_out(outFileName, logger);

  // Mechanical output file, which is input to the diffusion simulation
  const string stressFileName = "../Loading/CrackTip_out.hdf5";
  H5IO H5File_stress(stressFileName, logger);

  // Material: hydrogen mechanics interaction
  vector<BaseTrapping*> trapMatVec;
  trapMatVec.push_back(new MechTrap("2D", H5File_in, 1, logger));

  // Nodes ----------
  Nodes Nodes;
  Nodes.ReadNodes(H5File_in);
 
  // Finite elements: Trapping element type
  vector<BaseElemTrap*> TrapElemVec;
  TrapElemVec.push_back(new Quad4TH(H5File_in, Nodes, 1, logger));

  // Read initial nodal stresses from the input stress file.
  TrapElemVec[0]->ReadNodalStress(H5File_stress, 0);

  // Initialize the system -----------

  // Models: Trapping model
  TrappingModel model(TrapElemVec, H5File_in, logger);

  // Sets initial uniform concentration
  model.setUniformCon(0.003453);

  // Calculate the stiffness matrix of all elements
  model.CalcElemStiffMatx(TrapElemVec, trapMatVec);

  // Dirichlet BCs (has to use this function even if only zero flux BCs are imposed)
  model.InitializeBC(H5File_in);

  // Assemble global stiffness matrix
  model.Assemble(TrapElemVec);

  // Output Initial values
  model.WriteAvCon(TrapElemVec, H5File_out, 0);
	model.WriteOut(H5File_out, to_string(0));

  // Solver: Linear with direct solver
  LinearTransport linearSolver(model.getK(), logger, "DIRECT");

  logger.LoopMessage();

  // Solver loop -----------

  // Number of steps
  int nSteps = model.get_nSteps();
  // Frequency of full-field output
  int tOut = 1;

  // Solving diffusion while increasing the stress
  logger.log("  >>> Stress diffusion loop <<<\n", "INFO");

  for (int iStep=1; iStep<nSteps+1; iStep++){

    logger.StepIncrement(iStep);

    // Read stress for the step
    TrapElemVec[0]->ReadNodalStress(H5File_stress, iStep);
    
    // The stiffness matrix and KSP have to be updated for the new stress
    model.CalcElemStiffMatx(TrapElemVec, trapMatVec);
		model.Assemble(TrapElemVec);
		linearSolver.UpdateKSP(model.getK());

    // Sets BC and updates the RHS. Required even for zero flux BCs
    model.setBC();

    //Solve
    linearSolver.Solve(model.getX(), model.getF());

    // Write field output (in this case, each time increment)
    if (!(iStep%tOut)){
      
      model.WriteOut(H5File_out, to_string(iStep));
      model.WriteAvCon(TrapElemVec, H5File_out, iStep);

      logger.FieldOutput(iStep);
    }

  }

  // Solving diffusion while holding the stress
  logger.log("  >>> Hold loop <<<\n", "INFO");

  // Set new time increment
  model.Update_dt(TrapElemVec, 130);

  // NOTE: Since stress is being held constant, no need to update the stiffness matrix during the loop
  model.CalcElemStiffMatx(TrapElemVec, trapMatVec);
  model.Assemble(TrapElemVec);
  linearSolver.UpdateKSP(model.getK());

  for (int iStep=nSteps+1; iStep<nSteps+501; iStep++){

    logger.StepIncrement(iStep);

    // Update boundary conditions
    model.setBC();

    //Solve
    linearSolver.Solve(model.getX(), model.getF());

    // Write field output
    if (!(iStep%tOut)){
      
      model.WriteOut(H5File_out, to_string(iStep));
      model.WriteAvCon(TrapElemVec, H5File_out, iStep);

      logger.FieldOutput(iStep);
    }

  }

  // Print simulation exit message upon successful completion. 
  logger.ExitMessage();

  // Deallocate
  for (auto* elem : TrapElemVec) { // Elements
    delete elem;
  }
    for (auto* mat : trapMatVec) { // Material
    delete mat;
  }

  return 0;
}