#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "FiniteElements/Mechanics/Quad4.h"
#include "Materials/Mechanics/IsoHard.h"
#include "Models/MechModel.h"

using namespace std;

int main(int argc, char **argv){

  // use ./CrackTipQuad4 -snes_monitor to view covergence behavior

  // Initialize PETSc (always here!)
  PetscInitialize(&argc, &argv, NULL, NULL);

  // Read inputs -----------

  // Model name, same as `Simul`
  string modelName = "CrackTip";
  
  // Logger object for handling terminal user interface
  Logger logger(PETSC_COMM_WORLD, modelName+".log");
  logger.StartTimer();
  logger.IntroMessage();

  // Initialize I/O files
  const string infileName = modelName+"_in.hdf5";
  H5IO H5File_in(infileName, logger);

  const string outfileName = modelName+"_out.hdf5";
  H5IO H5File_out(outfileName, logger);

  // Material
  vector<BaseMechanics*> matVec;
  matVec.push_back(new IsoHard("2D", H5File_in, 1, logger));

  // Nodes
  Nodes Nodes;
  Nodes.ReadNodes(H5File_in);
 
  // Elements
  vector<BaseElemMech*> quad4ElemVec;
  quad4ElemVec.push_back(new Quad4(H5File_in, Nodes, 1, "ElastoPlastic", logger));

  // Initialize the sytem -----------

  // Model
  MechModel model(quad4ElemVec, H5File_in, logger, 1);
  
  // Dirichlet BCs
  model.InitializeDirichBC(H5File_in);

  // Calculate element stiffness matrix
  model.CalcElemStiffMatx(quad4ElemVec, matVec);

  // Assemble global stiffness matrix
  model.Assemble(quad4ElemVec);

  // Write initial conditions
  model.CalcNodVals(quad4ElemVec);
  model.WriteOut(quad4ElemVec, H5File_out, to_string(0));

  // Solver loop -----------

	logger.LoopMessage();

  // Solve the system
  int nSteps = model.get_nSteps();

  for (int iStep=1; iStep<nSteps+1; iStep++){

    logger.StepIncrement(iStep);

    // Apply BCs
    model.setDirichBC();

    // Non-linear solver
    model.SolveSNES(quad4ElemVec, matVec, iStep);
      
    // Map stresses and strains to nodes
    model.CalcNodVals(quad4ElemVec);

    // Write field output
    model.WriteOut(quad4ElemVec, H5File_out, to_string(iStep));

    logger.FieldOutput(iStep);

  }

  logger.ExitMessage();

  // Deallocate
  for (auto* elem : quad4ElemVec) { // Elements
    delete elem;
  }
    for (auto* mat : matVec) { // Material
    delete mat;
  }

  return 0;
}