#include <iostream>
#include <vector>

#include "H5IO.h"
#include "Nodes.h"
#include "FiniteElements/Mechanics/Quad4.h"
#include "FiniteElements/PFF/Quad4PFF.h"
#include "Materials/Mechanics/IsoHard.h"
#include "Materials/PFF/PFF.h"
#include "Models/MechModel.h"
#include "Models/PFFModel.h"
#include "Solvers/LinearMech.h"

using namespace std;

int main(int argc, char **argv){

  // Initialize PETSc (always here!)
  PetscInitialize(&argc, &argv, NULL, NULL);

  // Read inputs -----------

  // Model name
  string mechModelName = "Mech_SENT_EPTH";  // Mechanics
  string pffModelName = "PFF_SENT_EPTH";    // Phase-field fracture
  
  // Logger object for handling terminal user interface
  Logger logger(PETSC_COMM_WORLD, mechModelName+".log");
  logger.StartTimer();
  logger.IntroMessage();

  // Initialize I/O files  -------------------------

  const string mechInFileName = mechModelName+"_in.hdf5";
  H5IO mechH5File_in(mechInFileName, logger);

  const string mechOutFileName = mechModelName+"_out.hdf5";
  H5IO mechH5File_out(mechOutFileName, logger);

  const string pffInFileName = pffModelName+"_in.hdf5";
  H5IO pffH5File_in(pffInFileName, logger);

  const string pffOutFileName = pffModelName+"_out.hdf5";
  H5IO pffH5File_out(pffOutFileName, logger);

  // Nodes  -------------------------

  Nodes Nodes;
  Nodes.ReadNodes(mechH5File_in);

  // Material  -------------------------

  vector<BaseMechanics*> mechMatVec;
  mechMatVec.push_back(new IsoHard("2D", mechH5File_in, 1, logger));

  vector<BasePFF*> pffMatVec;
  pffMatVec.push_back(new PFF("2D", pffH5File_in, 1, logger));
 
  // Elements  -------------------------

  vector<BaseElemMech*> mechElemVec;
  mechElemVec.push_back(new Quad4(mechH5File_in, Nodes, 1, "ElastoPlastic", logger));

  vector<BaseElemPFF*> pffElemVec;
  pffElemVec.push_back(new Quad4PFF(pffH5File_in, Nodes, 1, logger));

  // ----------- Initialize the system -----------

  // Mechanics model -----------
  MechModel modelMechanics(mechElemVec, mechH5File_in, logger, "DIRECT", 1);
  
  // Dirichlet BCs
  modelMechanics.InitializeDirichBC(mechH5File_in);

  // Calculate element stiffness matrix
  modelMechanics.CalcElemStiffMatx(mechElemVec, mechMatVec);

  // Assemble global stiffness matrix
  modelMechanics.Assemble(mechElemVec);

  // Write initial conditions
  modelMechanics.CalcNodVals(mechElemVec);
  modelMechanics.WriteOut(mechElemVec, mechH5File_out, to_string(0));

  // PFF model -----------
  PFFModel modelPFF(pffElemVec, pffH5File_in, logger);

  // Sets a constant critical work energy density
  modelPFF.set_const_wc(pffElemVec, pffMatVec);

  // Calculate gPhi_d
  modelPFF.Calc_gPhi_d(pffElemVec);

  // Write initial conditions
  modelPFF.WriteOut(pffH5File_out, to_string(0));

  // Linear solver for PFF
  LinearMech pffLinearSolver(modelPFF.getK());

  // Solver loop -----------

	logger.LoopMessage();

  // Solve the system
  int nSteps = modelMechanics.get_nSteps();

  for (int iStep=1; iStep<nSteps+1; iStep++){

    logger.StepIncrement(iStep);

    // Mechanics -----------

    // Apply BCs
    modelMechanics.setDirichBC();

    // Solve mechanics
    modelMechanics.SolveSNES(mechElemVec, mechMatVec, iStep, pffElemVec);
      
    // Map stresses and strains to nodes
    modelMechanics.CalcNodVals(mechElemVec);

    // Write field output
    modelMechanics.WriteOut(mechElemVec, mechH5File_out, to_string(iStep));

    // PFF -----------

    // Spectral decomposition of elastic strain energy density
    modelPFF.CalcPsiSpectral(pffElemVec, mechElemVec, mechMatVec);

    // *******************

    // Chose a crack driving force by uncommenting ONLY one of the options below:

    // // Brittle crack driving force
    // modelPFF.CalcDrivForcB(pffElemVec);   

    // Elastoplastic crack driving force
    modelPFF.CalcDrivForcEP(pffElemVec, mechElemVec);

    // // Elastoplastic crack driving force with a threshold
    // modelPFF.CalcDrivForcEP_TH(pffElemVec, mechElemVec, 0.25);       

    // *******************
    
    // Build element stiffness matrix
    modelPFF.CalcElemStiffMatx(pffElemVec);
    // Assemble global stiffness matrix
    modelPFF.Assemble(pffElemVec);
    // Calculate RHS
    modelPFF.CalcFH(pffElemVec);

    // Solve PFF (linear)
    pffLinearSolver.Solve(modelPFF.getX(), modelPFF.getF());
    // Calculate gPhi_d
    modelPFF.Calc_gPhi_d(pffElemVec);
    // Write field output
    modelPFF.WriteOut(pffH5File_out, to_string(iStep));

    //     -----------

    logger.FieldOutput(iStep);

  }

  logger.ExitMessage();

  // Deallocate

  // Elements -----------
  for (auto* elem : mechElemVec) { // Elements
    delete elem;
  }

    for (auto* elem : pffElemVec) { // Elements
    delete elem;
  }

  // Materials -----------
  for (auto* mat : mechMatVec) { // Mech material
    delete mat;
  }

  for (auto* mat : pffMatVec) { // Pff material
    delete mat;
  }

  return 0;
}