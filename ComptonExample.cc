//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: exampleB1.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "B1DetectorConstruction.hh"
#include "B1DetectorConstruction_HPGe.hh"
#include "B1DetectorConstruction_GAGG.hh"
#include "B1DetectorConstruction_Shipping_Container.hh"
#include "B1DetectorConstruction_HPGe_Shipping_Container.hh"
#include "B1DetectorConstruction_GAGG_Shipping_Container.hh"
#include "B1ActionInitialization.hh"
#include <vector>

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"
#include "FTFP_BERT.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  std::vector<G4String> macros;
  G4bool interactive = false;
  G4bool shipping = false;
  G4int detectorNum = 0;
  G4int particleNum = 0;

  if ( argc == 1 )
  {
    interactive = true;
  }
  else
  {
    for (int i = 1; i < argc; i++)
    {
      G4String arg = argv[i];
      if (arg == "-i" || arg == "--interactive")
      {
        interactive = true;
      }
      else if (arg == "-d" || arg == "--detector")
      {
        detectorNum = std::stoi(argv[i+1]);
        i++;
      }
      else if (arg == "-p" || arg == "--particle")
      {
        particleNum = std::stoi(argv[i+1]);
        i++;
      }
      else if (arg == "-s" || arg == "--shipping")
      {
          shipping = true;
      }
      else
      {
        macros.push_back(arg);
      }
    }
  }
  // Choose the Random engine I HAVE CHANGED THIS - NOTE FROM JACK TRYING TO FIX LACK OF RANDOMNESS
  // G4Random::setTheEngine(new CLHEP::RanecuEngine);  < -----TONIES OLD CODE POSSIBLY PREVENTING RANDOMNESS
  //set the random engine
  
  
  /// New random number generator to generate independent trials 

CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
G4long seed = time(NULL);    //use system time as seed
CLHEP::HepRandom::setTheSeed(seed); // set the random seed
  
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // Detector construction
  //if (hpGe) runManager->SetUserInitialization(new B1DetectorConstruction_HPGe());
  //else runManager->SetUserInitialization(new B1DetectorConstruction());

  if (detectorNum==1)
  {
      runManager->SetUserInitialization(new B1DetectorConstruction());
  }
  
  else if (detectorNum==2)
  {
   runManager->SetUserInitialization(new B1DetectorConstruction_HPGe());   
  }
  else if (detectorNum==3)
  {
   runManager->SetUserInitialization(new B1DetectorConstruction_GAGG());   
  }
  else
  {
      runManager->SetUserInitialization(new B1DetectorConstruction());
  }
  
  if (shipping){
        if (detectorNum==1)
  {
      runManager->SetUserInitialization(new B1DetectorConstruction_Shipping_Container());
  }
  
  else if (detectorNum==2)
  {
   runManager->SetUserInitialization(new B1DetectorConstruction_HPGe_Shipping_Container());   
  }
  else if (detectorNum==3)
  {
   runManager->SetUserInitialization(new B1DetectorConstruction_GAGG_Shipping_Container());   
  }
  else
  {
      runManager->SetUserInitialization(new B1DetectorConstruction());
  }
  }
  //runManager->SetUserInitialization(new B1DetectorConstruction());
  
  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4RadioactiveDecayPhysics);
  physicsList->SetVerboseLevel(0);
  runManager->SetUserInitialization(physicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new B1ActionInitialization(particleNum, detectorNum));
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  G4UIExecutive* ui = nullptr;
  if (interactive)
    {
      G4cout << "Creating interactive UI session ...";
      ui = new G4UIExecutive(argc, argv);
  }
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //

  for (auto macro : macros)
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }

  if (interactive)
  {
    if (ui->IsGUI())
    {
      UImanager->ApplyCommand("/control/execute init_vis.mac");
    }
    else
    {
      UImanager->ApplyCommand("/run/initialize");
    }
    ui->SessionStart();
    delete ui;
  }
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
