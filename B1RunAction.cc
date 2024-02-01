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
// $Id: B1RunAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
//#include "G4ParameterManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <array> 
#include <cmath>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector <G4double> EnergyArrayRaw;
  

B1RunAction::B1RunAction()
: G4UserRunAction()
  //fEdep("Edep", 0.),
  //fEdep2("Edep2", 0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register parameter to the parameter manager
  //G4ParameterManager* parameterManager = G4ParameterManager::Instance();
  //parameterManager->RegisterParameter(fEdep);
  //parameterManager->RegisterParameter(fEdep2); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

 

  // reset parameters to their initial values
  //G4ParameterManager* parameterManager = G4ParameterManager::Instance();
  //parameterManager->Reset();
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

// smears the data as a detector would
G4double m = 10 * G4UniformRand();
G4double c = 5 * G4UniformRand();
std::vector <G4double> energiesLinearGain;
for(unsigned int i = 0; i<EnergyArrayRaw.size(); i++){
  energiesLinearGain.push_back(m*EnergyArrayRaw[i] + c);
};


//finds max and min value of the smeared energy
G4double maxE = energiesLinearGain[0];
for(unsigned int i = 0; i<energiesLinearGain.size(); i++){
  if(energiesLinearGain[i] > maxE){
    maxE = energiesLinearGain[i];
  };
};
G4double minE = energiesLinearGain[0];
for(unsigned int i = 0; i<energiesLinearGain.size(); i++){
  if(energiesLinearGain[i] < minE){
    minE = energiesLinearGain[i];
  };
};

//create the bins locations
G4double rangeE = maxE - minE;
int rangeEint = floor(rangeE)/2;
std::vector <int> Bins;
for(unsigned int i=0; i<rangeEint; i++){
  Bins.push_back(i*2);
};

//turns all the energies to their bin values
for(int i : Bins){
  int count =0;
for(int j : energiesLinearGain){
  if(abs(i-j)<1){
    energiesLinearGain[count] = i;
    ++ count;
  };};};


//create a vector of unique energies. 
std::vector <G4double> UniqueEnergies;
UniqueEnergies.push_back(energiesLinearGain[0]);
for(G4double i : energiesLinearGain){
  bool Exists = false;
  for(G4double j :UniqueEnergies){
    if(i==j){
      Exists = true;
    };
  };
  if(Exists == false){
    UniqueEnergies.push_back(i);
  };
};

//gets the count rate of each energy.
std::vector <G4double> CountRate;
for(G4double i : UniqueEnergies){
  int count3 = 0;
  for(G4double j : energiesLinearGain){
    if(i==j){
      ++ count3;
    };
  };
  CountRate.push_back(count3);
};
int zeroIndex = 0;
for(unsigned int i = 0; i < UniqueEnergies.size(); i++){
  if(UniqueEnergies[i]==0){
    zeroIndex = i;
    break;
  };
};
UniqueEnergies.erase(UniqueEnergies.begin() + zeroIndex);
CountRate.erase(CountRate.begin() + zeroIndex);

std::ofstream myFile;
    myFile.open("output1.csv");
  for(unsigned int i=0; i < UniqueEnergies.size(); i++){
    
    myFile << UniqueEnergies[i] <<" " << CountRate[i]<< "\n";
   
};
  myFile.close();
  
  // Merge parameters 
  //G4ParameterManager* parameterManager = G4ParameterManager::Instance();
  //parameterManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  //G4double edep  = fEdep.GetValue();
  //G4double edep2 = fEdep2.GetValue();
  
  //G4double rms = edep2 - edep*edep/nofEvents;
  //if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  
       
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " particles"
     << G4endl
    << "------------------------------------------------------------"
     << G4endl
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
//  fEdep  += edepG4doubl;
//  fEdep2 += edep*edepG4doubl;
};

void B1RunAction::AddEnergy(G4double energy)
{
  EnergyArrayRaw.push_back(energy);
};
  






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

