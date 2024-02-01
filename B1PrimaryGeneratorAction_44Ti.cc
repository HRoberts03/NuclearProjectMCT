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
// $Id: B1PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// 44Titanium

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  
  // positioning
  G4double theta = (2*M_PI*G4UniformRand());
  G4double phi = (M_PI*G4UniformRand());
  G4double x = sin(theta)*cos(phi);
  G4double y = sin(theta)*sin(phi);
  G4double z = cos(theta);
  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
  fParticleGun->SetParticleEnergy(0*eV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,15.85*cm)); // this is the position of the source in HPGe detector
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //from data set 1 
    
    
double intensities_sum =  0.34 + 5.6 + 11.1 + 1.20 + 0.61 + 93.0 + 96.8 + 0.092;

double norm_rand = intensities_sum*G4UniformRand();
   
    if(0<norm_rand && norm_rand< 0.34){
	  fParticleGun->SetParticleEnergy(0.4*keV);
    }
	  
    else if(0.34< norm_rand && norm_rand <  0.34+5.6){
	  fParticleGun->SetParticleEnergy(4.086*keV);
   }
    else if(0.34 + 5.6 <norm_rand && norm_rand< 0.34 + 5.6 + 11.1){
	  fParticleGun->SetParticleEnergy(4.091*keV);
     }

    else if(0.34 + 5.6 + 11.1< norm_rand && norm_rand < 0.34 + 5.6 + 11.1 + 1.20){
	  fParticleGun->SetParticleEnergy(4.461*keV);
}
	  
    else if(0.34 + 5.6 + 11.1 + 1.20 < norm_rand && norm_rand < 0.34 + 5.6 + 11.1 + 1.20 + 0.61){
	  fParticleGun->SetParticleEnergy(4.461*keV);
 }
    else if(0.34 + 5.6 + 11.1 + 1.20 + 0.61 < norm_rand && norm_rand < 0.34 + 5.6 + 11.1 + 1.20 + 0.61 + 93.0){
	  fParticleGun->SetParticleEnergy(67.8679*keV);
}
 
    else if(0.34 + 5.6 + 11.1 + 1.20 + 0.61 + 93.0 < norm_rand && norm_rand < 0.34 + 5.6 + 11.1 + 1.20 + 0.61 + 93.0 + 96.8){
	  fParticleGun->SetParticleEnergy(78.3234*keV);
 }
    else if(0.34 + 5.6 + 11.1 + 1.20 + 0.61 + 93.0 + 96.8 < norm_rand && norm_rand < 0.34 + 5.6 + 11.1 + 1.20 + 0.61 + 93.0 + 96.8 + 0.092){
	  fParticleGun->SetParticleEnergy(146.22*keV);
 }
 
  //source position
  G4double theta = (2*M_PI*G4UniformRand());
  G4double phi = (M_PI*G4UniformRand());
  G4double x = sin(theta)*cos(phi);
  G4double y = sin(theta)*sin(phi);
  G4double z = cos(theta);
 
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

