
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


B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{

  // default particle kinematic
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  
  // positioning
  G4double theta = (2*M_PI*G4UniformRand());
  G4double phi = (M_PI*G4UniformRand());
  G4double x = sin(theta)*cos(phi);
  G4double y = sin(theta)*sin(phi);
  G4double z = cos(theta);
  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
  fParticleGun->SetParticleEnergy(0*eV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
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
    
    
double intensities_sum =  7.5*std::pow(10, -4) + 0.0078 + 0.25 + 2.0601 + 0.0075 +  0.0076 + 99.85
+ 99.9826 + 0.00120 + 2*10*-6;

double norm_rand = intensities_sum*G4UniformRand();
   
    if(0<norm_rand && norm_rand< 7.5*std::pow(10, -4)){
	  fParticleGun->SetParticleEnergy(2158.57*keV);
    }
	  
   if(7.5*std::pow(10, -4)< norm_rand && norm_rand <  0.0078+ 7.5*std::pow(10, -4)){
	  fParticleGun->SetParticleEnergy(826.10*keV);
   }
     if(0.0078+ 7.5*std::pow(10, -4)<norm_rand && norm_rand<0.0078+ 7.5*std::pow(10, -4) + 0.25){
	  fParticleGun->SetParticleEnergy(1332.492*keV);
     }

  if(0.0078+ 7.5*std::pow(10, -4) + 0.25< norm_rand && norm_rand < 2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25){
	  fParticleGun->SetParticleEnergy(58.603*keV);
}
	  

// from data set 3
 if(2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25< norm_rand && norm_rand <   2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075){
	  fParticleGun->SetParticleEnergy(347.14*keV);
 }
	  
	  
if(2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 < norm_rand && norm_rand <   2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076){
	  fParticleGun->SetParticleEnergy(826.10*keV);
}
 
 if(2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076< norm_rand && norm_rand <2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076 +  99.85){
	  fParticleGun->SetParticleEnergy(1173.228*keV);
 }
 if(2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076 +  99.85< norm_rand && norm_rand < 2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076 +  99.85 +  99.9826){
	  fParticleGun->SetParticleEnergy(1332.492*keV);
 }
	  
if(2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076 +  99.85 +  99.9826< norm_rand && norm_rand < 2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076 +  99.85 +  99.9826 + 0.00120){
	  fParticleGun->SetParticleEnergy(2158.57*keV);
}

 
 if(2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076 +  99.85 +  99.9826 + 0.00120< norm_rand && norm_rand <   2.0601 + 0.0078+ 7.5*std::pow(10, -4) + 0.25 +0.0075 +0.0076 +  99.85 +  99.9826 + 0.00120 + 2*10*-6){
	  fParticleGun->SetParticleEnergy(2505.692*keV);
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

