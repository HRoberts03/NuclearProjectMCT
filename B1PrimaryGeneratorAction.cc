
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
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{

  // default particle kinematic
  G4int n_particle = 1000;
  fParticleGun  = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  
  // positioning
  G4double cosTheta = 2*G4UniformRand() - 1., phi = CLHEP::twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double ux = sinTheta*std::cos(phi),
           uy = sinTheta*std::sin(phi),
           uz = cosTheta;
  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  fParticleGun->SetParticleEnergy(0*eV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));  

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
    
    // Source Definitions

     //Co60 
    std::vector<double> Co60_intensities = {7.5 * std::pow(10, -4), 0.0078, 0.25, 2.0601, 0.0075,
                                             0.0076, 99.85, 99.9826, 0.00120, 2 * std::pow(10, -6)};
    std::vector<double> Co60_energies = {2158.57, 826.10, 1332.492, 58.603, 347.14,
                                             826.10, 826.10, 1332.492, 2158.57, 2505.692};
    size_t Co60_size = Co60_intensities.size();
    
    //Cs137
    std::vector<double> Cs137_intensities = {0.91,1.99,3.64,0.348,0.672,0.213,5.8*std::pow(10,-4),85.10};
    
    std::vector<double> Cs137_energies = {4.47,31.817,32.194,36.304,36.378,37.255,283.5,661.657};
    
    size_t Cs137_size = Cs137_intensities.size();
    
    //Ba133
    
    std::vector<double> Ba133_intensities = { 16.4,1.42,15.1,27.6,2.64,5.10,1.61,17.69,1.6*std::pow(10,-4),15.7,33.9,62.2,
    5.88,11.4,3.51,2.14,2.65,32.9,0.638,0.453,7.16,18.34,62.05,8.94};
    
    std::vector<double> Ba133_energies = {4.47,12.327,31.817,32.194,36.304,36.378,37.255,275.925,288,4.29,30.625,30.973,34.92,34.987,35.818,53.1622,79.6142,80.9979,160.6120,223.2368,276.3989,302.8508,356.0129,383.8485};
    
    size_t Ba133_size = Ba133_intensities.size();
    
    
    
    
    //generates source from vectors
    
    std::vector<double> cumulative_sums;
    double sum_value = 0;

    for (size_t i = 0; i < Co60_size; ++i) {
        double value = Co60_intensities[i];
        sum_value += value;
        cumulative_sums.push_back(sum_value);
    }

    double cumulative_sum = cumulative_sums.back();
    double norm_rand = cumulative_sum * G4UniformRand();

    for (size_t i = 0; i < Co60_size; ++i) {
        if (i == 0 && 0 < norm_rand && norm_rand < Co60_intensities[0]) {
            // Set the energy and exit the function
            fParticleGun->SetParticleEnergy(Co60_energies[0] * keV);
        } 
        
        else if (i > 0 && cumulative_sums[i - 1] < norm_rand && norm_rand < cumulative_sums[i]) {
            // Set the energy and exit the function
            fParticleGun->SetParticleEnergy(Co60_energies[i] * keV);
        }
    }


 
  // Isotropic Source
  
  //Code to cover total solid angle
  
 G4double cosTheta = 2*G4UniformRand() - 1., phi = CLHEP::twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double ux = sinTheta*std::cos(phi),
           uy = sinTheta*std::sin(phi),
           uz = cosTheta;

// Apply isotropic properties

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
  
