
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
#include <G4GeneralParticleSource.hh>

using namespace std;

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
{
  fGPS = new G4GeneralParticleSource();
  G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma"); // gets gamma from source list
  fGPS->SetParticleDefinition(gamma); // sets particle to gamma
  fGammaGPS = fGPS->GetCurrentSource();  // shortcut to access source
  fGammaGPS->GetPosDist()->SetPosDisType("Surface"); // sets the type of source
  fGammaGPS->GetPosDist()->SetPosDisType("Circle"); // sets the shape of source
  fGammaGPS->GetPosDist()->SetRadius(1.04*mm); 
  
  fGammaGPS->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.)); // sets particle position
  fGammaGPS->GetAngDist()->SetAngDistType("iso"); // isotropic source
  fGammaGPS->GetEneDist()->SetMonoEnergy(0*keV); // defualt energy value
}

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fGPS;
}

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
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
    
   
   
     // Calculate total intensity
    std::vector<double> cumulative_sums;
    double sum_value = 0;
    // appends the cumalitive sum for the specific emmision
    for (size_t i = 0; i < Co60_size; ++i) {
        double value = Co60_intensities[i];
        sum_value += value;
        cumulative_sums.push_back(sum_value);
    }
    // generates normalised random number for emission selection
    double cumulative_sum = cumulative_sums.back();
    double norm_rand = cumulative_sum * G4UniformRand();
    
    // source generation
   
    for (size_t i = 0; i < Co60_size; ++i) {
        if (i == 0 && 0 < norm_rand && norm_rand < Co60_intensities[0]) {
            // Set the energy and exit the function
            fGammaGPS->GetEneDist()->SetMonoEnergy(Co60_energies[0]*keV);

        } 
        
        else if (i > 0 && cumulative_sums[i - 1] < norm_rand && norm_rand < cumulative_sums[i]) {
            // Set the energy and exit the function
             fGammaGPS->GetEneDist()->SetMonoEnergy(Co60_energies[i]*keV);
        }
    }
  // generate the selected event
  fGPS->GeneratePrimaryVertex(anEvent);
}
