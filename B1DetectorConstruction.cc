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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "G4VisAttributes.hh"  // i put this

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 50*cm;
  G4double world_sizeZ  = 50*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld = new G4Box("World",0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  

	// Germanium crystal dimensions and material
	G4double crystal_rad = 2.545*cm;
	G4double crystal_hz = 1*cm;  // half (value of be2020)
	G4double cavity_rad = 0.5*cm;  // guess
	G4double cavity_hz = 0.15*cm;  // guess
    G4Material* crystal_mat = nist->FindOrBuildMaterial("G4_Ge");
    
    // Germanium crytal logical volume
    G4Tubs* crystal_shape = new G4Tubs("Crystal", 0*cm, crystal_rad, crystal_hz, 0*deg, 360*deg);
    
    G4Tubs* cavity_shape = new G4Tubs("Cavity", 0*cm, cavity_rad, cavity_hz, 0*deg, 360*deg);
    G4VSolid* subtract = new G4SubtractionSolid("Crystal-Cavity", crystal_shape, cavity_shape, 0, G4ThreeVector(0, 0, -0.84*cm));  //-0.85 cm is a guess from our cavity_hz guess
    G4LogicalVolume* crystal_log = new G4LogicalVolume(subtract, crystal_mat, "Crystal Log");
    
    // Aliminium tube parameters
    G4double tube_outer_rad = 3.35*cm;  // is this right, should be a bit bigger?
    G4double tube_hz = 12.12*cm;  // length is halved, not sure about this
    G4Material* tube_mat = nist->FindOrBuildMaterial("G4_Al"); // not sure about this

    // Create logical volume for the tube
	G4Tubs* tube_shape = new G4Tubs("Tube", 0*cm, tube_outer_rad, tube_hz, 0*deg, 360*deg);
	G4LogicalVolume* tube_log = new G4LogicalVolume(tube_shape, tube_mat, "TubeLog");
	
	// Window parameters
	G4double window_radius = 3.35*cm;  // same as tube
	G4double window_hz = 0.03*cm;  // half
	G4Material* window_mat = nist->FindOrBuildMaterial("G4_C");
	
	//Creates logical volume for the window
	G4Tubs* window_shape = new G4Tubs("Window", 0*mm, window_radius, window_hz, 0*deg, 360*deg);
	G4LogicalVolume* window_log = new G4LogicalVolume(window_shape, window_mat, "WindowLog");
	
	
	// Vacuum parameter
	G4double vacuum_radius = 3.25*cm;  // accounting for 0.1cm Al walls
	G4double vacuum_hz = 12.07*cm;
	G4Material* vacuum_mat = nist->FindOrBuildMaterial("G4_Galactic");
	
	//Creates logical volume for vacuum
	G4Tubs* vacuum_shape = new G4Tubs("Vacuum", 0*mm, vacuum_radius, vacuum_hz, 0*deg, 360*deg);
	G4LogicalVolume* vacuum_log = new G4LogicalVolume(vacuum_shape, vacuum_mat, "VacuumLog");
	
	// GlassA parameter
	G4double GlassA_x = 11*cm;
	G4double GlassA_y = 10*cm;
	G4double GlassA_z = 0.5*cm;
	G4double hole_rad = 4*cm;
	G4Material* glass_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
	
	//Create logical volume for GlassA
	G4Box* glassA_shape = new G4Box("GlassA", GlassA_x, GlassA_y, GlassA_z);
	
	G4Tubs* hole_shape = new G4Tubs("Hole", 0*cm, hole_rad, GlassA_z+1*cm, 0*deg, 360*deg);
    G4VSolid* glassA_new = new G4SubtractionSolid("Glass-Hole", glassA_shape, hole_shape);
    
	G4LogicalVolume* glassA_log = new G4LogicalVolume(glassA_new, glass_mat, "GlassALog");
	
	// GlassB parameter
	G4double GlassB_x = 0.5*cm;
	G4double GlassB_y = 10*cm;
	G4double GlassB_z = 12.9*cm;
	
	G4double small_x = 5.45*cm;
	G4double small_y = 10*cm;
	G4double small_z = 0.5*cm;
	
	//Create logical volume for GlassB
	G4Box* glassB_shape = new G4Box("GlassB", GlassB_x, GlassB_y, GlassB_z);
	
	G4Box* small_shape = new G4Box("Small", small_x, small_y, small_z);
	G4VSolid* glassB_new = new G4UnionSolid("GlassB_fr", glassB_shape, small_shape, 0, G4ThreeVector(-4.95*cm, 0, -13.4*cm));
	G4LogicalVolume* glassB_log = new G4LogicalVolume(glassB_new, glass_mat, "GlassBLog");
	
	// Set visualization attributes for the materials
    G4VisAttributes* crystalVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    crystal_log->SetVisAttributes(crystalVisAtt);
    
    G4VisAttributes* vacuumVisAtt = new G4VisAttributes(G4Colour());
    vacuum_log->SetVisAttributes(vacuumVisAtt);

    G4VisAttributes* tubeVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    tube_log->SetVisAttributes(tubeVisAtt);
    
    G4VisAttributes* windowVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
    window_log->SetVisAttributes(windowVisAtt);
    
    G4VisAttributes* glassVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    glassA_log->SetVisAttributes(glassVisAtt);
    
    G4VisAttributes* glassBVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    glassB_log->SetVisAttributes(glassBVisAtt);
    

	// Place
	
	G4RotationMatrix* GlassC_rot = new G4RotationMatrix();
	GlassC_rot->rotateZ(180*deg);
	
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), tube_log, "Al Tube", logicWorld, false, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0.02*cm), vacuum_log, "Vacuum", tube_log, false, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 10.52*cm), crystal_log, "Ge Crystal", vacuum_log, false, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 12.15*cm), window_log, "Window", logicWorld, false, checkOverlaps);  // z is tube length - window thickness
    new G4PVPlacement(0, G4ThreeVector(0, 0, 14.6*cm), glassA_log, "GlassA", logicWorld, false, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(-10.5*cm, 0, 1.2*cm), glassB_log, "GlassB", logicWorld, false, checkOverlaps);
    new G4PVPlacement(GlassC_rot, G4ThreeVector(10.5*cm, 0, 1.2*cm), glassB_log, "GlassC", logicWorld, false, checkOverlaps);


 //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
