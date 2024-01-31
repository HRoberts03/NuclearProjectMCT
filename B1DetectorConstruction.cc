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
#include <G4VisAttributes.hh>


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
  
  //Colours
  G4VisAttributes* green = new G4VisAttributes(G4Colour::Green());
  green->SetVisibility(true);
  green->SetForceSolid(true);
  G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());
  red->SetVisibility(true);
  red->SetForceSolid(true);
  G4VisAttributes* blue = new G4VisAttributes(G4Colour::Blue());
  blue->SetVisibility(true);
  blue->SetForceSolid(true);
  G4VisAttributes* yellow = new G4VisAttributes(G4Colour::Yellow());
  yellow->SetVisibility(true);
  yellow->SetForceSolid(true);
  
  G4double world_sizeXY = 1000*cm;
  G4double world_sizeZ  = 1000*cm;
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
  

    // Table Volume
    G4double table_sizeX = 100*cm;
    G4double table_sizeY = 10*cm;
    G4double table_sizeZ  = 100*cm;
    G4Material* table_mat = nist->FindOrBuildMaterial("G4_CELLULOSE_BUTYRATE");
    G4Box* tableWalls = new G4Box("table",0.5*table_sizeX, 0.5*table_sizeY, 0.5*table_sizeZ);     //its size
    G4LogicalVolume* table_log = new G4LogicalVolume(tableWalls, table_mat, "table"); //its name
    
                                   
    //walls
      G4double wall_sizeXY = 500*cm;
      G4double wall_sizeZ  = 500*cm;
      G4Material* wall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
  
     G4Box* solidWalls = new G4Box("Walls",0.5*wall_sizeXY, 0.5*wall_sizeXY, 0.5*wall_sizeZ);     //its size
      
    G4LogicalVolume* logicwalls = new G4LogicalVolume(solidWalls, wall_mat, "wall");    //its name
    
     
  G4ThreeVector wall_pos = G4ThreeVector(0,0,0);
  G4VPhysicalVolume* wall_physical = new G4PVPlacement(0, wall_pos, "wall", logicwalls, physWorld, false, 0, checkOverlaps);	

    // Air
  //
  G4double air_sizeXY = 450*cm;
  G4double air_sizeZ  = 450*cm;
  G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidAir = new G4Box("air",0.5*air_sizeXY, 0.5*air_sizeXY, 0.5*air_sizeZ);     //its size
      
  G4LogicalVolume* airWorld = new G4LogicalVolume(solidAir, air_mat, "inner"); //its name
  
  G4ThreeVector air_pos = G4ThreeVector(0,0,0);
  G4VPhysicalVolume* air_physical = new G4PVPlacement(0, air_pos, "air", airWorld, wall_physical, false, 0, checkOverlaps);	

	//face construction 
	G4double face_r = 5.71*cm;
	G4double back_hz =  2*mm;
	G4double face_hz =  0.5*mm;
	G4Tubs* face_shape = new G4Tubs("face", 0*mm, face_r, face_hz, 0*deg, 360*deg);
	G4LogicalVolume* face_log = new G4LogicalVolume(face_shape, nist->FindOrBuildMaterial("G4_Al"), "face");
	face_log->SetVisAttributes(red);


	// scintillator construction
	G4double scint_r = 5.51*cm;
	G4double scint_hz = 5.40*cm;
	G4Tubs* scintillator_shape = new G4Tubs("scintillator",0*mm,scint_r,scint_hz,0*deg,360*deg);
	G4LogicalVolume* scint_log = new G4LogicalVolume(scintillator_shape, nist->FindOrBuildMaterial("G4_SODIUM_IODIDE"), "scintillator");
    scint_log->SetVisAttributes(blue);
    // Shield Construction
    G4double shield_r = 5.71*cm;
	G4double shield_hz = 5.40*cm + face_hz + back_hz;
	G4Tubs* shield_shape = new G4Tubs("shield",0,shield_r,shield_hz,0*deg,360*deg);
	G4LogicalVolume* shield_log = new G4LogicalVolume(shield_shape, nist->FindOrBuildMaterial("G4_Al"), "shield");
	shield_log->SetVisAttributes(yellow);
	
	//Detector Back Construction
	G4double back_r = 5.71*cm;
	G4Tubs* back_shape = new G4Tubs("back", 0*mm, back_r, back_hz, 0*deg, 360*deg);
	G4LogicalVolume* back_log = new G4LogicalVolume(back_shape, nist->FindOrBuildMaterial("G4_Al"), "face");
	back_log->SetVisAttributes(green);


	
	// Apparatus Shifts
    G4double face_shift = 0*cm;
    G4double main_detector_shift = face_shift + 2*face_hz;
    G4double table_shift = shield_r + table_sizeY*0.5;
    G4double hz_shift = 20*cm;


	
/*
	G4Box* shape = new G4Box("shape", 1*mm,1*mm,1*mm);
	G4LogicalVolume* LV = new G4LogicalVolume(shape,world_mat,"LV");
		G4PVPlacement(0,G4ThreeVector(0,0*cm,0), LV,"PV", NaI_log, false,0,checkOverlaps);
*/

    

	// shield Placement
	G4ThreeVector shield_pos = G4ThreeVector(0,0,hz_shift + shield_hz  + face_hz);
	G4RotationMatrix* shield_rot = new G4RotationMatrix();
	shield_rot -> rotateX(0*deg);
	G4VPhysicalVolume* shield = new G4PVPlacement(0, shield_pos, "shield", shield_log, air_physical, false, 0, checkOverlaps);
    
    // Scintillator Placement
	G4ThreeVector scintillator_pos = G4ThreeVector(0,0,0);
	G4RotationMatrix* scintillator_rot = new G4RotationMatrix();
	scintillator_rot -> rotateX(0*deg);
	G4VPhysicalVolume* scintillator = new G4PVPlacement(0, scintillator_pos, "scintillator", scint_log, shield, false, 0, checkOverlaps);
		
	//Table placement
    G4ThreeVector table_pos = G4ThreeVector(0,-table_shift,0);
	G4RotationMatrix* table_rot = new G4RotationMatrix();
	table_rot->rotateX(0*deg);
    G4VPhysicalVolume* table = new G4PVPlacement(0, table_pos, "table", table_log, air_physical, false, 0, checkOverlaps);
    
    // back placement
    G4ThreeVector back_pos = G4ThreeVector(0,0,shield_hz);
	G4RotationMatrix* back_rot = new G4RotationMatrix();
	back_rot->rotateX(0*deg);
    G4VPhysicalVolume* back_physical = new G4PVPlacement(0, back_pos, "back", back_log, shield, false, 0, checkOverlaps);	
    
     //face placement
	G4ThreeVector face_pos = G4ThreeVector(0,0,-1*shield_hz);
	G4RotationMatrix* face_rot = new G4RotationMatrix();
	face_rot->rotateX(0*deg);
    G4VPhysicalVolume* face = new G4PVPlacement(0, face_pos, "face", face_log, shield, false, 0, checkOverlaps);	



 //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
