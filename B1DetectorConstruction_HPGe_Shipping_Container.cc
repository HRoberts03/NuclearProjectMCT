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

#include "B1DetectorConstruction_HPGe_Shipping_Container.hh"
#include "G4VisAttributes.hh"
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

B1DetectorConstruction_HPGe_Shipping_Container::B1DetectorConstruction_HPGe_Shipping_Container()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction_HPGe_Shipping_Container::~B1DetectorConstruction_HPGe_Shipping_Container()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction_HPGe_Shipping_Container::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 500*cm;
  G4double world_sizeZ  = 500*cm;
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
    
        //Defining Stainless Steel (Corten)
    G4Element* aluminium = nist->FindOrBuildElement("Al");
    G4Element* carbon = nist->FindOrBuildElement("C");
    G4Element* manganese = nist->FindOrBuildElement("Mn");
    G4Element* phosphorus = nist->FindOrBuildElement("P");
    G4Element* silicon = nist->FindOrBuildElement("S");
    G4Element* copper = nist->FindOrBuildElement("Cu");
    G4Element* chromium = nist->FindOrBuildElement("Cr");
    G4Element* nickel = nist->FindOrBuildElement("Ni");
    G4Element* iron = nist->FindOrBuildElement("Fe");
    
    G4Material* stainless_steel = new G4Material("stainless_steel", 8.04 * kg / m3, 9, kStateSolid);
    stainless_steel->AddElement(aluminium,  0.015/26.981539* perCent);
    stainless_steel->AddElement(carbon,  0.12/12* perCent);
    stainless_steel->AddElement(manganese,  0.2/54.938044* perCent);
    stainless_steel->AddElement(silicon,  0.030/28.0855* perCent);
    stainless_steel->AddElement(phosphorus,  0.07/30.973762* perCent);
    stainless_steel->AddElement(copper,  0.25/63.546* perCent);
    stainless_steel->AddElement(chromium,  0.5/51.9961* perCent);
    stainless_steel->AddElement(nickel,  0.65/58.6934* perCent);
    
    G4double percent_iron = 100 - (0.015/26.981539 + 0.12/12 + 0.2/54.938044 + 0.030/28.0855 + 0.07/30.973762 + 0.25/63.546 +  0.5/51.9961 + 0.65/58.6934);
    stainless_steel->AddElement(iron,  percent_iron* perCent);


    
    //Shipping Container Volume
    
    G4double container_sizeX = 6.05*m;
    G4double container_sizeZ = 2.44*m;
    G4double container_sizeY = 2.59*m;
    
    G4Material* container_mat = nist->FindOrBuildMaterial("stainless_steel");
    G4Box* containerWalls = new G4Box("container",0.5*container_sizeX, 0.5*container_sizeY, 0.5*container_sizeZ);     //its size
    
    // Cavity Volume
    
    G4double cavity_sizeX = 5.44*m;
    G4double cavity_sizeZ = 2.28*m;
    G4double cavity_sizeY = 2.26*m;
    
    G4Box* cavityWalls = new G4Box("cavity",0.5*cavity_sizeX, 0.5*cavity_sizeY, 0.5*cavity_sizeZ); 
   // G4LogicalVolume* cavity_log = new G4LogicalVolume(cavityWalls, "G4_AIR", "cavity"); //its name
    
    G4VSolid* container_total = new G4SubtractionSolid("container_total", containerWalls, cavityWalls);
    G4ThreeVector container_pos = G4ThreeVector(0,0,0);
    
    
    // Creating the Container w/o Wooden Floor
    
    // G4VPhysicalVolume* wall_physical = new G4PVPlacement(0, wall_pos, "wall", logicwalls, physWorld, false, 0, checkOverlaps);
    G4LogicalVolume* container_log = new G4LogicalVolume(container_total, container_mat, "container_total");
    G4VisAttributes* ContainerVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    ContainerVisAtt->SetForceSolid(false);
    container_log->SetVisAttributes(ContainerVisAtt);
    
    G4VPhysicalVolume* container_physical = new G4PVPlacement(0, container_pos, container_log, "container_total", logicWorld, false, 0, checkOverlaps);	
    
    // Wooden Floor Volume
    
    G4double table_sizeX = 544*cm;
    G4double table_sizeY = 16.5*cm;
    G4double table_sizeZ  = 228*cm;
    G4Material* table_mat = nist->FindOrBuildMaterial("G4_CELLULOSE_BUTYRATE");
    G4Box* tableWalls = new G4Box("table",0.5*table_sizeX, 0.5*table_sizeY, 0.5*table_sizeZ);     //its size
    G4LogicalVolume* table_log = new G4LogicalVolume(tableWalls, table_mat, "table"); //its name
    
    
    G4double table_shift = - (cavity_sizeY + table_sizeY)*0.5;
    G4ThreeVector table_pos = G4ThreeVector(0, table_shift,0);
    
    G4VPhysicalVolume* table_physical = new G4PVPlacement(0, table_pos, table_log, "table", container_log, false, 0, checkOverlaps);	
    G4VisAttributes* tableVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
    tableVisAtt->SetForceSolid(true);
    table_log->SetVisAttributes(tableVisAtt);
    
    // Materials
    G4Material* Al_mat = nist->FindOrBuildMaterial("G4_Al");
    G4Material* vacuum_mat = nist->FindOrBuildMaterial("G4_Galactic");
    G4Material* crystal_mat = nist->FindOrBuildMaterial("G4_Ge");
	G4Material* window_mat = nist->FindOrBuildMaterial("G4_C");
	G4Material* glass_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

    // Aliminium housing parameters
    G4double housing_rad = (67.00/2)*mm;
    G4double housing_hz = (243.00/2)*mm;

    // Vacuum parameters
	G4double vacuum_rad = (66.00/2)*mm;  // accounting for 0.5 mm Al walls
	G4double vacuum_hz = (241.90/2)*mm;  // accounting for 0.5 mm Al wall at the back and 0.60 mm window at the front

	// Germanium crystal parameters
	G4double crystal_rad = (50.90/2)*mm;
	G4double crystal_hz = (20.00/2)*mm;
	G4double cavity_rad = (10.00/2)*mm;  // guesstimate
	G4double cavity_hz = (3.00/2)*mm;  // guesstimate

    // Window parameters
	G4double window_rad = (67.00/2)*mm;  // same as crystal (guess)
	G4double window_hz = (0.60/2)*mm;

	// Top glass parametersta
	G4double top_glass_x = (220.00/2)*mm;
	G4double top_glass_y = (200.00/2)*mm;
	G4double top_glass_z = (10.00/2)*mm;
	G4double hole_rad = (80.00/2)*mm;
	
	// Side glass parameter
	G4double side_glass_x = (10.00/2)*mm;
	G4double side_glass_y = (200.00/2)*mm;
	G4double side_glass_z = (258.00/2)*mm;
	
	G4double glass_end_x = (109.00/2)*mm;
	G4double glass_end_y = (200.00/2)*mm;
	G4double glass_end_z = (10.00/2)*mm;
	
	// Lollipop parameter
	G4double lollipop_house_x = (29.50/2)*mm;
	G4double lollipop_house_y = (39.00/2)*mm;
	G4double lollipop_house_z = (9.00/2)*mm;
	
	G4double lollipop_leg_x = (9.00/2)*mm;
	G4double lollipop_leg_y = (139.00/2)*mm;
	G4double lollipop_leg_z = (5.00/2)*mm;
	
    // ALuminium housing logical volume
	G4Tubs* housing_shape = new G4Tubs("Housing_walls", 0, housing_rad, housing_hz, 0*deg, 360*deg);
	G4LogicalVolume* housing_log = new G4LogicalVolume(housing_shape, Al_mat, "Housing_Log");
	
    // Vacuum in housing logical volume
	G4Tubs* vacuum_shape = new G4Tubs("Vacuum", 0, vacuum_rad, vacuum_hz, 0*deg, 360*deg);
	G4LogicalVolume* vacuum_log = new G4LogicalVolume(vacuum_shape, vacuum_mat, "VacuumLog");
    
    // Germanium crytal logical volume
    G4Tubs* crystal_tub_shape = new G4Tubs("Crystal_Tub", 0, crystal_rad, crystal_hz, 0*deg, 360*deg);
    G4Tubs* cavity_shape = new G4Tubs("Cavity", 0, cavity_rad, cavity_hz, 0*deg, 360*deg);
    G4VSolid* crystal_shape = new G4SubtractionSolid("Crystal", crystal_tub_shape, cavity_shape, 0, G4ThreeVector(0, 0, -crystal_hz+cavity_hz));
    G4LogicalVolume* crystal_log = new G4LogicalVolume(crystal_shape, crystal_mat, "Crystal_Log");
	
	// Window logical volume
	G4Tubs* window_shape = new G4Tubs("Window", 0, window_rad, window_hz, 0*deg, 360*deg);
	G4LogicalVolume* window_log = new G4LogicalVolume(window_shape, window_mat, "Window_Log");
	
	// Top glass logical volume
	G4Box* glass_plate_shape = new G4Box("Glass_Pate", top_glass_x, top_glass_y, top_glass_z);
	G4Tubs* hole_shape = new G4Tubs("Hole", 0, hole_rad, top_glass_z+1*cm, 0*deg, 360*deg);
    G4VSolid* top_glass_shape = new G4SubtractionSolid("Top_Glass", glass_plate_shape, hole_shape);
	G4LogicalVolume* top_glass_log = new G4LogicalVolume(top_glass_shape, glass_mat, "Top_Glass_Log");
	
	// Side glass logical volume
	G4Box* side_plate_shape = new G4Box("Side_plate", side_glass_x, side_glass_y, side_glass_z);
	G4Box* end_shape = new G4Box("Glass_end", glass_end_x, glass_end_y, glass_end_z);
	G4VSolid* side_glass_shape = new G4UnionSolid("Side_glass", side_plate_shape, end_shape, 0, G4ThreeVector(-glass_end_x+side_glass_x, 0, -side_glass_z+glass_end_z));
	G4LogicalVolume* side_glass_log = new G4LogicalVolume(side_glass_shape, glass_mat, "Side_Glass_Log");
	
	// Lollipop logical volume
	G4Box* lollipop_housing_shape = new G4Box("Lollipop_housing", lollipop_house_x, lollipop_house_y, lollipop_house_z);
	G4Box* lollipop_leg_shape = new G4Box("Lollipop_leg", lollipop_leg_x, lollipop_leg_y, lollipop_leg_z);
	G4VSolid* lollipop_shape = new G4UnionSolid("Lollipop", lollipop_housing_shape, lollipop_leg_shape, 0, G4ThreeVector(0, lollipop_house_y+lollipop_leg_y, 0));
	G4LogicalVolume* lollipop_log = new G4LogicalVolume(lollipop_shape, glass_mat, "Lollipop_log");
	
	// Set visualization attributes for the materials
	G4VisAttributes* housingVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
	// housingVisAtt->SetForceSolid(true);
    housing_log->SetVisAttributes(housingVisAtt);
    
    G4VisAttributes* vacuumVisAtt = new G4VisAttributes(G4Colour());
    // vacuumVisAtt->SetForceSolid(true);
    vacuum_log->SetVisAttributes(vacuumVisAtt);
    
    G4VisAttributes* crystalVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    crystalVisAtt->SetForceSolid(true);
    crystal_log->SetVisAttributes(crystalVisAtt);
    
    G4VisAttributes* windowVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
    windowVisAtt->SetForceSolid(true);
    window_log->SetVisAttributes(windowVisAtt);
    
    G4VisAttributes* glassVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    glassVisAtt->SetForceSolid(true);
    top_glass_log->SetVisAttributes(glassVisAtt);
    side_glass_log->SetVisAttributes(glassVisAtt);
    
    G4VisAttributes* lollipopVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    lollipopVisAtt-> SetForceSolid(true);
    lollipop_log->SetVisAttributes(lollipopVisAtt);

	// Place
	G4double window_pos = housing_hz - window_hz;
	G4double crys_pos = window_pos - window_hz - 5.50*mm - crystal_hz;
	G4double top_glass_pos = housing_hz + 32.00*mm + top_glass_z;
	G4double side_glass_pos_x = top_glass_x + side_glass_x;
	G4double side_glass_pos_z = top_glass_pos-side_glass_z+top_glass_z;
    G4double lollipop_pos = top_glass_pos + 2*top_glass_z + lollipop_house_z;
    
    G4double offset = 200.00*mm;

	
	G4RotationMatrix* side_glass_rot = new G4RotationMatrix();
	side_glass_rot->rotateZ(180*deg);
	
    new G4PVPlacement(0, G4ThreeVector(0, 0, offset), housing_log, "Al_Tube", logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -0.05*mm), vacuum_log, "Vacuum", housing_log, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, crys_pos), crystal_log, "Crystal", vacuum_log, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, window_pos), window_log, "Window", housing_log, false, 0, checkOverlaps);  // z is tube length - window thickness
    new G4PVPlacement(0, G4ThreeVector(0, 0, offset + top_glass_pos), top_glass_log, "Top_Glass", logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(-side_glass_pos_x, 0, offset + side_glass_pos_z), side_glass_log, "Side_Glass_A", logicWorld, false, 0, checkOverlaps);
    new G4PVPlacement(side_glass_rot, G4ThreeVector(side_glass_pos_x, 0, offset + side_glass_pos_z), side_glass_log, "Side_Glass_B", logicWorld, false, 1, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, offset+lollipop_pos), lollipop_log, "Lollipop", logicWorld, false, 0, checkOverlaps);

 //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
