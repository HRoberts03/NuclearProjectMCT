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

#include "B1DetectorConstruction_Shipping_Container.hh"

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
#include <G4VisAttributes.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction_Shipping_Container::B1DetectorConstruction_Shipping_Container()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction_Shipping_Container::~B1DetectorConstruction_Shipping_Container()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction_Shipping_Container::Construct()
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
  //yellow->SetForceSolid(true);
  
  G4double world_sizeXY = 2000*cm;
  G4double world_sizeZ  = 2000*cm;
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
    
                                  
    // Air
  //
  G4double air_sizeXY = 1000*cm;
  G4double air_sizeZ  = 1000*cm;
  G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidAir = new G4Box("air",0.5*air_sizeXY, 0.5*air_sizeXY, 0.5*air_sizeZ);     //its size
      
  G4LogicalVolume* airWorld = new G4LogicalVolume(solidAir, air_mat, "inner"); //its name
  
  G4ThreeVector air_pos = G4ThreeVector(0,0,0);
  G4VPhysicalVolume* air_physical = new G4PVPlacement(0, air_pos, "air", airWorld, physWorld, false, 0, checkOverlaps);	
    
	//face construction 
	G4double face_r = 5.71*cm;
	G4double back_hz =  2*mm;
	G4double face_hz =  0.5*mm;
	G4Tubs* face_shape = new G4Tubs("face", 0*mm, face_r, face_hz, 0*deg, 360*deg);
	G4LogicalVolume* face_log = new G4LogicalVolume(face_shape, nist->FindOrBuildMaterial("G4_Al"), "face");
	face_log->SetVisAttributes(red);


	// scintillator construction
	
	// define scintillator crystal properties
	
	G4Material* NaI = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    G4Element* Thallium = nist->FindOrBuildElement("Tl");
	
	G4Material* NaTi = new G4Material("NaTi", 3.67*g/cm3, 2, kStateSolid);
    NaTi->AddMaterial(NaI, 99.5*perCent);
    NaTi->AddElement(Thallium, 0.5*perCent);
    
    // define the crystal logical volume
	G4double scint_r = 5.51*cm;
	G4double scint_hz = 5.40*cm;
	G4Tubs* scintillator_shape = new G4Tubs("scintillator",0*mm,scint_r,scint_hz,0*deg,360*deg);
	G4LogicalVolume* scint_log = new G4LogicalVolume(scintillator_shape, NaTi, "scintillator");
    scint_log->SetVisAttributes(blue);
    // Shield Construction
    G4double shield_r = 5.71*cm;
	G4double shield_hz = 5.40*cm;
	G4Tubs* shield_shape = new G4Tubs("shield",0,shield_r,shield_hz,0*deg,360*deg);
	G4LogicalVolume* shield_log = new G4LogicalVolume(shield_shape, nist->FindOrBuildMaterial("G4_Al"), "shield");
	shield_log->SetVisAttributes(yellow);
	
	//Detector Back Construction
	G4double back_r = 5.71*cm;
	G4Tubs* back_shape = new G4Tubs("back", 0*mm, back_r, back_hz, 0*deg, 360*deg);
	G4LogicalVolume* back_log = new G4LogicalVolume(back_shape, nist->FindOrBuildMaterial("G4_Al"), "face");
	back_log->SetVisAttributes(green);


	
	// Apparatus Shifts
   
    G4double hz_shift = 10*cm;
    G4double shipping_shift = 1.25*m;
    
    

	// shield Placement
	G4ThreeVector shield_pos = G4ThreeVector(0,0,hz_shift+shield_hz+face_hz+shipping_shift);
	G4RotationMatrix* shield_rot = new G4RotationMatrix();
	shield_rot -> rotateX(0*deg);
	G4VPhysicalVolume* shield = new G4PVPlacement(0, shield_pos, "shield", shield_log, air_physical, false, 0, checkOverlaps);
    
    // Scintillator Placement
	G4ThreeVector scintillator_pos = G4ThreeVector(0,0,shipping_shift);
	G4RotationMatrix* scintillator_rot = new G4RotationMatrix();
	scintillator_rot -> rotateX(0*deg);
     new G4PVPlacement(0, scintillator_pos, "Crystal", scint_log, shield, false, 0, checkOverlaps);

    
    // back placement
    G4ThreeVector back_pos = G4ThreeVector(0,0,shield_hz+back_hz+shipping_shift);
	G4RotationMatrix* back_rot = new G4RotationMatrix();
	back_rot->rotateX(0*deg);
     new G4PVPlacement(0, back_pos, "back", back_log, shield, false, 0, checkOverlaps);	
    
    // face placement
	G4ThreeVector face_pos = G4ThreeVector(0,0,-1*shield_hz-face_hz+shipping_shift);
	G4RotationMatrix* face_rot = new G4RotationMatrix();
	face_rot->rotateX(0*deg);
    new G4PVPlacement(0, face_pos, "face", face_log, shield, false, 0, checkOverlaps);	



 //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
