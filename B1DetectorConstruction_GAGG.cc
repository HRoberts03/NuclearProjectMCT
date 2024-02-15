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
/// \Constructing GAGG detector

#include "B1DetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialTable.hh"
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
  

    // Define GAGG molecule
    G4Element* aluminium = nist->FindOrBuildElement("Al");
    G4Element* gadolinium = nist->FindOrBuildElement("Gd");
    G4Element* gallium = nist->FindOrBuildElement("Ga");
    G4Element* oxygen = nist->FindOrBuildElement("O");
    G4Element* cerium = nist->FindOrBuildElement("Ce");

    G4Material* GAGG = new G4Material("GAGG", 6.67 * g / cm3, 4, kStateSolid);
    GAGG->AddElement(gadolinium,  3);
    GAGG->AddElement(aluminium, 2);
    GAGG->AddElement(gallium, 2);
    GAGG->AddElement(oxygen, 12);
    
    // Define GAGGCe crystal
    G4Material* GAGGCe = new G4Material("GAGGCe", 6.67*g/cm3, 2, kStateSolid);
    GAGGCe->AddMaterial(GAGG, 99*perCent);
    GAGGCe->AddElement(cerium, 1*perCent);

    // Define entrance window material
    G4Material* window_mat = nist->FindOrBuildMaterial("G4_MYLAR");

    // Define Mu metal
    G4Element* nickel = nist->FindOrBuildElement("Ni");
    G4Element* molybdenum = nist->FindOrBuildElement("Mo");
    G4Element* manganese = nist->FindOrBuildElement("Mn");
    G4Element* silicon = nist->FindOrBuildElement("Si");
    G4Element* cobalt = nist->FindOrBuildElement("Co");
    G4Element* copper = nist->FindOrBuildElement("Cu");
    G4Element* chromium = nist->FindOrBuildElement("Cr");
    G4Element* phosphorous = nist->FindOrBuildElement("P");
    G4Element* sulfur = nist->FindOrBuildElement("S");
    G4Element* carbon = nist->FindOrBuildElement("C");
    G4Element* iron = nist->FindOrBuildElement("Fe");

    G4Material* mu_metal = new G4Material("Mu_metal", 8.7 * g / cm3, 11, kStateSolid);
    mu_metal->AddElement(nickel,  80* perCent);
    mu_metal->AddElement(molybdenum, 4.5* perCent);
    mu_metal->AddElement(manganese, 0.8* perCent);
    mu_metal->AddElement(silicon, 0.5* perCent);
    mu_metal->AddElement(cobalt, 0.5 * perCent);
    mu_metal->AddElement(copper, 0.3 * perCent);
    mu_metal->AddElement(chromium, 0.3 * perCent);
    mu_metal->AddElement(phosphorous, 0.02 * perCent);
    mu_metal->AddElement(sulfur, 0.01 * perCent);
    mu_metal->AddElement(carbon, 0.05 * perCent);
    mu_metal->AddElement(iron, 13.02* perCent);

    // Define PMT material
    G4Material* pmt_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    
    G4MaterialPropertiesTable* pmt_properties = new G4MaterialPropertiesTable();
    pmt_properties->AddConstProperty("EFFICIENCY", 0.26);
    pmt_mat->SetMaterialPropertiesTable(pmt_properties);
    
    // Define vacuum material
    G4Material* vacuum_mat = nist->FindOrBuildMaterial("G4_Galactic");
    
    // Define aluminium material
    G4Material* Al_mat = nist->FindOrBuildMaterial("G4_Al");
    
    // Define air material
    G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
    
    // Air parameters
    G4double air_rad = 45.00/2.*mm;
    G4double air_hz = 140.00/2.*mm;
    
    // Crystal parameters
    G4double crystal_rad = 18.00/2.*mm;
    G4double crystal_hz = 25.00/2.*mm;

    // Crystal housing (top-bottom of these volumes are the "front-back" of the detector)
    G4double crystalAl_outer_rad = (23.00/2.)*mm;
    G4double crystalAl_hz = (20.50/2.)*mm;

    // Entrance window parameters
    G4double window_rad = (18.00/2.)*mm;
    G4double window_hz = (18.00/2.)*um;
    
    // Al junction between magnetic metal shield and crystal housing
    G4double junction_outer_rad = 39.00/2*mm;
    G4double junction_hz = 5.50/2.*mm;  // measured

    // Magnetic Sheild parameters
    G4double shield_inner_rad = ((39.00/2)-0.64)*mm;  // accounting for thickness of shield walls (0.64mm)
    G4double shield_outer_rad = 39.00/2*mm;
    G4double shield_hz = 87.00/2*mm;  // measured
    
    // Photomultiplier tube parameters
    G4double pmt_rad = 28.50/2*mm;
    G4double pmt_hz = 60.00/2*mm;

    // Vacuum inside pmt parameters
    G4double pmt_vacuum_rad = 27.50/2*mm;
    G4double pmt_vacuum_hz = 59.00/2*mm;

    // Bottom Aluminium Housing
    G4double BottomAl_inner_rad = 38.00/2*mm;
    G4double BottomAl_outer_rad = 39.00/2*mm;
    G4double BottomAl_hz = 25.00/2*mm;  // ((134-(87+20.5+1))-0.5), 0.5 subtracted to acocunt for thickness of end piece
    G4double BottomAl_end_hz = 0.50/2*mm;
    
    // Air logical volume
    G4Tubs* air_shape = new G4Tubs("Air_shape", 0, air_rad, air_hz, 0*deg, 360*deg);
    G4LogicalVolume* air_log = new G4LogicalVolume(air_shape, air_mat, "Air_Log");

    // Crystal housing logical volume
    G4Tubs* crys_housing = new G4Tubs("Crys_Housing", crystal_rad, crystalAl_outer_rad, crystalAl_hz, 0*deg, 360*deg);
    G4LogicalVolume* crys_housing_log = new G4LogicalVolume(crys_housing, Al_mat, "Crys_Housing_Log");
    
    // Entrance window logical volume
    G4Tubs* window_shape = new G4Tubs("Window", 0, window_rad, window_hz, 0*deg, 360*deg);
    G4LogicalVolume* window_log = new G4LogicalVolume(window_shape, window_mat, "Mylar_Entrance_Window");
    
    // Crystal logical volume
    G4Tubs* crystal_shape = new G4Tubs("Crystal", 0, crystal_rad, crystal_hz, 0*deg, 360*deg);
    G4LogicalVolume* crystal_log = new G4LogicalVolume(crystal_shape, GAGGCe, "GAGG_Ce_Crystal");
    
    // Al junction logical volume
    G4Tubs* junction_shape = new G4Tubs("Junction", crystal_rad, junction_outer_rad, junction_hz, 0*deg, 360*deg);
    G4LogicalVolume* junction_log = new G4LogicalVolume(junction_shape, Al_mat, "Junction");

    // Magnetic Shield logical volume
    G4Tubs* shield_shape = new G4Tubs("Shield", shield_inner_rad, shield_outer_rad, shield_hz, 0*deg, 360*deg);
    G4LogicalVolume* shield_log = new G4LogicalVolume(shield_shape, mu_metal, "Metal_Shield");

    // Photomultiplier tube logical volume
    G4Tubs* pmt_shape = new G4Tubs("PMT_walls", 0, pmt_rad, pmt_hz, 0*deg, 360*deg);
    G4LogicalVolume* pmt_log = new G4LogicalVolume(pmt_shape, pmt_mat, "PMT");

    // Vacuum inside pmt logical volume
    G4Tubs* vacuum_shape = new G4Tubs("Vacuum_shape", 0, pmt_vacuum_rad, pmt_vacuum_hz, 0*deg, 360*deg);
    G4LogicalVolume* vacuum_log = new G4LogicalVolume(vacuum_shape, vacuum_mat, "Vacuum");

    // Bottom Al housing logical volume
    G4Tubs* al_housing_walls_shape = new G4Tubs("Bottom_Al_Walls", BottomAl_inner_rad, BottomAl_outer_rad, BottomAl_hz, 0*deg, 360*deg);
    G4Tubs* al_housing_end_shape = new G4Tubs("Bottom_Al_end", 0, BottomAl_outer_rad, BottomAl_end_hz, 0*deg, 360*deg);
    G4VSolid* bottom_al_shape = new G4UnionSolid("Bottom_Al_housing", al_housing_walls_shape, al_housing_end_shape, 0, G4ThreeVector(0., 0., (-BottomAl_hz-BottomAl_end_hz)*mm));
    G4LogicalVolume* al_bottom_hous_log = new G4LogicalVolume(bottom_al_shape, Al_mat, "Al_bottom");

    // Set visualization attributes for the materials
    G4VisAttributes* crys_housingVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    crys_housingVisAtt->SetForceSolid(false);
    crys_housing_log->SetVisAttributes(crys_housingVisAtt);
    
    G4VisAttributes* windowVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    windowVisAtt->SetForceSolid(true);
    window_log->SetVisAttributes(windowVisAtt);

    G4VisAttributes* crystalVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    crystalVisAtt->SetForceSolid(true);
    crystal_log->SetVisAttributes(crystalVisAtt);
    
    G4VisAttributes* junctionVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    junctionVisAtt->SetForceSolid(true);
    junction_log->SetVisAttributes(junctionVisAtt);
    
    G4VisAttributes* shieldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    // shieldVisAtt->SetForceSolid(true);
    shield_log->SetVisAttributes(shieldVisAtt);
    
    G4VisAttributes* pmtVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    // pmtVisAtt->SetForceSolid(true);
    pmt_log->SetVisAttributes(pmtVisAtt);

    G4VisAttributes* vacuumVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    vacuumVisAtt->SetForceSolid(true);
    vacuum_log->SetVisAttributes(vacuumVisAtt);
    
    G4VisAttributes* al_bottomVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    al_bottomVisAtt->SetForceSolid(true);
    al_bottom_hous_log->SetVisAttributes(al_bottomVisAtt);
    
    G4VisAttributes* airVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    airVisAtt->SetForceWireframe(true);
    air_log->SetVisAttributes(airVisAtt);


	// Place
	G4double housing_shift = (air_hz-crystalAl_hz);
	G4double crystal_pos = housing_shift + crystalAl_hz - crystal_hz -1.00*mm;
	G4double window_pos = crystal_pos + crystal_hz + window_hz;
	G4double junction_pos = housing_shift-(crystalAl_hz+junction_hz);
	G4double shield_pos = junction_pos-(junction_hz+shield_hz);
	G4double pmt_pos = crystal_pos-crystal_hz-pmt_hz;
	G4double bottom_pos = shield_pos-(shield_hz+BottomAl_hz);

	new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), air_log, "Air_Parent", logicWorld, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(0., 0., crystal_pos), crystal_log, "Crystal", air_log, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(0., 0., housing_shift), crys_housing_log, "Crystal_Housing", air_log, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0., 0., window_pos), window_log, "Window", air_log, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(0., 0., junction_pos), junction_log, "Junction", air_log, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(0., 0., shield_pos), shield_log, "Magnetic_Shield", air_log, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(0., 0., pmt_pos), pmt_log, "PMT", air_log, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vacuum_log, "Vacuum_in_PMT", pmt_log, false, 0, checkOverlaps);
	new G4PVPlacement(0, G4ThreeVector(0., 0., bottom_pos), al_bottom_hous_log, "Bottom_Housing", air_log, false, 0, checkOverlaps);

    //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
