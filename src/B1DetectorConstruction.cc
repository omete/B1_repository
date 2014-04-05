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
// $Id: B1DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
// Added by O.Mete
#include "G4ParticleGun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 200*mm, env_sizeZ = 50*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.0*env_sizeXY;
  G4double world_sizeZ  = 1.0*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeXY, world_sizeXY, world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        env_sizeXY, env_sizeXY, env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Plasma Section
  // 
  // Material definition by using NIST database 
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Li");

  // ---> Material definition by hand
  // //////////////////////////////////////////////////////////////
  G4double a; 
  G4double z;
  G4double np;
  G4double density;
  G4double avo_num;
  G4double vol_p;  // plasma volume 
  
  // Calculate the gas density
  //a=85.4678*g/mole; // Rb
  //z=37; // Rb
  a=6.941*g/mole; // Li
  z = 3; // Li
  avo_num = 6.02e23;
  np = 6e14*cm-3;       // plasma density, number of ions
  //vol_p = 3141.6*cm3;

  density = (np*a)/avo_num;

  new G4Material("mygas", z, a, density, kStateGas);
     //mygas can always be updated for Rb properties in gas state.

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // Get materials
  G4Material* plasma = G4Material::GetMaterial("mygas");


  // //////////////////////////////////////////////////////////////
  // ---> Material add until here

    
  // //////////////////////////////////////////////////////////////
  // ---> Start defining plasma section here in subsections
 
  // Cylinder section 1 /////////////////////////////////////////////
    
  G4ThreeVector pos1 = G4ThreeVector(0, 0, 25.*m);
        
  G4double shape1_rmin =  0.*mm;
  G4double shape1_rmax =  100.*mm;
  G4double shape1_hz = 25.*m;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Tubs* solidShape1 = new G4Tubs("Shape1", 
                                    shape1_rmin, 
                                    shape1_rmax, 
                                    shape1_hz,
                                    shape1_phimin, 
                                    shape1_phimax);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                            plasma,          //its material 
                           "Shape1");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
    // Cylinder section 2 /////////////////////////////////////////////
    /*
    G4ThreeVector pos2 = G4ThreeVector(0, 0, 75.*m);
    
    G4Tubs* solidShape2 = new G4Tubs("Shape2",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape2 =
    new G4LogicalVolume(solidShape2,         //its solid
                        plasma,          //its material
                        "Shape2");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos2,                    //at position
                      logicShape2,             //its logical volume
                      "Shape2",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    // Cylinder section 3 /////////////////////////////////////////////
    G4ThreeVector pos3 = G4ThreeVector(0, 0, 125.*m);
    
    G4Tubs* solidShape3 = new G4Tubs("Shape3",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape3 =
    new G4LogicalVolume(solidShape3,         //its solid
                        plasma,          //its material
                        "Shape3");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos3,                    //at position
                      logicShape3,             //its logical volume
                      "Shape3",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    // Cylinder section 4 /////////////////////////////////////////////
    G4ThreeVector pos4 = G4ThreeVector(0, 0, 175.*m);

    G4Tubs* solidShape4 = new G4Tubs("Shape4",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape4 =
    new G4LogicalVolume(solidShape4,         //its solid
                        plasma,          //its material
                        "Shape4");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos4,                    //at position
                      logicShape4,             //its logical volume
                      "Shape4",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    // Cylinder section 5 /////////////////////////////////////////////
    G4ThreeVector pos5 = G4ThreeVector(0, 0, 225.*m);
    
    G4Tubs* solidShape5 = new G4Tubs("Shape5",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape5 =
    new G4LogicalVolume(solidShape5,         //its solid
                        plasma,          //its material
                        "Shape5");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos5,                    //at position
                      logicShape5,             //its logical volume
                      "Shape5",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

    
    // Cylinder section 6 /////////////////////////////////////////////
    G4ThreeVector pos6 = G4ThreeVector(0, 0, 275.*m);
    
    G4Tubs* solidShape6 = new G4Tubs("Shape6",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape6 =
    new G4LogicalVolume(solidShape6,         //its solid
                        plasma,          //its material
                        "Shape6");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos6,                    //at position
                      logicShape6,             //its logical volume
                      "Shape6",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

    // Cylinder section 7 /////////////////////////////////////////////
    G4ThreeVector pos7 = G4ThreeVector(0, 0, 325.*m);
    
    G4Tubs* solidShape7 = new G4Tubs("Shape7",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape7 =
    new G4LogicalVolume(solidShape7,         //its solid
                        plasma,          //its material
                        "Shape7");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos7,                    //at position
                      logicShape7,             //its logical volume
                      "Shape7",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

    // Cylinder section 8 /////////////////////////////////////////////
    G4ThreeVector pos8 = G4ThreeVector(0, 0, 375.*m);
    
    G4Tubs* solidShape8 = new G4Tubs("Shape8",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape8 =
    new G4LogicalVolume(solidShape8,         //its solid
                        plasma,          //its material
                        "Shape8");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos8,                    //at position
                      logicShape8,             //its logical volume
                      "Shape8",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    // Cylinder section 9 /////////////////////////////////////////////
    G4ThreeVector pos9 = G4ThreeVector(0, 0, 425.*m);
    
    G4Tubs* solidShape9 = new G4Tubs("Shape9",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape9 =
    new G4LogicalVolume(solidShape9,         //its solid
                        plasma,          //its material
                        "Shape9");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos9,                    //at position
                      logicShape9,             //its logical volume
                      "Shape9",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

    // Cylinder section 10 /////////////////////////////////////////////
    G4ThreeVector pos10 = G4ThreeVector(0, 0, 475.*m);
    
    G4Tubs* solidShape10 = new G4Tubs("Shape10",
                                     shape1_rmin,
                                     shape1_rmax,
                                     shape1_hz,
                                     shape1_phimin,
                                     shape1_phimax);
    
    G4LogicalVolume* logicShape10 =
    new G4LogicalVolume(solidShape10,         //its solid
                        plasma,          //its material
                        "Shape10");           //its name
    
    new G4PVPlacement(0,                       //no rotation
                      pos10,                    //at position
                      logicShape10,             //its logical volume
                      "Shape10",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  */
    
  fScoringVolume = logicShape1;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
