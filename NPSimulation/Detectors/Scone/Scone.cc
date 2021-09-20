/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : March 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Scone simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
#include <iostream>
#include <string>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SubtractionSolid.hh"

// NPTool header
#include "Scone.hh"
#include "CalorimeterScorers.hh"
#include "ProcessScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Scone_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.*MeV;
  const double ResoTime = 0.5*ns ;
  const double ResoEnergy = 0.001*MeV ;
  const double XSection = 25.1*mm ; 
  const double YSection = 25.6*mm ;
  const double LengthR1 = 1000*mm ;
  const double LengthR2 = 500*mm ;
  const double Length2x2 = 400*mm ;
  const double GdThickness = 25*um ;
  const double DistanceBetweenBar = 1*mm ;
  const string Material = "EJ200";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Scone Specific Method
Scone::Scone(){
  InitializeMaterial();
  m_Event = new TSconeData() ;
  m_SconeScorer = 0;
  m_GdScorer = 0;
  m_FCScorer = 0;

  m_SquareDetector = 0;
  m_FissionChamberVolume = 0;

  m_BuildRing1 = 1;
  m_BuildRing2 = 1;
  m_BuildFissionChamber = 0;

  m_NumberOfInnerDetector = 16;
  m_NumberOfRing1Detector = 8;
  m_NumberOfRing2Detector = 16;

  m_Assembly = 0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));   
  m_Vis2x2 = new G4VisAttributes(G4Colour(0, 0.5, 0.6, 1.0));   
  m_Vis6x6R1 = new G4VisAttributes(G4Colour(0.2, 1, 1, 1));   
  m_Vis6x6R2 = new G4VisAttributes(G4Colour(0.2, 0.8, 0.6, 1));   
  m_VisGd = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 1));   
}

Scone::~Scone(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Scone::Build2x2Assembly(int DetNumber){
  
  G4String Name = "2x2block_" + to_string(DetNumber);

  //--- Full assembly definition ---//
  double Xsize = 2*Scone_NS::XSection + 2*2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar;
  double Ysize = 2*Scone_NS::YSection + 2*2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar;
  double Zsize = Scone_NS::Length2x2;

  G4Box* solidMotherVolume = new G4Box("2x2block", 0.5*Xsize, 0.5*Ysize, 0.5*Zsize);
  m_2x2Assembly = new G4LogicalVolume(solidMotherVolume, m_MaterialVaccuum, Name, 0, 0, 0);

  G4VisAttributes * MotherVolumeVisAtt = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
  MotherVolumeVisAtt->SetForceWireframe(true);
  //m_2x2Assembly->SetVisAttributes(G4VisAttributes::Invisible);
  m_2x2Assembly->SetVisAttributes(MotherVolumeVisAtt);

  //--- One bar definitation ---//
  G4Box* box = new G4Box("Scone_Box",Scone_NS::XSection*0.5,
      Scone_NS::YSection*0.5,Scone_NS::Length2x2*0.5);
  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Scone_NS::Material);

  m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_Scone_Box_",0,0,0);
  m_SquareDetector->SetVisAttributes(m_Vis2x2);
  m_SquareDetector->SetSensitiveDetector(m_SconeScorer);

  //--- Gd Layer around the plastic bar ---//
  G4Material* Gd_Material = MaterialManager::getInstance()->GetMaterialFromLibrary("Gd");
  G4Box* box_tmp = new G4Box("Box_tmp",(Scone_NS::XSection+Scone_NS::GdThickness)*0.5, (Scone_NS::YSection+Scone_NS::GdThickness)*0.5, Scone_NS::Length2x2*0.49);
  G4VSolid* Gd_layer = (G4VSolid*) new G4SubtractionSolid("Gd_layer", box_tmp, box,0,G4ThreeVector(0,0,0));

  G4LogicalVolume* Gd_layer_volume = new G4LogicalVolume(Gd_layer, Gd_Material,"Gd_layer_volume",0,0,0);
  Gd_layer_volume->SetVisAttributes(m_VisGd);
  Gd_layer_volume->SetSensitiveDetector(m_GdScorer);

  double posX = 0;
  double posY = 0;
  double posX_orig = -0.5*(Scone_NS::XSection + Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar);
  double posY_orig = -0.5*(Scone_NS::YSection + Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar);
  double BarNumber = 0;
  G4ThreeVector Tv;
  Tv.setX(0); Tv.setY(0); Tv.setZ(0);
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      BarNumber++;

      posX= posX_orig+i*(Scone_NS::XSection + 2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar);
      posY= posY_orig+j*(Scone_NS::YSection + 2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar);
      Tv.setX(posX);
      Tv.setY(posY);
      
      new G4PVPlacement(new G4RotationMatrix(0,0,0),
          Tv,
          m_SquareDetector,"PlasticBar",
          m_2x2Assembly, false, BarNumber);
      
      new G4PVPlacement(new G4RotationMatrix(0,0,0),
          Tv,
          Gd_layer_volume,"GdLayer",
          m_2x2Assembly, false, 0);
    }
  }
  return m_2x2Assembly;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Scone::Build6x6Assembly(int DetNumber, double plastic_length){

  G4String Name = "6x6block_" + to_string(DetNumber);
  
  //--- Full assembly definition ---//
  double Xsize = 6*Scone_NS::XSection + 6*2*Scone_NS::GdThickness + 5*Scone_NS::DistanceBetweenBar;
  double Ysize = 6*Scone_NS::YSection + 6*2*Scone_NS::GdThickness + 5*Scone_NS::DistanceBetweenBar;
  double Zsize = plastic_length;

  G4Box* solidMotherVolume = new G4Box("6x6block", 0.5*Xsize, 0.5*Ysize, 0.5*Zsize);
  m_6x6Assembly = new G4LogicalVolume(solidMotherVolume, m_MaterialVaccuum, Name, 0, 0, 0);

  G4VisAttributes * MotherVolumeVisAtt = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
  MotherVolumeVisAtt->SetForceWireframe(true);
  m_6x6Assembly->SetVisAttributes(MotherVolumeVisAtt);

  //--- One bar definitation ---//
  G4Box* box = new G4Box("Scone_Box",Scone_NS::XSection*0.5,
      Scone_NS::YSection*0.5,plastic_length*0.5);
  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Scone_NS::Material);
  m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_Scone_Box",0,0,0);
  if(plastic_length == Scone_NS::LengthR2)
    m_SquareDetector->SetVisAttributes(m_Vis6x6R2);
  else 
    m_SquareDetector->SetVisAttributes(m_Vis6x6R1);
  m_SquareDetector->SetSensitiveDetector(m_SconeScorer);

  //--- Gd Layer around the plastic bar ---//
  G4Material* Gd_Material = MaterialManager::getInstance()->GetMaterialFromLibrary("Gd");
  G4Box* box_tmp = new G4Box("Box_tmp",(Scone_NS::XSection+Scone_NS::GdThickness)*0.5, (Scone_NS::YSection+Scone_NS::GdThickness)*0.5, plastic_length*0.49);
  G4VSolid* Gd_layer = (G4VSolid*) new G4SubtractionSolid("Gd_layer",box_tmp,box,0,G4ThreeVector(0,0,0));

  G4LogicalVolume* Gd_layer_volume = new G4LogicalVolume(Gd_layer, Gd_Material,"Gd_layer_volume",0,0,0);
  Gd_layer_volume->SetVisAttributes(m_VisGd);
  Gd_layer_volume->SetSensitiveDetector(m_GdScorer);

  double posX = 0;
  double posY = 0;
  double posX_orig=-2.5*Scone_NS::XSection-4.5*Scone_NS::GdThickness-2.5*Scone_NS::DistanceBetweenBar;
  double posY_orig=-2.5*Scone_NS::YSection-4.5*Scone_NS::GdThickness-2.5*Scone_NS::DistanceBetweenBar;
  double BarNumber = 0;
  G4ThreeVector Tv;
  Tv.setX(0); Tv.setY(0); Tv.setZ(0);
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      BarNumber++;

      posX= posX_orig+i*(Scone_NS::XSection + 2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar);
      posY= posY_orig+j*(Scone_NS::YSection + 2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar);  
      Tv.setX(posX);
      Tv.setY(posY);
      
      new G4PVPlacement(new G4RotationMatrix(0,0,0),
          Tv,
          m_SquareDetector,"PlasticBar",
          m_6x6Assembly, false, BarNumber);
 
      new G4PVPlacement(new G4RotationMatrix(0,0,0),
          Tv,
          Gd_layer_volume,"GdLayer",
          m_6x6Assembly, false, 0);
    }
  }
  return m_6x6Assembly;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Scone::BuildSquareDetector(){
  if(!m_SquareDetector){
    G4Box* box = new G4Box("Scone_Box",Scone_NS::XSection*0.5,
        Scone_NS::YSection*0.5,Scone_NS::LengthR1*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Scone_NS::Material);
    m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_Scone_Box",0,0,0);
    m_SquareDetector->SetVisAttributes(m_VisSquare);
    m_SquareDetector->SetSensitiveDetector(m_SconeScorer);
  }
  return m_SquareDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in Detecor
void Scone::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Scone");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Ring1","Ring2","FissionChamber"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Scone " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      m_BuildRing1 = blocks[i]->GetInt("Ring1");
      m_BuildRing2 = blocks[i]->GetInt("Ring2");
      m_BuildFissionChamber = blocks[i]->GetInt("FissionChamber");
      AddDetector(Pos);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Scone::ConstructDetector(G4LogicalVolume* world){
  Build2x2Block(world);

  if(m_BuildRing1==1)
    BuildRing1(world);

  if(m_BuildRing2==1)
    BuildRing2(world);

  if(m_BuildFissionChamber==1){
    G4RotationMatrix* Rot_FC = new G4RotationMatrix(0,0,0);
    G4ThreeVector Pos_FC = G4ThreeVector(0,0,0);
    BuildFissionChamber()->MakeImprint(world,Pos_FC,Rot_FC,0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Scone::BuildFissionChamber(){
  if(!m_FissionChamberVolume){
    m_FissionChamberVolume = new G4AssemblyVolume();

    G4RotationMatrix* Rv = new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0); Tv.setY(0); Tv.setZ(0);

    // Gas Volume //
    double gas_box_width = 5.*cm;
    double gas_box_height = 5.*cm;
    double gas_box_length = 10.*cm;
    
    G4Box* gas_box = new G4Box("gas_box",0.5*gas_box_width,0.5*gas_box_height,0.5*gas_box_length);
    
    G4Material* CF4 = MaterialManager::getInstance()->GetGasFromLibrary("CF4", 1.1*bar, 293*kelvin);
    G4LogicalVolume* gas_box_logic = new G4LogicalVolume(gas_box, CF4, "gas_box_logic",0,0,0);
  
    G4VisAttributes* Vis_gas = new G4VisAttributes(G4Colour(0.776, 0.662, 0.662, 0.5));
    gas_box_logic->SetVisAttributes(Vis_gas);
    gas_box_logic->SetSensitiveDetector(m_FCScorer);

    m_FissionChamberVolume->AddPlacedVolume(gas_box_logic, Tv, Rv);

    // Al // 
    double Al_box_width = 5.*cm;
    double Al_box_height = 5.*cm;
    double Al_box_length = 0.15*um;

    G4Box* box1 = new G4Box("box1", 0.5*Al_box_width, 0.5*Al_box_height, 0.5*Al_box_length);
    
    G4Material* Al_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4LogicalVolume* Al_box_logic = new G4LogicalVolume(box1, Al_material, "Al_box_logic",0,0,0);
    
    G4VisAttributes* VisAl = new G4VisAttributes(G4Colour(0.839, 0.803, 0.803, 1));
    Al_box_logic->SetVisAttributes(VisAl);

    Tv.setZ(-0.5*gas_box_length-0.5*Al_box_length);
    //m_FissionChamberVolume->AddPlacedVolume(Al_box_logic,Tv,Rv);
    Tv.setZ(0.5*gas_box_length+0.5*Al_box_length);
    //m_FissionChamberVolume->AddPlacedVolume(Al_box_logic,Tv,Rv);
  
    // Kapton // 
    double kapton_box_width = 5.*cm;
    double kapton_box_height = 5.*cm;
    double kapton_box_length = 50*um;

    G4Box* box2 = new G4Box("box2", 0.5*kapton_box_width, 0.5*kapton_box_height, 0.5*kapton_box_length);
    
    G4Material* kapton_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Kapton");
    G4LogicalVolume* kapton_box_logic = new G4LogicalVolume(box2, kapton_material, "kapton_box_logic",0,0,0);
    
    G4VisAttributes* VisKapton = new G4VisAttributes(G4Colour(0.539, 0.503, 0.503, 1));
    kapton_box_logic->SetVisAttributes(VisKapton);

    Tv.setZ(-0.5*gas_box_length-0.5*Al_box_length-0.5*kapton_box_length);
    //m_FissionChamberVolume->AddPlacedVolume(Al_box_logic,Tv,Rv);
    Tv.setZ(0.5*gas_box_length+0.5*Al_box_length+0.5*kapton_box_length);
    //m_FissionChamberVolume->AddPlacedVolume(kapton_box_logic,Tv,Rv);

  }
  return m_FissionChamberVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::Build2x2Block(G4LogicalVolume* world){

  G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
  double posX_orig = 0;
  double posX_neg = posX_orig - 2*(Scone_NS::XSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);
  double posX_pos = posX_orig + 2*(Scone_NS::XSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);
  double posY_orig = 0;
  double posY_neg = posY_orig - 2*(Scone_NS::YSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);
  double posY_pos = posY_orig + 2*(Scone_NS::YSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);

  double posX[16] = {posX_orig, posX_orig, posX_neg, posX_neg, posX_neg, posX_pos, posX_pos, posX_pos, posX_orig, posX_orig, posX_neg, posX_neg, posX_neg, posX_pos, posX_pos, posX_pos};
  double posY[16] = {posY_neg, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_orig, posY_pos};
  double posZ[16] = {-300,-300,-300,-300,-300,-300,-300,-300,300,300,300,300,300,300,300,300};

  G4ThreeVector Det_pos;

  for(int i=0; i<m_NumberOfInnerDetector; i++){
  //for(int i=0; i<8; i++){
    m_Assembly++;
    G4String Name = "2x2block_" + to_string(m_Assembly);
    Det_pos = G4ThreeVector(posX[i], posY[i], posZ[i]);

    new G4PVPlacement(G4Transform3D(*Rot, Det_pos), 
        Build2x2Assembly(m_Assembly), 
        Name, world, 
        false, m_Assembly);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::BuildRing1(G4LogicalVolume* world){
  G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
  double posX_orig = 0;
  double posX_neg = posX_orig - 6*(Scone_NS::XSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);
  double posX_pos = posX_orig + 6*(Scone_NS::XSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);
  double posY_orig = 0;
  double posY_neg = posY_orig - 6*(Scone_NS::YSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);
  double posY_pos = posY_orig + 6*(Scone_NS::YSection + Scone_NS::DistanceBetweenBar + 2*Scone_NS::GdThickness);

  double posX[8] = {posX_orig, posX_orig, posX_neg, posX_neg, posX_neg, posX_pos, posX_pos, posX_pos};
  double posY[8] = {posY_neg, posY_pos, posY_neg, posY_orig, posY_pos, posY_neg, posY_orig, posY_pos};
  double posZ[8] = {0,0,0,0,0,0,0,0};

  G4ThreeVector Det_pos;

  for(int i=0; i<m_NumberOfRing1Detector; i++)
  {
    m_Assembly++;
    G4String Name = "6x6block_Ring1_" + to_string(m_Assembly);
    Det_pos = G4ThreeVector(posX[i], posY[i], posZ[i]);
    
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos),
        Build6x6Assembly(m_Assembly, Scone_NS::LengthR1),
        Name, world,
        false, m_Assembly);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::BuildRing2(G4LogicalVolume* world){
  G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
  double shiftX = Scone_NS::XSection + 2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar;
  double shiftY = Scone_NS::YSection + 2*Scone_NS::GdThickness + Scone_NS::DistanceBetweenBar;
  double posX_orig = 0;
  double posX_left = posX_orig - 12*shiftX;
  double posX_right = posX_orig + 12*shiftX;
  double posY_orig = 0;
  double posY_down = posY_orig - 12*shiftY;
  double posY_up = posY_orig + 12*shiftY;

  double posX[16];
  posX[0]= posX_left;
  posX[1]= posX_left;
  posX[2]= posX_left;
  posX[3]= posX_left;
  posX[4]= posX_left;
  posX[5]= posX_right;
  posX[6]= posX_right;
  posX[7]= posX_right;
  posX[8]= posX_right;
  posX[9]= posX_right;
  posX[10]= posX_left+6*shiftX;
  posX[11]= posX_left+12*shiftX;
  posX[12]= posX_left+18*shiftX;
  posX[13]= posX_left+6*shiftX;
  posX[14]= posX_left+12*shiftX;
  posX[15]= posX_left+18*shiftX;

  double posY[16];
  posY[0]= posY_down;
  posY[1]= posY_down+6*shiftY;
  posY[2]= posY_down+12*shiftY;
  posY[3]= posY_down+18*shiftY;
  posY[4]= posY_down+24*shiftY;
  posY[5]= posY_down;
  posY[6]= posY_down+6*shiftY;
  posY[7]= posY_down+12*shiftY;
  posY[8]= posY_down+18*shiftY;
  posY[9]= posY_down+24*shiftY; 
  posY[10]= posY_down;
  posY[11]= posY_down;
  posY[12]= posY_down;
  posY[13]= posY_up;
  posY[14]= posY_up;
  posY[15]= posY_up;

  G4ThreeVector Det_pos;
  G4String Name;
  for(int i=0; i<m_NumberOfRing2Detector; i++)
  {
    m_Assembly++;
    Name = "6x6block_Ring2_" + to_string(m_Assembly);
    Det_pos = G4ThreeVector(posX[i], posY[i], 0);
  
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos),
        Build6x6Assembly(m_Assembly, Scone_NS::LengthR2),
        Name, world,
        false, m_Assembly);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Scone::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Scone")){
    pTree->Branch("Scone", "TSconeData", &m_Event) ;
  }
  pTree->SetBranchAddress("Scone", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Scone::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_SconeScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Scone_NS::ResoEnergy);
    if(Energy>Scone_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Scone_NS::ResoTime);
      int DetectorNbr = level[0];
      int PlasticNbr = level[1];
      m_Event->SetEnergy(DetectorNbr,PlasticNbr,Energy);
      m_Event->SetTime(DetectorNbr,PlasticNbr,Time); 
    }
  }
  Scorer->clear();

  ///////////
  // Process scorer for plastic bar
  ProcessScorers::PS_Process* Process_scorer = (ProcessScorers::PS_Process*) m_SconeScorer->GetPrimitive(2);
  unsigned int ProcessMult = Process_scorer->GetMult();
  bool kPlasticCapture = false;
  double PlasticCaptureTime = 0;
  for(unsigned int i=0; i<ProcessMult; i++){
    string process_name = Process_scorer->GetProcessName(i);
    if(process_name=="nCapture"){
      kPlasticCapture = true;
      PlasticCaptureTime = Process_scorer->GetProcessTime(i);
    }
  }
  vector<double> gamma_energy;
  gamma_energy = Process_scorer->GetGammaEnergy();
  for(unsigned int i=0; i< gamma_energy.size(); i++){
    m_Event->SetGammaEnergy(gamma_energy[i]);
  }
  vector<double> proton_energy;
  vector<double> proton_time;
  proton_energy = Process_scorer->GetProtonEnergy();
  proton_time = Process_scorer->GetProtonTime();
  for(unsigned int i=0; i<proton_energy.size(); i++){
    m_Event->SetProtonEnergy(proton_energy[i]);
    m_Event->SetProtonTime(proton_time[i]);
  }

  //Process_scorer->clear();

  ///////////
  // Process scorer for Gd
  ProcessScorers::PS_Process* GdProcess_scorer = (ProcessScorers::PS_Process*) m_GdScorer->GetPrimitive(0);
  ProcessMult = GdProcess_scorer->GetMult();
  bool kGdCapture = false;
  double GdCaptureTime = 0;
  for(unsigned int i=0; i<ProcessMult; i++){
    string process_name = GdProcess_scorer->GetProcessName(i);
    if(process_name=="nCapture"){
      kGdCapture = true;
      GdCaptureTime = GdProcess_scorer->GetProcessTime(i);
    }
  }
  if(kGdCapture){
    m_Event->SetCapture(2);
    m_Event->SetCaptureTime(GdCaptureTime);
  }
 
  else if(kPlasticCapture){
    m_Event->SetCapture(1);
    m_Event->SetCaptureTime(PlasticCaptureTime);
  }
  
  else m_Event->SetCapture(0);
  GdProcess_scorer->clear();

  ///////////
  // Process scorer for fission chamber
  if(m_BuildFissionChamber==1){
    ProcessScorers::PS_Process* FCProcess_scorer = (ProcessScorers::PS_Process*) m_FCScorer->GetPrimitive(0);
    vector<int> FC_process = FCProcess_scorer->GetFCProcess();
    for(unsigned int i=0; i<FC_process.size(); i++){
      m_Event->SetFCProcess(FC_process[i]);
    }
    FCProcess_scorer->clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Scone::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_SconeScorer = CheckScorer("SconeScorer",already_exist);
  m_GdScorer = CheckScorer("GdScorer",already_exist);
  m_FCScorer = CheckScorer("FCScorer",already_exist);

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; 
  level.push_back(1);
  level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord,0) ;
  G4VPrimitiveScorer* Process= new ProcessScorers::PS_Process("Process",0) ;
  //and register it to the multifunctionnal detector
  m_SconeScorer->RegisterPrimitive(Calorimeter);
  m_SconeScorer->RegisterPrimitive(Interaction);
  m_SconeScorer->RegisterPrimitive(Process);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SconeScorer) ;

  G4VPrimitiveScorer* GdProcess= new ProcessScorers::PS_Process("GdProcess",0) ;
  m_GdScorer->RegisterPrimitive(GdProcess);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_GdScorer) ;

  G4VPrimitiveScorer* FCProcess= new ProcessScorers::PS_Process("FCProcess",0) ;
  m_FCScorer->RegisterPrimitive(FCProcess);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_FCScorer) ;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Scone::InitializeMaterial(){
  m_MaterialVaccuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Scone::Construct(){
  return  (NPS::VDetector*) new Scone();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Scone{
    public:
      proxy_nps_Scone(){
        NPS::DetectorFactory::getInstance()->AddToken("Scone","Scone");
        NPS::DetectorFactory::getInstance()->AddDetector("Scone",Scone::Construct);
      }
  };

  proxy_nps_Scone p_nps_Scone;
}
