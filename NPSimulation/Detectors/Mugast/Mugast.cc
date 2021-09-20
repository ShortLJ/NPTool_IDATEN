/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2019                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Mugast simulation                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Trap.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
// NPTool header
#include "Mugast.hh"
#include "MugastReverseMap.h"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
#include "NPCore.h"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Mugast_NS{
  // Resolution
  const G4double SigmaTime    = 0.212765957 ;// = 500ps
  const G4double SigmaEnergy  = 0.019      ;// 
//  const G4double TimeOffset   = 500         ;// 500 ns stop

  // Threshold
  const G4double EnergyThreshold = 1 * MeV;

  // Geometry
  //const G4double AluStripThickness = 0.4*micrometer ;
  const G4double SiliconThickness  = 300*micrometer ;

  // Square

  // const G4double SquareBase          = 88*mm;
  //  const G4double SquareHeight        = 87*mm;

  const G4double SquareBase          = 91.716*mm;
  const G4double SquareHeight        = 94.916*mm;
  // const G4double SquareHeight        = 194.916*mm;
  const G4double SquareLength        = 1*cm;
  // Trapezoid 
  const G4double TrapezoidBaseLarge     = 91.48*mm; 
  const G4double TrapezoidBaseSmall     = 26*mm; 
  const G4double TrapezoidHeight        = 104.688*mm;
  const G4double TrapezoidLength        = 1*cm;
  //Annular 
   const G4int NbrRingStrips  = 16;
   const G4int NbrSectorStrips = 16;
   const G4int NbQuadrant      = 4;
   const G4double WaferOutterRadius = 50*mm;
   const G4double WaferInnerRadius  = 23*mm;
   const G4double WaferThickness    = 500*micrometer;
   const G4double WaferRCut         = 45.5*mm; 
   const G4double ActiveWaferOutterRadius = 48*mm;
   const G4double ActiveWaferInnerRadius  = 24*mm;
   const G4double ActiveWaferRCut     = 44.5*mm; 
   const G4double PCBPointsX[8]={-40,40,60,60,40,-40,-60,-60};
   const G4double PCBPointsY[8]={60,60,40,-40,-60,-60,-40,40};
   const G4double PCBThickness=3.2*mm;
   //const G4double PCBInnerRadius=0*mm;



}
using namespace Mugast_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Mugast Specific Method
Mugast::Mugast(){
  m_Event = new TMugastData() ;
  m_SquareScorer= 0;
  m_TrapezoidScorer= 0;
  m_AnnularScorer= 0;
  m_SquareDetector = 0;
  m_TrapezoidDetector = 0;
  m_AnnularDetector = 0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));   
  m_VisCylinder = new G4VisAttributes(G4Colour(0, 0, 1, 0.5));   

}

Mugast::~Mugast(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Mugast::AddDetector(int DetectorNumber,string Shape,G4ThreeVector PX1_Y1 ,G4ThreeVector PX1_Y128 ,G4ThreeVector PX128_Y1,G4ThreeVector PX128_Y128){
  m_X1_Y1.push_back(PX1_Y1); // Top Left Corner Position Vector
  m_X1_Y128.push_back(PX1_Y128); // Bottom Left Corner Position Vector
  m_X128_Y1.push_back(PX128_Y1); // Bottom Right Corner Position Vector
  m_X128_Y128.push_back(PX128_Y128); // Center Corner Position Vector
  m_DetectorNumber.push_back(DetectorNumber);
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Mugast::BuildSquareDetector(){
  if(!m_SquareDetector){
    G4String Name = "MugastSquare";

    G4Box*           solidSquare = new G4Box(Name, 0.5*SquareBase, 0.5*SquareHeight, 0.5*SquareLength);
    G4LogicalVolume* logicSquare = new G4LogicalVolume(solidSquare,
        MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"), 
        Name, 0, 0, 0);

    G4VisAttributes* SquareVisAtt = new G4VisAttributes(G4Colour(0.90, 0.90, 0.90)); 
    SquareVisAtt->SetForceWireframe(true); 
    logicSquare->SetVisAttributes(SquareVisAtt);

    // Silicon detector itself
    G4ThreeVector  positionFirstStage = G4ThreeVector(0, 0, 0);

    G4Box*           solidFirstStage = new G4Box("solidFirstStage", 0.5*SquareBase, 0.5*SquareHeight, 0.5*SiliconThickness);
    G4LogicalVolume* logicFirstStage = new G4LogicalVolume(solidFirstStage, 
        MaterialManager::getInstance()->GetMaterialFromLibrary("Si"), 
        "logicFirstStage", 
        0, 0, 0);

    new G4PVPlacement(0, positionFirstStage, logicFirstStage, Name + "_FirstStage", logicSquare, false, 0);
    m_SquareDetector = logicSquare;
    // Set First Stage sensible
    logicFirstStage->SetSensitiveDetector(m_SquareScorer);

    ///Visualisation of FirstStage Strip
    G4VisAttributes* FirstStageVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));   
    logicFirstStage->SetVisAttributes(FirstStageVisAtt);
  }

  return m_SquareDetector;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Mugast::BuildTrapezoidDetector(){
  if(!m_TrapezoidDetector){
    G4String Name = "MugastTrapezoid";

    // Definition of the volume containing the sensitive detector
    G4Trap* solidTrapezoid = new G4Trap(Name, 
        TrapezoidLength*0.5, 0*deg, 0*deg, 
        TrapezoidHeight*0.5,  TrapezoidBaseSmall*0.5,TrapezoidBaseLarge*0.5, 0*deg, 
        TrapezoidHeight*0.5,  TrapezoidBaseSmall*0.5,TrapezoidBaseLarge*0.5, 0*deg);
    G4LogicalVolume* logicTrapezoid = new G4LogicalVolume(solidTrapezoid, 
        MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"), 
        Name, 
        0, 0, 0);

    G4VisAttributes* TrapezoideVisAtt = new G4VisAttributes(G4Colour(0.90, 0.90, 0.90));
    TrapezoideVisAtt->SetForceWireframe(true); 
    logicTrapezoid->SetVisAttributes(TrapezoideVisAtt);

    // Silicon detector itself
    G4ThreeVector  positionFirstStage = G4ThreeVector(0, 0, 0);

    G4Trap* solidFirstStage = new G4Trap("solidFirstStage", 
        0.5*SiliconThickness, 0*deg, 0*deg, 
        TrapezoidHeight*0.5,  TrapezoidBaseSmall*0.5,TrapezoidBaseLarge*0.5, 0*deg, 
        TrapezoidHeight*0.5,  TrapezoidBaseSmall*0.5,TrapezoidBaseLarge*0.5, 0*deg);
    G4LogicalVolume* logicFirstStage = new G4LogicalVolume(solidFirstStage, 
        MaterialManager::getInstance()->GetMaterialFromLibrary("Si"),
        "logicFirstStage", 
        0, 0, 0);

    new G4PVPlacement(0,
        positionFirstStage,
        logicFirstStage,
        Name + "_FirstStage",
        logicTrapezoid,
        false,
        0);

    // Set First Stage sensible
    logicFirstStage->SetSensitiveDetector(m_TrapezoidScorer);

    ///Visualisation of FirstStage Strip
    G4VisAttributes* FirstStageVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));  
    logicFirstStage->SetVisAttributes(FirstStageVisAtt);

    m_TrapezoidDetector=logicTrapezoid;
  }

  return m_TrapezoidDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Mugast::BuildAnnularDetector(){

  if(!m_AnnularDetector){
    G4Material* Silicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    G4Material* Vacuum  = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    G4Material* PCB     = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
    ////////////////////////////////////////////////////////////////
    ////////////// Starting Volume Definition //////////////////////
    ////////////////////////////////////////////////////////////////
    // Name of the module
    G4String Name = "MugastAnnular";

    // Building the PCB
    // The PCB is a simple extruded volume from 8 reference points
    vector<G4TwoVector> polygon;
    for(unsigned int i = 0 ; i < 8 ; i++){
      G4TwoVector Point(PCBPointsX[i],PCBPointsY[i]);
      polygon.push_back(Point);
    }

    // Master volume containing all the detector
    G4ExtrudedSolid* solidAnnularS1 = new G4ExtrudedSolid(Name,
        polygon,
        PCBThickness*0.5,
        G4TwoVector(0,0),1,
        G4TwoVector(0,0),1);

    // Definition of the volume containing the sensitive detector
    m_AnnularDetector= new G4LogicalVolume(solidAnnularS1, Vacuum, Name, 0, 0, 0);
    m_AnnularDetector->SetVisAttributes(G4VisAttributes::Invisible);

    // PCB Base
    G4ExtrudedSolid* solidPCBBase = new G4ExtrudedSolid("PCBBase",
        polygon,
        PCBThickness*0.5,
        G4TwoVector(0,0),1,
        G4TwoVector(0,0),1);   

    // Wafer Shape to be substracted to the PCB
    G4Tubs* solidWaferShapeBase = new G4Tubs("WaferShape", 
        0,
        WaferOutterRadius,
        PCBThickness,
        0*deg, 
        360*deg); 

    // A no rotation matrix is always handy ;)
    G4RotationMatrix* norotation = new  G4RotationMatrix(); 
    // Rotation of the box that make the Si cut                       
    G4RotationMatrix* cutrotation = new  G4RotationMatrix(G4ThreeVector(0,0,1),45*deg);                        
    G4ThreeVector cutposition1(80*mm+WaferRCut,0,0); cutposition1.setPhi(45*deg);
    G4Transform3D transform1(*cutrotation,cutposition1);

    G4Box* solidCutout = new G4Box("cuttout",80*mm,80*mm,80*mm); 

    G4SubtractionSolid* solidWaferShape1 = new G4SubtractionSolid("WaferShape1",
        solidWaferShapeBase,
        solidCutout,
        transform1);


    G4ThreeVector cutposition2(-80*mm-WaferRCut,0,0); cutposition2.setPhi(-135*deg);
    G4Transform3D transform2(*cutrotation,cutposition2);
    G4SubtractionSolid* solidWaferShape = new G4SubtractionSolid("WaferShape",
        solidWaferShape1,
        solidCutout,
        transform2);


    // PCB final
    G4SubtractionSolid* solidPCB = new G4SubtractionSolid("MugastAnnular_PCB1",
        solidPCBBase,
        solidWaferShape);

    G4LogicalVolume* logicPCB = new G4LogicalVolume(solidPCB, PCB, "MugastAnnular_PCB", 0, 0, 0);

    new G4PVPlacement(G4Transform3D(*norotation, G4ThreeVector()),
        logicPCB,
        "MugastAnnular_PCB",
        m_AnnularDetector,
        false,
        0);

    G4VisAttributes* PCBVisAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2)) ;
    logicPCB->SetVisAttributes(PCBVisAtt);


    // Wafer itself
    G4Tubs* solidWaferBase = new G4Tubs("Wafer", 
        WaferInnerRadius,
        WaferOutterRadius,
        0.5*WaferThickness,
        0*deg, 
        360*deg); 

    G4SubtractionSolid* solidWafer1 = new G4SubtractionSolid("Wafer1",
        solidWaferBase,
        solidCutout,
        transform1);

    G4SubtractionSolid* solidWafer = new G4SubtractionSolid("Wafer",
        solidWafer1,
        solidCutout,
        transform2);

    G4LogicalVolume* logicWafer = new G4LogicalVolume(solidWafer, Silicon, "MugastAnnular_Wafer", 0, 0, 0);
    new G4PVPlacement(G4Transform3D(*norotation, G4ThreeVector()),
        logicWafer,
        "MugastAnnular_Wafer",
        m_AnnularDetector,
        false,
        0);

    G4VisAttributes* SiVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3)) ;
    logicWafer->SetVisAttributes(SiVisAtt); 

    // Active Wafer
    G4Tubs* solidActiveWaferBase = new G4Tubs("ActiveWafer", 
        ActiveWaferInnerRadius,
        ActiveWaferOutterRadius,
        0.5*WaferThickness,
        0*deg, 
        360*deg); 


    G4ThreeVector activecutposition1(80*mm+ActiveWaferRCut,0,0); activecutposition1.setPhi(45*deg);
    G4Transform3D activetransform1(*cutrotation,activecutposition1);

    G4SubtractionSolid* solidActiveWafer1 = new G4SubtractionSolid("ActiveWafer1",
        solidActiveWaferBase,
        solidCutout,
        activetransform1);

    G4ThreeVector activecutposition2(-80*mm-ActiveWaferRCut,0,0); activecutposition2.setPhi(-135*deg);
    G4Transform3D activetransform2(*cutrotation,activecutposition2);

    G4SubtractionSolid* solidActiveWafer = new G4SubtractionSolid("ActiveWafer",
        solidActiveWafer1,
        solidCutout,
        activetransform2);

    G4LogicalVolume* logicActiveWafer = new G4LogicalVolume(solidActiveWafer, Silicon, "MugastAnnular_ActiveWafer", 0, 0, 0);
    new G4PVPlacement(G4Transform3D(*norotation, G4ThreeVector()),
        logicActiveWafer,
        "MugastAnnular_ActiveWafer",
        logicWafer,
        false,
        0);

    logicActiveWafer->SetVisAttributes(SiVisAtt);

    // Set Silicon strip sensible
    logicActiveWafer->SetSensitiveDetector(m_AnnularScorer);
  }
  return m_AnnularDetector;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Mugast::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Mugast");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector found" << endl;

  // Cartesian Case
  vector<string> cart = {"DetectorNumber","X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128"};
  vector<string> annular = {"DetectorNumber","Center"};

  for (unsigned int i = 0; i < blocks.size(); i++) {

    string shape = blocks[i]->GetMainValue();

    if (NPOptionManager::getInstance()->GetVerboseLevel())
      cout << endl << "//// Mugast detector " << shape << endl;

    if (blocks[i]->HasTokenList(cart)&& (shape=="Square"|| shape=="Trapezoid")) {
      int DetectorNumber = blocks[i]->GetInt("DetectorNumber");
      G4ThreeVector A
        = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1", "mm"));
      G4ThreeVector B
        = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y128", "mm"));
      G4ThreeVector C
        = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y1", "mm"));
      G4ThreeVector D
        = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y128", "mm"));

      AddDetector(DetectorNumber,shape,A,B,C,D);
    }
    else if (blocks[i]->HasTokenList(annular)&& (shape=="Annular")) {
      int DetectorNumber = blocks[i]->GetInt("DetectorNumber");
      G4ThreeVector A
        = NPS::ConvertVector(blocks[i]->GetTVector3("Center", "mm"));
      AddDetector(DetectorNumber,shape,A);
    }

    else {
      cout << "WARNING: Missing token for Mugast blocks, check your input "
        "file"
        << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After Detecor
// onstruction::AddDetector Method
void Mugast::ConstructDetector(G4LogicalVolume* world){

  for (unsigned short i = 0 ; i < m_DetectorNumber.size() ; i++) {
    if(m_Shape[i]=="Square"){
      G4RotationMatrix* rot    = NULL                   ;
      G4ThreeVector     pos    = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     u      = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     v      = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     w      = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     Center = G4ThreeVector(0, 0, 0) ;
      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      u = m_X128_Y1[i] - m_X1_Y1[i];
      u = u.unit();

      v = (m_X1_Y128[i] + m_X128_Y128[i] - m_X1_Y1[i] - m_X128_Y1[i]);
      v = v.unit();

      w = u.cross(v);
      w = w.unit();

      Center = (m_X1_Y1[i] + m_X1_Y128[i] + m_X128_Y1[i] + m_X128_Y128[i]) / 4;

      // Passage Matrix from Lab Referential to Telescope Referential
      rot = new G4RotationMatrix(u, v, w);
      // translation to place Telescope
      pos = w * SiliconThickness* 0.5 + Center;
      new G4PVPlacement(G4Transform3D(*rot, pos), BuildSquareDetector(), "MugastSquare", world, false, m_DetectorNumber[i]);
    }
    else if(m_Shape[i]=="Trapezoid"){
      G4RotationMatrix* rot    = NULL                   ;
      G4ThreeVector     pos    = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     u      = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     v      = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     w      = G4ThreeVector(0, 0, 0) ;
      G4ThreeVector     Center = G4ThreeVector(0, 0, 0) ;
      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      u = m_X128_Y1[i] - m_X1_Y1[i];
      u = u.unit();

      v = (m_X1_Y128[i] + m_X128_Y128[i] - m_X1_Y1[i] - m_X128_Y1[i] );
      v = v.unit();

      w = u.cross(v);
      w = w.unit();

      Center = (m_X1_Y1[i] + m_X1_Y128[i] + m_X128_Y1[i] + m_X128_Y128[i]) / 4;

      // Passage Matrix from Lab Referential to Telescope Referential
      rot = new G4RotationMatrix(u, v, w);
      rot->rotate(180*deg,w);
      // translation to place Telescope
      pos = w * SiliconThickness* 0.5 + Center;
      new G4PVPlacement(G4Transform3D(*rot, pos), BuildTrapezoidDetector(), "MugastTrapezoid", world, false, m_DetectorNumber[i]);
    }

    else if(m_Shape[i]=="Annular"){
      G4RotationMatrix* rot = new G4RotationMatrix();
      new G4PVPlacement(G4Transform3D(*rot,m_X1_Y1[i]), BuildAnnularDetector(), "MugastAnnular", world, false, m_DetectorNumber[i]);
      }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Mugast::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Mugast")){
    pTree->Branch("Mugast", "TMugastData", &m_Event) ;
  }
  pTree->SetBranchAddress("Mugast", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Mugast::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Square
  DSSDScorers::PS_Rectangle* SquareScorer = (DSSDScorers::PS_Rectangle*) m_SquareScorer->GetPrimitive(0);


  // Loop on the Square map
  unsigned int sizeFront= SquareScorer->GetLengthMult();

  for (unsigned int i=0 ; i<sizeFront ; i++){

    double Energy = RandGauss::shoot(SquareScorer->GetEnergyLength(i), SigmaEnergy);

    if(Energy>EnergyThreshold){
      double Time       = RandGauss::shoot(SquareScorer->GetTimeLength(i), SigmaTime);
      int DetNbr        = SquareScorer->GetDetectorLength(i);
      int StripFront    = SquareScorer->GetStripLength(i);

      m_Event->SetDSSDXE(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseSquareX[StripFront-1],
          NPL::EnergyToADC(Energy, 0, 63, 8192, 16384));

      m_Event->SetDSSDXT(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseSquareX[StripFront-1],
          NPL::EnergyToADC(Time, 0, 1000, 8192, 16384));

    }
  } 

  unsigned int sizeBack= SquareScorer->GetWidthMult();
  for (unsigned int i=0 ; i<sizeBack ; i++){
    double Energy = RandGauss::shoot(SquareScorer->GetEnergyWidth(i), SigmaEnergy);

    if(Energy>EnergyThreshold){
      double Time       = RandGauss::shoot(SquareScorer->GetTimeWidth(i), SigmaTime);
      int DetNbr        = SquareScorer->GetDetectorWidth(i);
      int StripBack     = SquareScorer->GetStripWidth(i);

      m_Event->SetDSSDYE(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseSquareY[StripBack-1],
          NPL::EnergyToADC(Energy, 0, 63, 8192, 0));

      m_Event->SetDSSDYT(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseSquareY[StripBack-1],
          NPL::EnergyToADC(Time, 0, 1000, 8192, 16384));
    }
  }
  // clear map for next event
  SquareScorer->clear();

  ///////////
  // Trapezoid
  DSSDScorers::PS_Rectangle* TrapezoidScorer = (DSSDScorers::PS_Rectangle*) m_TrapezoidScorer->GetPrimitive(0);


  // Loop on the Trapezoid map
  sizeFront= TrapezoidScorer->GetLengthMult();

  for (unsigned int i=0 ; i<sizeFront ; i++){

    double Energy = RandGauss::shoot(TrapezoidScorer->GetEnergyLength(i), SigmaEnergy);

    if(Energy>EnergyThreshold){
      double Time       = RandGauss::shoot(TrapezoidScorer->GetTimeLength(i), SigmaTime);
      int DetNbr        = TrapezoidScorer->GetDetectorLength(i);
      int StripFront    = TrapezoidScorer->GetStripLength(i);

      m_Event->SetDSSDXE(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseTrapezeX[StripFront-1],
          NPL::EnergyToADC(Energy, 0, 63, 8192, 16384));

      m_Event->SetDSSDXT(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseTrapezeX[StripFront-1],
          NPL::EnergyToADC(Time, 0, 1000, 8192, 16384));

    }
  } 

  sizeBack= TrapezoidScorer->GetWidthMult();
  for (unsigned int i=0 ; i<sizeBack ; i++){

    double Energy = RandGauss::shoot(TrapezoidScorer->GetEnergyWidth(i), SigmaEnergy);

    if(Energy>EnergyThreshold){
      double Time       = RandGauss::shoot(TrapezoidScorer->GetTimeWidth(i), SigmaTime);
      int DetNbr        = TrapezoidScorer->GetDetectorWidth(i);
      int StripBack     = 128-TrapezoidScorer->GetStripWidth(i)+1;

      m_Event->SetDSSDYE(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseTrapezeY[StripBack-1],
          NPL::EnergyToADC(Energy, 0, 63, 8192, 0));

      m_Event->SetDSSDYT(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseTrapezeY[StripBack-1],
          NPL::EnergyToADC(Time, 0, 1000, 8192, 16384));
    }
  }
  // clear map for next event
  TrapezoidScorer->clear();

  ///////////
  // Annular
  DSSDScorers::PS_Annular* AnnularScorer = (DSSDScorers::PS_Annular*) m_AnnularScorer->GetPrimitive(0);


  // Loop on the Annular map
  sizeFront= AnnularScorer->GetRingMult();
  unsigned int sizeQuadrant= AnnularScorer->GetQuadrantMult();

  for (unsigned int i=0 ; i<sizeFront ; i++){

    double Energy = RandGauss::shoot(AnnularScorer->GetEnergyRing(i), SigmaEnergy);

    if(Energy>EnergyThreshold){
      double Time       = RandGauss::shoot(AnnularScorer->GetTimeRing(i), SigmaTime);
      unsigned int DetNbr        = AnnularScorer->GetDetectorRing(i);
      unsigned int StripFront    = AnnularScorer->GetStripRing(i);
   
      // Check for associated Quadrant strip
      int StripQuadrant = 0;
      for(unsigned int q = 0 ; q < sizeQuadrant ; q++){
        if(AnnularScorer->GetDetectorQuadrant(q)==DetNbr){
          StripQuadrant = AnnularScorer->GetStripQuadrant(q)-1;
          break;
          }
      }
      StripFront=StripFront+StripQuadrant*NbrRingStrips;
      m_Event->SetDSSDXE(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseAnnularX[StripFront-1],
          NPL::EnergyToADC(Energy, 0, 63, 8192, 16384));

      m_Event->SetDSSDXT(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseAnnularX[StripFront-1],
          NPL::EnergyToADC(Time, 0, 1000, 8192, 16384));

    }
  } 

  sizeBack= AnnularScorer->GetSectorMult();
  for (unsigned int i=0 ; i<sizeBack ; i++){

    double Energy = RandGauss::shoot(AnnularScorer->GetEnergySector(i), SigmaEnergy);

    if(Energy>EnergyThreshold){
      double Time       = RandGauss::shoot(AnnularScorer->GetTimeSector(i), SigmaTime);
      int DetNbr        = AnnularScorer->GetDetectorSector(i);
      int StripBack     = AnnularScorer->GetStripSector(i);
      m_Event->SetDSSDYE(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseAnnularY[StripBack-1],
          NPL::EnergyToADC(Energy, 0, 63, 8192, 0));

      m_Event->SetDSSDYT(MG_DetectorType::MG_NOCHANGE,DetNbr,
          MUGAST_MAP::ReverseAnnularY[StripBack-1],
          NPL::EnergyToADC(Time, 0, 1000, 8192, 16384));
    }
  }
  // clear map for next event
  AnnularScorer->clear();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Mugast::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_SquareScorer= CheckScorer("SquareScorer",already_exist) ;
  m_TrapezoidScorer= CheckScorer("TrapezoidScorer",already_exist) ;
  m_AnnularScorer= CheckScorer("AnnularScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  G4VPrimitiveScorer* SquareScorer     = new DSSDScorers::PS_Rectangle("MugastSquareScorer",1,
      SquareBase,
      SquareHeight,
      128,
      128);


  G4VPrimitiveScorer* TrapezoidScorer  = new DSSDScorers::PS_Rectangle("MugastTrapezoidScorer",1,
      TrapezoidBaseLarge,
      TrapezoidHeight,
      128,
      128);

  G4VPrimitiveScorer* AnnularScorer = new  DSSDScorers::PS_Annular("MugastAnnularScorer",
        2,
        ActiveWaferInnerRadius,
        ActiveWaferOutterRadius,
        -8*22.5*deg, //MUST2 campaign 2009, See: Phd Sandra Giron
        +8*22.5*deg,
        NbrRingStrips,
        NbrSectorStrips,
        NbQuadrant);


  G4VPrimitiveScorer* InteractionS= new InteractionScorers::PS_Interactions("InteractionS",ms_InterCoord, 0) ;
  G4VPrimitiveScorer* InteractionT= new InteractionScorers::PS_Interactions("InteractionT",ms_InterCoord, 0) ;
  G4VPrimitiveScorer* InteractionA= new InteractionScorers::PS_Interactions("InteractionA",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_SquareScorer->RegisterPrimitive(SquareScorer);
  m_SquareScorer->RegisterPrimitive(InteractionS);
  m_TrapezoidScorer->RegisterPrimitive(TrapezoidScorer);
  m_TrapezoidScorer->RegisterPrimitive(InteractionT);
  m_AnnularScorer->RegisterPrimitive(AnnularScorer);
  m_AnnularScorer->RegisterPrimitive(InteractionA);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_SquareScorer) ;
  G4SDManager::GetSDMpointer()->AddNewDetector(m_TrapezoidScorer) ;
  G4SDManager::GetSDMpointer()->AddNewDetector(m_AnnularScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Mugast::Construct(){
  return  (NPS::VDetector*) new Mugast();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Mugast{
    public:
      proxy_nps_Mugast(){
        NPS::DetectorFactory::getInstance()->AddToken("Mugast","Mugast");
        NPS::DetectorFactory::getInstance()->AddDetector("Mugast",Mugast::Construct);
      }
  };

  proxy_nps_Mugast p_nps_Mugast;
}
