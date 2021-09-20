/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe the Tiara Silicon array                              *
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
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4TwoVector.hh"
//G4 sensitive
#include "G4SDManager.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4PVDivision.hh"
#include "G4SubtractionSolid.hh"

// NPS
#include "Tiara.hh"
#include "NPSDetectorFactory.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "SiliconScorers.hh"

// NPL
#include "NPOptionManager.h"
#include "RootOutput.h"
using namespace TIARA;

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Tiara::Tiara(){
  InitializeMaterial();
  m_EventBarrel = new TTiaraBarrelData();
  m_EventHyball = new TTiaraHyballData();

  // Dark Grey
  SiliconVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3)) ;
  // Green
  PCBVisAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2)) ;
  // Gold Yellow
  PADVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.2)) ;
  // Light Grey
  FrameVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
  // Light Blue
  GuardRingVisAtt = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1)) ;

  m_boolChamber = false;
  m_boolInner = false;
  m_boolOuter = false;
  m_InnerBarrelScorer = 0 ;
  m_OuterBarrelScorer = 0 ;
  m_HyballScorer = 0 ;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Tiara::~Tiara(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Tiara::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("Tiara","Barrel");

  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> token = {"InnerBarrel","OuterBarrel","Chamber"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Barrel " << i+1 <<  endl;
      m_boolInner   = blocks[i]->GetInt("InnerBarrel");
      m_boolOuter   = blocks[i]->GetInt("OuterBarrel");
      m_boolChamber = blocks[i]->GetInt("Chamber");
    }

    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  blocks.clear();
//by Shuya 180426
 // blocks = parser.GetAllBlocksWithToken("TiaraHyballWedge");
  blocks = parser.GetAllBlocksWithTokenAndValue("Tiara","Hyball");

  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  token.clear();
  token = {"Z","R","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Hyball Wedge" << i+1 <<  endl;
      double Z = blocks[i]->GetDouble("Z","mm");
      double R = blocks[i]->GetDouble("R","mm");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      m_HyballZ.push_back(Z);
      m_HyballR.push_back(R);
      m_HyballPhi.push_back(Phi);
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
void Tiara::ConstructDetector(G4LogicalVolume* world){

  if(m_boolChamber)
    ConstructChamber(world);

  if(m_boolInner)
    ConstructInnerBarrel(world);

  if(m_boolOuter)
    ConstructOuterBarrel(world);

  if(m_HyballZ.size())
    ConstructHyball(world);
}
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Tiara::ReadSensitive(const G4Event* event){
  m_EventBarrel->Clear();
  m_EventHyball->Clear();

  // InnerBarrel //
  NPS::HitsMap<G4double*>* InnerBarrelHitMap;
  std::map<G4int, G4double**>::iterator InnerBarrel_itr;
  G4int InnerBarrelCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("Tiara_InnerBarrelScorer/InnerBarrel");
  InnerBarrelHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(InnerBarrelCollectionID));

  // Loop on the InnerBarrel map
  for (InnerBarrel_itr = InnerBarrelHitMap->GetMap()->begin() ; InnerBarrel_itr != InnerBarrelHitMap->GetMap()->end() ; InnerBarrel_itr++){
    G4double* Info = *(InnerBarrel_itr->second);

    // Downstream Energy
  //by Shuya 180504. Note I changed this unit from MeV to keV because Barrel data are given in keV in real data.
    //double ED = RandGauss::shoot(Info[0],ResoEnergyInnerBarrel);
    double ED = RandGauss::shoot(Info[0]/keV,ResoEnergyInnerBarrel);
    if(ED>EnergyThreshold){
      m_EventBarrel->SetFrontDownstreamE(Info[3],Info[4],ED);
      m_EventBarrel->SetFrontDownstreamT(Info[3],Info[4],Info[2]);
    }

    // Upstream Energy
  //by Shuya 180504. Note I changed this unit from MeV to keV because Barrel data are given in keV in real data.
    //double EU = RandGauss::shoot(Info[1],ResoEnergyInnerBarrel);
    double EU = RandGauss::shoot(Info[1]/keV,ResoEnergyInnerBarrel);
    if(EU>EnergyThreshold){
      m_EventBarrel->SetFrontUpstreamE(Info[3],Info[4],EU);
      m_EventBarrel->SetFrontUpstreamT(Info[3],Info[4],Info[2]);
	//cout << "BARREL " << Info[3] << " STRIP " << Info[4] << " " << Info[9]/deg << " " << Info[5] << " " << Info[6] << " " << Info[7] << endl;
    }

    // Back Energy
  //by Shuya 180504. Note I changed this unit from MeV to keV because Barrel data are given in keV in real data.
    //double EB = RandGauss::shoot(Info[1]+Info[0],ResoEnergyInnerBarrel);
    double EB = RandGauss::shoot(Info[1]/keV+Info[0]*keV,ResoEnergyInnerBarrel);
    if(EB>EnergyThreshold){
      m_EventBarrel->SetBackE(Info[3],EB);
      m_EventBarrel->SetBackT(Info[3],Info[2]);
    }
        // Interaction Coordinates
    ms_InterCoord->SetDetectedPositionX(Info[5]) ;
    ms_InterCoord->SetDetectedPositionY(Info[6]) ;
    ms_InterCoord->SetDetectedPositionZ(Info[7]) ;
    ms_InterCoord->SetDetectedAngleTheta(Info[8]/deg) ;
    ms_InterCoord->SetDetectedAnglePhi(Info[9]/deg) ;
  }
  // Clear Map for next event
  InnerBarrelHitMap->clear();

  // OuterBarrel //
  NPS::HitsMap<G4double*>* OuterBarrelHitMap;
  std::map<G4int, G4double**>::iterator OuterBarrel_itr;
  G4int OuterBarrelCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("Tiara_OuterBarrelScorer/OuterBarrel");
  OuterBarrelHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(OuterBarrelCollectionID));

  // Loop on the OuterBarrel map
  for (OuterBarrel_itr = OuterBarrelHitMap->GetMap()->begin() ; OuterBarrel_itr != OuterBarrelHitMap->GetMap()->end() ; OuterBarrel_itr++){
    G4double* Info = *(OuterBarrel_itr->second);

    double E = RandGauss::shoot(Info[0]/keV,ResoEnergyOuterBarrel);
    if(E>EnergyThreshold){
      m_EventBarrel->SetOuterE(Info[7],Info[9],E);
      m_EventBarrel->SetOuterT(Info[7],Info[9],Info[1]);
    }
        // Interaction Coordinates
    ms_InterCoord->SetDetectedPositionX(Info[5]) ;
    ms_InterCoord->SetDetectedPositionY(Info[6]) ;
    ms_InterCoord->SetDetectedPositionZ(Info[7]) ;
    ms_InterCoord->SetDetectedAngleTheta(Info[8]/deg) ;
    ms_InterCoord->SetDetectedAnglePhi(Info[9]/deg) ;
  }
  // Clear Map for next event
  OuterBarrelHitMap->clear();


  // Hyball //
  NPS::HitsMap<G4double*>* HyballHitMap;
  std::map<G4int, G4double**>::iterator Hyball_itr;
  G4int HyballCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("Tiara_HyballScorer/Hyball");
  HyballHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(HyballCollectionID));

  // Loop on the Hyball map
  for (Hyball_itr = HyballHitMap->GetMap()->begin() ; Hyball_itr != HyballHitMap->GetMap()->end() ; Hyball_itr++){
    G4double* Info = *(Hyball_itr->second);

    // Front Energy
    double EF = RandGauss::shoot(Info[0],ResoEnergyHyball);
    if(EF>EnergyThreshold){
      int RingNumber=Info[8];
	//by Shuya 171009
      RingNumber=abs(RingNumber);
      //RingNumber=abs(RingNumber-17);
      Info[8]=RingNumber;
      m_EventHyball->SetRingE(Info[7],Info[8],EF);
      m_EventHyball->SetRingT(Info[7],Info[8],Info[1]);
    }

    // Back Energy
//by Shuya 171009. Infor[1] is Timing data... no make sense.
    //double EB = RandGauss::shoot(Info[1]+Info[0],ResoEnergyHyball);
    double EB = RandGauss::shoot(Info[0],ResoEnergyHyball);
    if(EB>EnergyThreshold){
//by Shuya 171012.
	double m_axis = -100.0;
	double m_phi = Info[6];
	if(Info[6]<0)	m_phi = Info[6]+2.0*M_PI;

	if(Info[7]==1)	m_axis = 210.0/180.0*M_PI;
	else if(Info[7]==2)	m_axis = 150.0/180.0*M_PI;
	else if(Info[7]==3)	m_axis = 90.0/180.0*M_PI;
	else if(Info[7]==4)	m_axis = 30.0/180.0*M_PI;
	else if(Info[7]==5)	m_axis = 330.0/180.0*M_PI;
	else if(Info[7]==6)	m_axis = 270.0/180.0*M_PI;

    	double m_StripPitchSector_Tiara = (0.5*HYBALL_ActiveWafer_Angle-(-0.5*HYBALL_ActiveWafer_Angle))/HYBALL_NumberOfRadialStrip;
	Info[9] = (int)((m_phi - (m_axis - 0.5*HYBALL_ActiveWafer_Angle)) / m_StripPitchSector_Tiara ) + 1 ;

//by Shuya 171009
      //m_EventHyball->SetSectorE(Info[7],Info[9],EF);
      m_EventHyball->SetSectorE(Info[7],Info[9],EB);
      m_EventHyball->SetSectorT(Info[7],Info[9],Info[1]);
    }
    // Interaction Coordinates
//by Shuya 171009
/*
    ms_InterCoord->SetDetectedPositionX(Info[5]) ;
    ms_InterCoord->SetDetectedPositionY(Info[6]) ;
    ms_InterCoord->SetDetectedPositionZ(Info[7]) ;
    ms_InterCoord->SetDetectedAngleTheta(Info[8]/deg) ;
    ms_InterCoord->SetDetectedAnglePhi(Info[9]/deg) ;
*/
    ms_InterCoord->SetDetectedPositionX(Info[2]) ;
    ms_InterCoord->SetDetectedPositionY(Info[3]) ;
    ms_InterCoord->SetDetectedPositionZ(Info[4]) ;
    ms_InterCoord->SetDetectedAngleTheta(Info[5]/deg) ;
    ms_InterCoord->SetDetectedAnglePhi(Info[6]/deg) ;
  }
  // Clear Map for next event
  HyballHitMap->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tiara::InitializeScorers(){
  //Look for previous definition of the scorer (geometry reload)
  bool already_exist = false;
  m_InnerBarrelScorer = CheckScorer("Tiara_InnerBarrelScorer",already_exist);
  m_OuterBarrelScorer = CheckScorer("Tiara_OuterBarrelScorer",already_exist);
  m_HyballScorer = CheckScorer("Tiara_HyballScorer",already_exist);

  // if the scorer were created previously nothing else need to be made
  if(already_exist) return;

  G4VPrimitiveScorer* InnerBarrel = new SILICONSCORERS::PS_Silicon_Resistive("InnerBarrel",1,
      INNERBARREL_ActiveWafer_Length,
      INNERBARREL_ActiveWafer_Width,
      INNERBARREL_NumberOfStrip);

  m_InnerBarrelScorer->RegisterPrimitive(InnerBarrel);

  G4VPrimitiveScorer* OuterBarrel = new SILICONSCORERS::PS_Silicon_Rectangle("OuterBarrel",1,
      INNERBARREL_ActiveWafer_Length,
      INNERBARREL_ActiveWafer_Width,
      1,
      OUTERBARREL_NumberOfStrip);

  m_OuterBarrelScorer->RegisterPrimitive(OuterBarrel);

  G4VPrimitiveScorer* Hyball= new SILICONSCORERS::PS_Silicon_Annular("Hyball",1,
      HYBALL_ActiveWafer_InnerRadius,
      HYBALL_ActiveWafer_OuterRadius,
      -0.5*HYBALL_ActiveWafer_Angle,0.5*HYBALL_ActiveWafer_Angle,
      HYBALL_NumberOfAnnularStrip,
      HYBALL_NumberOfRadialStrip);

  m_HyballScorer->RegisterPrimitive(Hyball);

  //   Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_InnerBarrelScorer) ;
  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterBarrelScorer) ;
  G4SDManager::GetSDMpointer()->AddNewDetector(m_HyballScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tiara::InitializeRootOutput(){
  TTree *pTree = RootOutput::getInstance()->GetTree();
  if(!pTree->FindBranch("TiaraBarrel")){
    pTree->Branch("TiaraBarrel", "TTiaraBarrelData", &m_EventBarrel) ;
  }

  if(!pTree->FindBranch("TiaraHyball")){
    pTree->Branch("TiaraHyball", "TTiaraHyballData", &m_EventHyball) ;
  }
  // This insure that the object are correctly bind in case of
  // a redifinition of the geometry in the simulation
  pTree->SetBranchAddress("TiaraBarrel", &m_EventBarrel) ;
  pTree->SetBranchAddress("TiaraHyball", &m_EventHyball) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tiara::ConstructInnerBarrel(G4LogicalVolume* world){
  // Tiara Barrel
  // The Barrel is made of 8 identical resistive strip detector
  // The PCB is made from a G4ExtrudeSolid, because it has beveled edge
  // the pcb is a substracted volume
  // the wafer goes into the hole, but a 15mm part is still covered by some PCB
  // the whole things is design so the local reference is the one of the wafer

  // Start by making a full pcb
  // We start by the definition of the point forming a PCB cross section
  vector<G4TwoVector> PCBCrossSection;
  double l1 = INNERBARREL_PCB_Thickness*0.5/tan(INNERBARREL_PCB_Bevel1_Theta);
  double l2 = INNERBARREL_PCB_Thickness*0.5/tan(INNERBARREL_PCB_Bevel2_Theta);

  PCBCrossSection.push_back(G4TwoVector(INNERBARREL_PCB_Width/2.-l2,-INNERBARREL_PCB_Thickness*0.5));
  PCBCrossSection.push_back(G4TwoVector(INNERBARREL_PCB_Width/2.,0));
  PCBCrossSection.push_back(G4TwoVector(INNERBARREL_PCB_Width/2.-l1,INNERBARREL_PCB_Thickness*0.5));

  PCBCrossSection.push_back(G4TwoVector(-INNERBARREL_PCB_Width/2.+l1,INNERBARREL_PCB_Thickness*0.5));
  PCBCrossSection.push_back(G4TwoVector(-INNERBARREL_PCB_Width/2.,0));
  PCBCrossSection.push_back(G4TwoVector(-INNERBARREL_PCB_Width/2.+l2,-INNERBARREL_PCB_Thickness*0.5));

  G4ExtrudedSolid* PCBFull =
    new G4ExtrudedSolid("PCBFull",
        PCBCrossSection,
        INNERBARREL_PCB_Length/2.,
        G4TwoVector(0,0),1,
        G4TwoVector(0,0),1);

  // A box having Wafer dimension but thicker than the PCB
  // Will be used to remove material from the PCB to have space for the wafer

  // Calculate the hole shift within the PCB
  G4ThreeVector HoleShift = G4ThreeVector(
      0,
      0,
      INNERBARREL_PCB_Offset-(INNERBARREL_PCB_Length/2-INNERBARREL_PCB_HoleLength/2));

  G4Box*  HoleShape = new G4Box("HoleShape",
      INNERBARREL_ActiveWafer_Width/2.,
      INNERBARREL_PCB_Thickness/2.+0.1*mm,
      INNERBARREL_PCB_HoleLength/2.);

  G4Box*  WaferShape = new G4Box("WaferShape",
      INNERBARREL_InertWafer_Width/2.,
      INNERBARREL_PCB_Thickness/2.,
      INNERBARREL_InertWafer_Length/2.);

  // The Silicon Wafer itself
  G4Box*  InertWaferFull = new G4Box("InertWaferFull",
      INNERBARREL_InertWafer_Width/2.,
      INNERBARREL_InertWafer_Thickness/2.,
      INNERBARREL_InertWafer_Length/2.);

  G4Box*  ActiveWafer = new G4Box("ActiveWafer",
      INNERBARREL_ActiveWafer_Width/2.,
      INNERBARREL_ActiveWafer_Thickness/2.,
      INNERBARREL_ActiveWafer_Length/2.);

  G4Box*  DeadLayer = new G4Box("DeadLayer",
      INNERBARREL_ActiveWafer_Width/2.,
      INNERBARREL_ActiveWafer_DeadLayerThickness/2.,
      INNERBARREL_ActiveWafer_Length/2.);

  G4Box*  ActiveWaferShape = new G4Box("ActiveWaferShape",
      INNERBARREL_ActiveWafer_Width/2.,
      INNERBARREL_PCB_Thickness/2.,
      INNERBARREL_ActiveWafer_Length/2.);


  // Substracting the hole Shape from the Stock PCB
  G4SubtractionSolid* PCB_1 = new G4SubtractionSolid("PCB_1", PCBFull, HoleShape,
      new G4RotationMatrix,HoleShift);

  // Substracting the wafer space from the Stock PCB
  G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_1, WaferShape,
      new G4RotationMatrix,
      G4ThreeVector(0,INNERBARREL_PCB_Thickness/2.-INNERBARREL_PCB_WaferDepth,0));

  // Substract active part from inert part of the Wafer
  G4SubtractionSolid* InertWafer = new G4SubtractionSolid("InertWafer", InertWaferFull, ActiveWaferShape,
      new G4RotationMatrix,
      G4ThreeVector(0,0,0));

  // Master Volume that encompass everything else
  G4LogicalVolume* logicBarrelDetector =
    new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector", 0, 0, 0);
  logicBarrelDetector->SetVisAttributes(G4VisAttributes::Invisible);

  // Sub Volume PCB
  G4LogicalVolume* logicPCB =
    new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
  logicPCB->SetVisAttributes(PCBVisAtt);

  // Sub Volume Wafer
  G4LogicalVolume* logicInertWafer =
    new G4LogicalVolume(InertWafer,m_MaterialSilicon,"logicInertWafer", 0, 0, 0);
  logicInertWafer->SetVisAttributes(GuardRingVisAtt);

  G4LogicalVolume* logicActiveWafer =
    new G4LogicalVolume(ActiveWafer,m_MaterialSilicon,"logicActiveWafer", 0, 0, 0);
  logicActiveWafer->SetVisAttributes(SiliconVisAtt);

  G4LogicalVolume* logicDeadLayer =
    new G4LogicalVolume(DeadLayer,m_MaterialSilicon,"logicActiveWaferDeadLayer", 0, 0, 0);
  logicDeadLayer->SetVisAttributes(SiliconVisAtt);

  // Set the sensitive volume
  logicActiveWafer->SetSensitiveDetector(m_InnerBarrelScorer);

  // Place the sub volumes in the master volume
  // Last argument is the detector number, used in the scorer to get the
  // revelant information
  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      G4ThreeVector(0,0,0),
      logicPCB,"Tiara_Barrel_PCB",logicBarrelDetector,
      false,0);

  G4ThreeVector WaferPosition(0,
      0.5*(INNERBARREL_PCB_Thickness+INNERBARREL_InertWafer_Thickness)-INNERBARREL_PCB_WaferDepth
      ,0);

  G4ThreeVector DeadLayerPositionF = WaferPosition + G4ThreeVector(0,-INNERBARREL_ActiveWafer_Thickness*0.5-INNERBARREL_ActiveWafer_DeadLayerThickness*0.5,0);

  G4ThreeVector DeadLayerPositionB = WaferPosition + G4ThreeVector(0,INNERBARREL_ActiveWafer_Thickness*0.5+INNERBARREL_ActiveWafer_DeadLayerThickness*0.5,0);



  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      WaferPosition,
      logicActiveWafer,"Barrel_Wafer",
      logicBarrelDetector,false,0);

  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      DeadLayerPositionF,
      logicDeadLayer,"Barrel_WaferDeadLayerFront",
      logicBarrelDetector,false,0);

  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      DeadLayerPositionB,
      logicDeadLayer,"Barrel_WaferDeadLayerBack",
      logicBarrelDetector,false,0);


  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      WaferPosition,
      logicInertWafer,"Barrel_Wafer_GuardRing",
      logicBarrelDetector,false,0);

  // The Distance from target is given by half the lenght of a detector
  // plus the length of a detector inclined by 45 deg.
  G4double DistanceFromTarget = INNERBARREL_PCB_Width*(0.5+sin(45*deg)) ;
  for( unsigned int i = 0; i < 8; i ++){
    // The following build the barrel, with detector one at the top
    // and going clowise looking upstrea

    // Detector are rotate by 45deg with each other
    // Detector 3 [i=2] is perpendicular to positive y-axis, Detector 5 [i=4] is perpendicular to positive x-axis,
    G4RotationMatrix* DetectorRotation =
//by Shuya 180521. This seems to give the right directio of Barrel geometry. Tentative. Wait to be comfirmed. 
      //new G4RotationMatrix(0*deg,0*deg,(360+(2-i)*45)*deg);
      new G4RotationMatrix(0*deg,0*deg,(360+(2-i)*-45)*deg);

    // There center is also rotated by 45deg
    G4ThreeVector DetectorPosition(0,DistanceFromTarget,0);
//by Shuya 180521. This seems to give the right directio of Barrel geometry. Tentative. Wait to be comfirmed. 
    //DetectorPosition.rotate((360+(2-i)*-45)*deg,G4ThreeVector(0,0,1));
    DetectorPosition.rotate((360+(2-i)*45)*deg,G4ThreeVector(0,0,1));

    // Place the Master volume with its two daugther volume at the final place
    new G4PVPlacement(G4Transform3D(*DetectorRotation,DetectorPosition),
        logicBarrelDetector,"Tiara_InnerBarrel_Detector",
        world,false,i+1);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tiara::ConstructOuterBarrel(G4LogicalVolume* world){
  // Tiara Outer Barrel
  // The outer Barrel is identical to the inner barrel but wider in terms of
  // geometry. It feature four non resistive strip on the front face

  // Start by making a full pcb
  // We start by the definition of the point forming a PCB cross section
  vector<G4TwoVector> PCBCrossSection;
  double l1 = OUTERBARREL_PCB_Thickness*0.5/tan(OUTERBARREL_PCB_Bevel1_Theta);
  double l2 = OUTERBARREL_PCB_Thickness*0.5/tan(OUTERBARREL_PCB_Bevel2_Theta);

  PCBCrossSection.push_back(G4TwoVector(OUTERBARREL_PCB_Width/2.-l2,-OUTERBARREL_PCB_Thickness*0.5));
  PCBCrossSection.push_back(G4TwoVector(OUTERBARREL_PCB_Width/2.,0));
  PCBCrossSection.push_back(G4TwoVector(OUTERBARREL_PCB_Width/2.-l1,OUTERBARREL_PCB_Thickness*0.5));

  PCBCrossSection.push_back(G4TwoVector(-OUTERBARREL_PCB_Width/2.+l1,OUTERBARREL_PCB_Thickness*0.5));
  PCBCrossSection.push_back(G4TwoVector(-OUTERBARREL_PCB_Width/2.,0));
  PCBCrossSection.push_back(G4TwoVector(-OUTERBARREL_PCB_Width/2.+l2,-OUTERBARREL_PCB_Thickness*0.5));

  G4ExtrudedSolid* PCBFull =
    new G4ExtrudedSolid("PCBFull",
        PCBCrossSection,
        OUTERBARREL_PCB_Length/2.,
        G4TwoVector(0,0),1,
        G4TwoVector(0,0),1);

  // A box having Wafer dimension but thicker than the PCB
  // Will be used to remove material from the PCB to have space for the wafer

  // Calculate the hole shift within the PCB
  G4ThreeVector HoleShift = G4ThreeVector(
      0,
      0,
      OUTERBARREL_PCB_Offset-(OUTERBARREL_PCB_Length/2-OUTERBARREL_PCB_HoleLength/2));

  G4Box*  HoleShape = new G4Box("HoleShape",
      OUTERBARREL_ActiveWafer_Width/2.,
      OUTERBARREL_PCB_Thickness/2.+0.1*mm,
      OUTERBARREL_PCB_HoleLength/2.);

  G4Box*  WaferShape = new G4Box("WaferShape",
      OUTERBARREL_InertWafer_Width/2.,
      OUTERBARREL_PCB_Thickness/2.,
      OUTERBARREL_InertWafer_Length/2.);

  // The Silicon Wafer itself
  G4Box*  InertWaferFull = new G4Box("InertWaferFull",
      OUTERBARREL_InertWafer_Width/2.,
      OUTERBARREL_ActiveWafer_Thickness/2.,
      OUTERBARREL_InertWafer_Length/2.);

  G4Box*  ActiveWafer = new G4Box("ActiveWafer",
      OUTERBARREL_ActiveWafer_Width/2.,
      OUTERBARREL_ActiveWafer_Thickness/2.,
      OUTERBARREL_ActiveWafer_Length/2.);

  G4Box*  ActiveWaferShape = new G4Box("ActiveWaferShape",
      OUTERBARREL_ActiveWafer_Width/2.,
      OUTERBARREL_PCB_Thickness/2.,
      OUTERBARREL_ActiveWafer_Length/2.);


  // Substracting the hole Shape from the Stock PCB
  G4SubtractionSolid* PCB_1 = new G4SubtractionSolid("PCB_1", PCBFull, HoleShape,
      new G4RotationMatrix,HoleShift);

  // Substracting the wafer space from the Stock PCB
  G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_1, WaferShape,
      new G4RotationMatrix,
      G4ThreeVector(0,OUTERBARREL_PCB_Thickness/2.-OUTERBARREL_PCB_WaferDepth,0));

  // Substract active part from inert part of the Wafer
  G4SubtractionSolid* InertWafer = new G4SubtractionSolid("InertWafer", InertWaferFull, ActiveWaferShape,
      new G4RotationMatrix,
      G4ThreeVector(0,0,0));

  // Master Volume that encompass everything else
  G4LogicalVolume* logicBarrelDetector =
    new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector", 0, 0, 0);
  logicBarrelDetector->SetVisAttributes(G4VisAttributes::Invisible);

  // Sub Volume PCB
  G4LogicalVolume* logicPCB =
    new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
  logicPCB->SetVisAttributes(PCBVisAtt);

  // Sub Volume Wafer
  G4LogicalVolume* logicInertWafer =
    new G4LogicalVolume(InertWafer,m_MaterialSilicon,"logicInertWafer", 0, 0, 0);
  logicInertWafer->SetVisAttributes(GuardRingVisAtt);

  G4LogicalVolume* logicActiveWafer =
    new G4LogicalVolume(ActiveWafer,m_MaterialSilicon,"logicActiveWafer", 0, 0, 0);
  logicActiveWafer->SetVisAttributes(SiliconVisAtt);

  // Set the sensitive detector
  logicActiveWafer->SetSensitiveDetector(m_OuterBarrelScorer);

  // The Distance from target is given by half the lenght of a detector
  // plus the length of a detector inclined by 45 deg.
  G4double DistanceFromTarget = OUTERBARREL_PCB_Width*(0.5+sin(45*deg)) ;

  // Place the sub volumes in the master volume
  // Last argument is the detector number, used in the scorer to get the
  // revelant information
  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      G4ThreeVector(0,0,0),
      logicPCB,"Tiara_OuterBarrel_PCB",logicBarrelDetector,
      false,0);


  G4ThreeVector WaferPosition(0,
      0.5*(OUTERBARREL_PCB_Thickness+OUTERBARREL_ActiveWafer_Thickness)-OUTERBARREL_PCB_WaferDepth
      ,0);

  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      WaferPosition,
      logicActiveWafer,"OuterBarrel_Wafer",
      logicBarrelDetector,false,0);

  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      WaferPosition,
      logicInertWafer,"OuterBarrel_Wafer_GuardRing",
      logicBarrelDetector,false,0);


  // The following build the barrel, with detector one at the top
  // and going clowise looking upstrea
  for( unsigned int i = 0; i < 8; i ++){
    // Detector are rotate by 45deg with each other
    G4RotationMatrix* DetectorRotation =
      new G4RotationMatrix(0*deg,0*deg,i*45*deg);

    // There center is also rotated by 45deg
    G4ThreeVector DetectorPosition(0,DistanceFromTarget,0);
    DetectorPosition.rotate(i*45*deg,G4ThreeVector(0,0,-1));

    // Place the Master volume with its daugthers at the final place
    new G4PVPlacement(G4Transform3D(*DetectorRotation,DetectorPosition),
        logicBarrelDetector,"Tiara_OuterBarrel_Detector",
        world,false,i+1);
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tiara::ConstructHyball(G4LogicalVolume* world){
  vector<G4TwoVector> PCBCrossSection;
  PCBCrossSection.push_back(G4TwoVector(28.108*mm ,-14.551*mm));
  PCBCrossSection.push_back(G4TwoVector(128.808*mm,-66.683*mm));
  PCBCrossSection.push_back(G4TwoVector(163.618*mm,-30.343*mm));
  PCBCrossSection.push_back(G4TwoVector(163.618*mm, 30.941*mm));
  PCBCrossSection.push_back(G4TwoVector(125.718*mm, 73.677*mm));
  PCBCrossSection.push_back(G4TwoVector(28.108*mm , 16.473*mm));

  G4ExtrudedSolid* PCBFull =
    new G4ExtrudedSolid("PCBFull",
        PCBCrossSection,
        0.5*HYBALL_PCB_THICKNESS,
        G4TwoVector(0,0),1,
        G4TwoVector(0,0),1);

  vector<G4TwoVector> WaferCrossSection;
  WaferCrossSection.push_back(G4TwoVector(29.108*mm ,-13.943*mm ));
  WaferCrossSection.push_back(G4TwoVector(123.022*mm,-62.561*mm ));
  WaferCrossSection.push_back(G4TwoVector(137.00*mm ,-24.157*mm ));
  WaferCrossSection.push_back(G4TwoVector(137.00*mm , 24.157*mm ));
  WaferCrossSection.push_back(G4TwoVector(122.677*mm, 63.508*mm ));
  WaferCrossSection.push_back(G4TwoVector(29.108*mm , 15.069*mm));

  G4ExtrudedSolid* WaferFull =
    new G4ExtrudedSolid("WaferFull",
        WaferCrossSection,
        0.5*HYBALL_ActiveWafer_Thickness,
        G4TwoVector(0,0),1,
        G4TwoVector(0,0),1);

  G4ExtrudedSolid* WaferShape =
    new G4ExtrudedSolid("WaferShape",
        WaferCrossSection,
        0.6*HYBALL_PCB_THICKNESS,
        G4TwoVector(0,0),1,
        G4TwoVector(0,0),1);

  // Active Wafer
  G4Tubs* ActiveWafer =
    new G4Tubs("HyballActiveWafer",HYBALL_ActiveWafer_InnerRadius,
        HYBALL_ActiveWafer_OuterRadius,0.5*HYBALL_ActiveWafer_Thickness,
        -0.5*HYBALL_ActiveWafer_Angle,HYBALL_ActiveWafer_Angle);

  G4Tubs* ActiveWaferShape =
    new G4Tubs("HyballActiveWaferShape",HYBALL_ActiveWafer_InnerRadius,
        HYBALL_ActiveWafer_OuterRadius,0.6*HYBALL_ActiveWafer_Thickness,
        -0.5*HYBALL_ActiveWafer_Angle,HYBALL_ActiveWafer_Angle);

 //by Shuya 180219. Hyball DeadLayer
  G4Tubs* Hyball_DeadLayer =
    new G4Tubs("HyballDeadLayer",HYBALL_ActiveWafer_InnerRadius,
        HYBALL_ActiveWafer_OuterRadius,0.5*HYBALL_DeadLayer_Thickness,
        -0.5*HYBALL_ActiveWafer_Angle,HYBALL_ActiveWafer_Angle);


  // Substract Active Wafer from Wafer
  G4SubtractionSolid* InertWafer = new G4SubtractionSolid("Hyball_InertWafer", WaferFull, ActiveWaferShape,
      new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0));


  // Substract Wafer shape from PCB
  G4SubtractionSolid* PCB = new G4SubtractionSolid("Hyball_PCB", PCBFull, WaferShape,
      new G4RotationMatrix,G4ThreeVector(0,0,0));

  // Logic Volume //
  // Logic Mother Volume
  G4LogicalVolume* logicHyball =
    new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicHyball", 0, 0, 0);
  logicHyball->SetVisAttributes(G4VisAttributes::Invisible);

  //  logic PCB
  G4LogicalVolume* logicPCB =
    new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
  logicPCB->SetVisAttributes(PCBVisAtt);

  // logic Inert Wafer
  G4LogicalVolume* logicIW =
    new G4LogicalVolume(InertWafer,m_MaterialSilicon,"logicIW", 0, 0, 0);
  logicIW->SetVisAttributes(GuardRingVisAtt);

  // logic Active Wafer
  G4LogicalVolume* logicAW =
    new G4LogicalVolume(ActiveWafer,m_MaterialSilicon,"logicAW", 0, 0, 0);
  logicAW->SetVisAttributes(SiliconVisAtt);
  logicAW->SetSensitiveDetector(m_HyballScorer);

  //by Shuya 180219
  // logic Hyball DeadLayer
  G4LogicalVolume* logicDL =
    new G4LogicalVolume(Hyball_DeadLayer,m_MaterialSilicon,"logicDL", 0, 0, 0);
  logicDL->SetVisAttributes(SiliconVisAtt);

  // Place all the Piece in the mother volume
  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      G4ThreeVector(0,0,0),
      logicPCB,"Hyball_PCB",
      logicHyball,false,0);

  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      G4ThreeVector(0,0,0),
      logicIW,"Hyball_InertWafer",
      logicHyball,false,0);

  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      G4ThreeVector(0,0,0),
      logicAW,"Hyball_ActiveWafer",
      logicHyball,false,0);

  //by Shuya 180219. Note beacause Hyball is placed in negative direction (such as -147 mm) in world coordinate, deadlayer (front side) positions must be added in positive direction.
  new G4PVPlacement(new G4RotationMatrix(0,0,0),
      G4ThreeVector(0,0,0.5*HYBALL_ActiveWafer_Thickness+0.5*HYBALL_DeadLayer_Thickness),
      logicDL,"Hyball_DeadLayer",
      logicHyball,false,0);

  for(unsigned int i = 0 ; i < m_HyballZ.size() ; i++){
    // Place mother volume
    new G4PVPlacement(new G4RotationMatrix(0,0,m_HyballPhi[i]),
        G4ThreeVector(0,0,m_HyballZ[i]),
        logicHyball,"Hyball",
        world,false,i+1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tiara::ConstructChamber(G4LogicalVolume* world){
  // Vaccum Chamber of Tiara
  // The chamber is made of a central cylinder surrounding the barrel Si
  // Two Cone that expeand out of the central cylinder to let room for Exogam
  // Two outer cylinder surrounding Hyball
  // Hyball is hold on a back plate that close the Diabolo Shaped Chamber

  // Making the Chamber //
  // We make the individual pieces, starting from the inside to the outside
  // Then we merge them together using the a G4AdditionSolid
  // The whole chamber is then placed

  //  Central Tube
  G4Tubs* solidCentralTube =
    new G4Tubs("TiaraChamberCentralTube",CHAMBER_CentralTube_Inner_Radius,
        CHAMBER_CentralTube_Outer_Radius,CHAMBER_CentralTube_Length/2.,
        0*deg,360*deg);

  // Forward-Backward Cones
  G4Cons* solidOuterCone =
    new G4Cons("TiaraChamberOuterCone",CHAMBER_CentralTube_Inner_Radius,
        CHAMBER_CentralTube_Outer_Radius,CHAMBER_OuterCylinder_Inner_Radius,
        CHAMBER_OuterCylinder_Outer_Radius,CHAMBER_OuterCone_Length/2.,
        0*deg,360*deg);

  // Outer Cylinder
  G4Tubs* solidOuterCylinder =
    new G4Tubs("TiaraChamberOuterCylinder",CHAMBER_OuterCylinder_Inner_Radius,
        CHAMBER_OuterCylinder_Outer_Radius,CHAMBER_OuterCylinder_Length/2.,
        0*deg,360*deg);

  // Add the volume together
  G4UnionSolid* solidTiaraChamberStep1 =
    new G4UnionSolid("TiaraChamber", solidCentralTube, solidOuterCone,
        new G4RotationMatrix,
        G4ThreeVector(0,0,CHAMBER_OuterCone_Z_Pos));

  G4UnionSolid* solidTiaraChamberStep2 =
    new G4UnionSolid("TiaraChamber", solidTiaraChamberStep1, solidOuterCone,
        new G4RotationMatrix(0,180*deg,0),
        G4ThreeVector(0,0,-CHAMBER_OuterCone_Z_Pos));

  G4UnionSolid* solidTiaraChamberStep3 =
    new G4UnionSolid("TiaraChamber", solidTiaraChamberStep2, solidOuterCylinder,
        new G4RotationMatrix,
        G4ThreeVector(0,0,CHAMBER_OuterCylinder_Z_Pos));

  G4UnionSolid* solidTiaraChamberStep4 =
    new G4UnionSolid("TiaraChamber", solidTiaraChamberStep3, solidOuterCylinder,
        new G4RotationMatrix,
        G4ThreeVector(0,0,-CHAMBER_OuterCylinder_Z_Pos));

  // Create Logic Volume
  G4LogicalVolume* logicTiaraChamber =
    new G4LogicalVolume(solidTiaraChamberStep4,m_MaterialAl,"logicTiaraChamber", 0, 0, 0);

  // Visual Attribute
  G4VisAttributes* ChamberVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.4,0.5,0.2));

  logicTiaraChamber->SetVisAttributes(ChamberVisAtt);

  // Place the whole chamber
  new G4PVPlacement(new G4RotationMatrix(0,0,0), G4ThreeVector(0,0,0),
      logicTiaraChamber,"TiaraChamber",world,false,0);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tiara::InitializeMaterial(){
  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialAl = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Tiara::Construct(){
  return (NPS::VDetector*) new Tiara();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_nps_tiara{
  public:
    proxy_nps_tiara(){
      NPS::DetectorFactory::getInstance()->AddToken("Tiara","Tiara");
      NPS::DetectorFactory::getInstance()->AddDetector("Tiara",Tiara::Construct);
    }
};

proxy_nps_tiara p_nps_tiara;
}
