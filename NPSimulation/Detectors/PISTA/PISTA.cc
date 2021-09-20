/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  PISTA simulation                             *
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
#include "G4Trap.hh"
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

// NPTool header
#include "PISTA.hh"
#include "DSSDScorers.hh"
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
namespace PISTA_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 0.2*ns ;
  const double ResoEnergy = 0.015*MeV ;
  const double DE_ResoEnergy = 0.015*MeV ;

  // Trapezoid dimension
  const double TrapezoidBaseLarge = 74.1*mm;
  //const double TrapezoidBaseLarge = 78.1*mm;
  const double TrapezoidBaseSmall = 39.3*mm;
  //const double TrapezoidBaseSmall = 43.3*mm;
  const double TrapezoidHeight = 57.8*mm;
  //const double TrapezoidHeight = 61.8*mm;
  const double TrapezoidLength = 1*cm;
  const double FirstStageThickness = 100*um;
  const double SecondStageThickness = 1*mm;
  const double DistanceBetweenSi = 5*mm;
  //const double FirstStageNbrOfStrips = 97;
  //const double SecondStageNbrOfStrips = 122;
}
using namespace PISTA_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// PISTA Specific Method
PISTA::PISTA(){
  m_Event = new TPISTAData() ;
  m_FirstStageScorer = 0;
  m_SecondStageScorer = 0;
  m_TrapezoidDetector = 0;
}

PISTA::~PISTA(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PISTA::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PISTA::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* PISTA::BuildTrapezoidDetector(){
  if(!m_TrapezoidDetector){
    // Definittion of the volume containing the sensitive detectors
    G4Trap* solidTrapezoid = new G4Trap("PISTA",
        TrapezoidLength*0.5, 0*deg, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg, 
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg);
    G4LogicalVolume* logicTrapezoid = new G4LogicalVolume(solidTrapezoid,
        MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"),
        "PISTA",
        0,0,0);

    G4VisAttributes* TrapezoidVisAtt = new G4VisAttributes(G4Colour(0,0,0,0.5));
    TrapezoidVisAtt->SetForceWireframe(true);
    logicTrapezoid->SetVisAttributes(TrapezoidVisAtt);

    // First stage silicon detector
    G4ThreeVector positionFirstStage = G4ThreeVector(0,0,-4*mm);

    G4Trap* solidFirstStage = new G4Trap("solidFirstSatge",
        FirstStageThickness*0.5, 0*deg, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg);
    G4LogicalVolume* logicFirstStage = new G4LogicalVolume(solidFirstStage,
        MaterialManager::getInstance()->GetMaterialFromLibrary("Si"),
        "logicFirstStage",
        0,0,0);
    new G4PVPlacement(0,
        positionFirstStage,
        logicFirstStage,
        "PISTA_FirstStage",
        logicTrapezoid,
        false,
        0);
    // Set First Stage sensitive
    logicFirstStage->SetSensitiveDetector(m_FirstStageScorer);

    // Visualisation of First Stage strips
    //G4VisAttributes* FirstStageVisAtt = new G4VisAttributes(G4Colour(0.2,0.8,0.5));
    G4VisAttributes* FirstStageVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.55,0.8));
    logicFirstStage->SetVisAttributes(FirstStageVisAtt);

    //////
    // Second stage silicon detector
    G4ThreeVector positionSecondStage = G4ThreeVector(0,0,-0.5*TrapezoidLength+DistanceBetweenSi);

    G4Trap* solidSecondStage = new G4Trap("solidSecondSatge",
        SecondStageThickness*0.5, 0*deg, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg,
        TrapezoidHeight*0.5, TrapezoidBaseLarge*0.5, TrapezoidBaseSmall*0.5, 0*deg);
    G4LogicalVolume* logicSecondStage = new G4LogicalVolume(solidSecondStage,
        MaterialManager::getInstance()->GetMaterialFromLibrary("Si"),
        "logicSecondStage",
        0,0,0);
    new G4PVPlacement(0,
        positionSecondStage,
        logicSecondStage,
        "PISTA_SecondStage",
        logicTrapezoid,
        false,
        0);
    // Set Second Stage sensitive
    logicSecondStage->SetSensitiveDetector(m_SecondStageScorer);

    // Visualisation of Second Stage strips
    G4VisAttributes* SecondStageVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    logicSecondStage->SetVisAttributes(SecondStageVisAtt);


    m_TrapezoidDetector = logicTrapezoid;
  }
  return m_TrapezoidDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void PISTA::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("PISTA");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetector(R,Theta,Phi);
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
void PISTA::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = TrapezoidHeight*0.5 + m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*PISTA_NS::TrapezoidLength*0.5;
    // Building Detector reference frame
    G4double ii = cos(m_Theta[i]) * cos(m_Phi[i]);
    G4double jj = cos(m_Theta[i]) * sin(m_Phi[i]);
    G4double kk = -sin(m_Theta[i]);
    G4ThreeVector Y(ii,jj,kk);
    G4ThreeVector w = Det_pos.unit();
    G4ThreeVector u = w.cross(Y);
    G4ThreeVector v = w.cross(u);
    v = v.unit();
    u = u.unit();

    G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
        BuildTrapezoidDetector(),
        "PISTA",world,false,i+1);

  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void PISTA::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("PISTA")){
    pTree->Branch("PISTA", "TPISTAData", &m_Event) ;
  }
  pTree->SetBranchAddress("PISTA", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void PISTA::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // First Stage scorer
  DSSDScorers::PS_Rectangle* FirstStageScorer= (DSSDScorers::PS_Rectangle*) m_FirstStageScorer->GetPrimitive(0);

  unsigned int sizeFront = FirstStageScorer->GetLengthMult(); 
  for(unsigned int i = 0 ; i < sizeFront ; i++){
    double Energy = RandGauss::shoot(FirstStageScorer->GetEnergyLength(i), DE_ResoEnergy);   
    if(Energy>EnergyThreshold){
      double Time = RandGauss::shoot(FirstStageScorer->GetTimeLength(i), ResoTime);
      int DetNbr  = FirstStageScorer->GetDetectorLength(i);
      int StripFront = FirstStageScorer->GetStripLength(i);
      m_Event->SetFirstStageXE(DetNbr, StripFront, Energy);
      m_Event->SetFirstStageXT(DetNbr, StripFront, Time);
    }
  }
  unsigned int sizeBack = FirstStageScorer->GetWidthMult(); 
  for(unsigned int i = 0 ; i < sizeBack ; i++){
    double Energy = RandGauss::shoot(FirstStageScorer->GetEnergyWidth(i), DE_ResoEnergy);   
    if(Energy>EnergyThreshold){
      double Time = RandGauss::shoot(FirstStageScorer->GetTimeWidth(i), ResoTime);
      int DetNbr  = FirstStageScorer->GetDetectorWidth(i);
      int StripBack = FirstStageScorer->GetStripWidth(i);
      m_Event->SetFirstStageYE(DetNbr, StripBack, Energy);
      m_Event->SetFirstStageYT(DetNbr, StripBack, Time);
    }
  }
  FirstStageScorer->clear();

  ///////////
  // Second Stage scorer
  DSSDScorers::PS_Rectangle* SecondStageScorer= (DSSDScorers::PS_Rectangle*) m_SecondStageScorer->GetPrimitive(0);

  unsigned int sizeFrontSecondStage = SecondStageScorer->GetLengthMult(); 
  for(unsigned int i = 0 ; i < sizeFrontSecondStage ; i++){
    double Energy = RandGauss::shoot(SecondStageScorer->GetEnergyLength(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      double Time = RandGauss::shoot(SecondStageScorer->GetTimeLength(i), ResoTime);
      int DetNbr  = SecondStageScorer->GetDetectorLength(i);
      int StripFront = SecondStageScorer->GetStripLength(i);
      m_Event->SetSecondStageXE(DetNbr, StripFront, Energy);
      m_Event->SetSecondStageXT(DetNbr, StripFront, Time);
    }
  }
  unsigned int sizeBackSecondStage = SecondStageScorer->GetWidthMult(); 
  for(unsigned int i = 0 ; i < sizeBackSecondStage ; i++){
    double Energy = RandGauss::shoot(SecondStageScorer->GetEnergyWidth(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      double Time = RandGauss::shoot(SecondStageScorer->GetTimeWidth(i), ResoTime);
      int DetNbr  = SecondStageScorer->GetDetectorWidth(i);
      int StripBack = SecondStageScorer->GetStripWidth(i);
      m_Event->SetSecondStageYE(DetNbr, StripBack, Energy);
      m_Event->SetSecondStageYT(DetNbr, StripBack, Time);
    }
  }
  SecondStageScorer->clear();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void PISTA::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_FirstStageScorer = CheckScorer("FirstStageScorer",already_exist) ;
  m_SecondStageScorer = CheckScorer("SecondStageScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  G4VPrimitiveScorer* FirstStageScorer = new DSSDScorers::PS_Rectangle("FirstStageScorer",1,
      TrapezoidBaseLarge,
      TrapezoidHeight,
      1,97);
  G4VPrimitiveScorer* SecondStageScorer = new DSSDScorers::PS_Rectangle("SecondStageScorer",1,
      TrapezoidBaseLarge,
      TrapezoidHeight,
      62,1);

  G4VPrimitiveScorer* InteractionFirstStage = new InteractionScorers::PS_Interactions("InteractionFirstStage",ms_InterCoord,0);
  G4VPrimitiveScorer* InteractionSecondStage = new InteractionScorers::PS_Interactions("InteractionSecondStage",ms_InterCoord,0);

  // Register it to the multifunctionnal detector
  m_FirstStageScorer->RegisterPrimitive(FirstStageScorer);
  m_FirstStageScorer->RegisterPrimitive(InteractionFirstStage);
  m_SecondStageScorer->RegisterPrimitive(SecondStageScorer);
  m_SecondStageScorer->RegisterPrimitive(InteractionSecondStage);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_FirstStageScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SecondStageScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* PISTA::Construct(){
  return  (NPS::VDetector*) new PISTA();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_PISTA{
    public:
      proxy_nps_PISTA(){
        NPS::DetectorFactory::getInstance()->AddToken("PISTA","PISTA");
        NPS::DetectorFactory::getInstance()->AddDetector("PISTA",PISTA::Construct);
      }
  };

  proxy_nps_PISTA p_nps_PISTA;
}
