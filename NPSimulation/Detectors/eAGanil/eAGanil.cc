/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  eAGanil simulation                                  *
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
#include "eAGanil.hh"
#include "CalorimeterScorers.hh"
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
namespace eAGanil_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 4.5*ns ;
  const double ResoEnergy = 1.0*MeV ;
  const double Thickness = 500*mm ;
  const string Material = "Pb";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// eAGanil Specific Method
eAGanil::eAGanil(){
  m_Event = new TeAGanilData() ;
  m_eAGanilScorer = 0;

  m_Length=0;

  // RGB Color + Transparency
  m_VisDetector = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));   
  m_VisTrap = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3, 0.5));   
}

eAGanil::~eAGanil(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void eAGanil::AddDetector(double  R, double  Theta, double  Phi, double EntranceWidth,double EntranceHeigh,double MR){
  m_SpecR.push_back(R);
  m_SpecTheta.push_back(Theta);
  m_SpecPhi.push_back(Phi);
  m_SpecEntranceWidth.push_back(EntranceWidth);
  m_SpecEntranceHeigh.push_back(EntranceHeigh);
  m_SpecMomentumResolution.push_back(MR);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void eAGanil::SetTrap(double Length, double  InnerRadius, double  OuterRadius, double BladesThickness, int NumberOfBlades,double  Phi, double WindowsThickness){
  m_Length=Length;
  m_InnerRadius=InnerRadius; 
  m_OuterRadius=OuterRadius;
  m_BladesThickness=BladesThickness;
  m_NumberOfBlades=NumberOfBlades;
  m_Phi=Phi;       
  m_WindowsThickness=WindowsThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* eAGanil::BuildDetector(unsigned int i){
  G4Box* box = new G4Box("eAGanil_Box",m_SpecEntranceHeigh[i]*0.5,
      m_SpecEntranceWidth[i]*0.5,eAGanil_NS::Thickness*0.5);

  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(eAGanil_NS::Material);
  G4LogicalVolume* Detector = new G4LogicalVolume(box,DetectorMaterial,"logic_eAGanil_spec",0,0,0);
  Detector->SetVisAttributes(m_VisDetector);
  Detector->SetSensitiveDetector(m_eAGanilScorer);
  return Detector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* eAGanil::BuildTrap(){
  G4Tubs* tubs = new G4Tubs("eAGanil_Trap",m_InnerRadius*0.9,m_OuterRadius*1.1,m_Length*0.5,0,360*deg );

  G4Material* VacuumMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4LogicalVolume* Trap = new G4LogicalVolume(tubs,VacuumMaterial,"logic_eAGanil_trap",0,0,0);
  Trap->SetVisAttributes(G4VisAttributes::Invisible);


  G4Box* box = new G4Box("eAGanil_Blades",
      m_BladesThickness*0.5,
      (m_OuterRadius-m_InnerRadius)*0.5,
      m_Length*0.5
      );

  G4Material* TrapMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4LogicalVolume* Blades = new G4LogicalVolume(box,TrapMaterial,"logic_eAGanil_trap",0,0,0);
  Blades->SetVisAttributes(m_VisTrap);
  G4RotationMatrix* Rot = new G4RotationMatrix();
  G4ThreeVector Pos(0,(m_OuterRadius-m_InnerRadius)*0.5+m_InnerRadius,0);
  for(unsigned int i = 0 ; i < m_NumberOfBlades ; i++){
    Rot->rotateZ(360.*deg/m_NumberOfBlades);
    Pos.rotateZ(360.*deg/m_NumberOfBlades);
   new G4PVPlacement(G4Transform3D(*Rot,Pos),
        Blades,
        "eAGanil",Trap,false,1);
    }
 
  return Trap;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void eAGanil::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("eAGanil","Spectrometer");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> sphe = {"R","Theta","Phi","EntranceWidth","EntranceHeight","MomentumResolution"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  eAGanil " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      double EW = blocks[i]->GetDouble("EntranceWidth","cm");
      double EH = blocks[i]->GetDouble("EntranceHeight","cm");
      double MR = blocks[i]->GetDouble("MomentumResolution","void"); 
      AddDetector(R,Theta,Phi,EW,EH,MR);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  blocks = parser.GetAllBlocksWithTokenAndValue("eAGanil","Trap");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> trap= {"Length","InnerRadius", "OuterRadius","BladesThickness","NumberOfBlades","Phi","WindowsThickness"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(trap)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  eAGanil " << i+1 <<  endl;

      double L = blocks[i]->GetDouble("Length","mm");
      double iR = blocks[i]->GetDouble("InnerRadius","mm");
      double oR = blocks[i]->GetDouble("OuterRadius","mm");
      double fT = blocks[i]->GetDouble("BladesThickness","mm");
      int    nF = blocks[i]->GetInt("NumberOfBlades");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      double wT = blocks[i]->GetDouble("WindowsThickness","mm");
      SetTrap(L,iR,oR,fT,nF,Phi,wT);
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
void eAGanil::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_SpecR.size() ; i++) {
    G4double wX = m_SpecR[i] * sin(m_SpecTheta[i] ) * cos(m_SpecPhi[i] ) ;
    G4double wY = m_SpecR[i] * sin(m_SpecTheta[i] ) * sin(m_SpecPhi[i] ) ;
    G4double wZ = m_SpecR[i] * cos(m_SpecTheta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*eAGanil_NS::Thickness*0.5;
    // Building Detector reference frame
    G4double ii = cos(m_SpecTheta[i]) * cos(m_SpecPhi[i]);
    G4double jj = cos(m_SpecTheta[i]) * sin(m_SpecPhi[i]);
    G4double kk = -sin(m_SpecTheta[i]);
    G4ThreeVector Y(ii,jj,kk);
    G4ThreeVector w = Det_pos.unit();
    G4ThreeVector u = w.cross(Y);
    G4ThreeVector v = w.cross(u);
    v = v.unit();
    u = u.unit();

    G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
        BuildDetector(i),
        "eAGanil",world,false,i+1);
  }
  if(m_Length){
    G4RotationMatrix* Rot = new G4RotationMatrix();
    Rot->rotateZ(m_Phi);
    
    new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
        BuildTrap(),
        "eAGanil",world,false,1);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void eAGanil::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("eAGanil")){
    pTree->Branch("eAGanil", "TeAGanilData", &m_Event) ;
  }
  pTree->SetBranchAddress("eAGanil", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void eAGanil::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_eAGanilScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),eAGanil_NS::ResoEnergy);
    if(Energy>eAGanil_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),eAGanil_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void eAGanil::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_eAGanilScorer = CheckScorer("eAGanilScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_eAGanilScorer->RegisterPrimitive(Calorimeter);
  m_eAGanilScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_eAGanilScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* eAGanil::Construct(){
  return  (NPS::VDetector*) new eAGanil();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_eAGanil{
    public:
      proxy_nps_eAGanil(){
        NPS::DetectorFactory::getInstance()->AddToken("eAGanil","eAGanil");
        NPS::DetectorFactory::getInstance()->AddDetector("eAGanil",eAGanil::Construct);
      }
  };

  proxy_nps_eAGanil p_nps_eAGanil;
}
