/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Nebula simulation                                   *
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
#include "Nebula.hh"
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
namespace Nebula_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime     = 4.5*ns ;
  const double ResoEnergy   = 1.0*MeV ;
  const double ModuleWidth  = 120*mm ;
  const double ModuleLength = 120*mm ;
  const double ModuleHeight = 1800*mm ;
  const double InterModule  = 1*mm ;
  const double VetoWidth    = 320*mm ;
  const double VetoLength   = 10*mm ;
  const double VetoHeight   = 1900*mm ;
  const double InterVeto    = 1*mm ;
  const int    VetoPerWall  = 12;
  const double WallToVeto   = 10*cm;

  const string Material = "BC400";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Nebula Specific Method
Nebula::Nebula(){
  m_Event = new TNebulaData() ;
  m_ModuleScorer = 0;
  m_VetoScorer = 0;
  m_Module = 0;
  m_Veto = 0;


  // RGB Color + Transparency
  m_VisModule = new G4VisAttributes(G4Colour(0.3, 0.3, 1, 1));   
  m_VisVeto   = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.2));   
  m_VisPMT    = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1, 1));   
  m_VisFrame  = new G4VisAttributes(G4Colour(0, 0.3, 1, 0.5));   

}

Nebula::~Nebula(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Nebula::AddWall(G4ThreeVector Pos, int NbrModule, bool Veto, bool Frame){
  // Convert the Pos value to R theta Phi as Spherical coordinate is easier in G4 
  m_Pos.push_back(Pos);
  m_NbrModule.push_back(NbrModule);
  m_HasVeto.push_back(Veto);
  m_HasFrame.push_back(Frame);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Nebula::BuildModule(){
  if(!m_Module){
    G4Box* box = new G4Box("Nebula_Module",Nebula_NS::ModuleWidth*0.5,
        Nebula_NS::ModuleHeight*0.5,Nebula_NS::ModuleLength*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Nebula_NS::Material);
    m_Module = new G4LogicalVolume(box,DetectorMaterial,"logic_Nebula_Module",0,0,0);
    m_Module->SetVisAttributes(m_VisModule);
    m_Module->SetSensitiveDetector(m_ModuleScorer);
  }
  return m_Module;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Nebula::BuildVeto(){
  if(!m_Veto){
    G4Box* box = new G4Box("Nebula_Veto",Nebula_NS::VetoWidth*0.5,
        Nebula_NS::VetoHeight*0.5,Nebula_NS::VetoLength*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Nebula_NS::Material);
    m_Veto = new G4LogicalVolume(box,DetectorMaterial,"logic_Nebula_Veto",0,0,0);
    m_Veto->SetVisAttributes(m_VisVeto);
    m_Veto->SetSensitiveDetector(m_VetoScorer);
  }
  return m_Veto;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Nebula::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Nebula");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"Pos","NumberOfModule","Veto","Frame"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Nebula " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("Pos","mm"));
      int NbrModule = blocks[i]->GetInt("NumberOfModule");
      bool Veto = blocks[i]->GetInt("Veto");
      bool Frame= blocks[i]->GetInt("Frame");
      AddWall(Pos,NbrModule,Veto,Frame);
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
void Nebula::ConstructDetector(G4LogicalVolume* world){
  unsigned int nbrM = 1 ;
  unsigned int nbrV = 1 ;
  for (unsigned short i = 0 ; i < m_Pos.size() ; i++) {
    for (unsigned short m = 0 ; m < m_NbrModule[i] ; m++) {
      G4RotationMatrix* Rot = new G4RotationMatrix();
      double offset = (Nebula_NS::ModuleWidth+Nebula_NS::InterModule)*(-m_NbrModule[i]*0.5+m)+Nebula_NS::ModuleWidth*0.5;
      G4ThreeVector Offset(offset,0,0);
      new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
          BuildModule(),
          "NebulaModule",world,false,nbrM++);
     }

    if(m_HasVeto[i]){
      for (unsigned short m = 0 ; m < Nebula_NS::VetoPerWall ; m++) {
        G4RotationMatrix* Rot = new G4RotationMatrix();
        double offset = (Nebula_NS::VetoWidth+Nebula_NS::InterVeto)*(-Nebula_NS::VetoPerWall*0.5+m)+Nebula_NS::VetoWidth*0.5;
        G4ThreeVector Offset(offset,0,-Nebula_NS::WallToVeto);
        new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
          BuildVeto(),
          "NebulaVeto",world,false,nbrV++);
     }


      }
     }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Nebula::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Nebula")){
    pTree->Branch("Nebula", "TNebulaData", &m_Event) ;
  }
  pTree->SetBranchAddress("Nebula", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Nebula::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Module scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_ModuleScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Nebula_NS::ResoEnergy);
    if(Energy>Nebula_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Nebula_NS::ResoTime);
      int DetectorNbr = level[0];
      //m_Event->SetEnergy(DetectorNbr,Energy);
      //m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Nebula::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_ModuleScorer = CheckScorer("NebulaModuleScorer",already_exist) ;
  m_VetoScorer = CheckScorer("NebulaVetoScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialise
  // Module 
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* ModuleCalorimeter= new CalorimeterScorers::PS_Calorimeter("ModuleCalorimeter",level, 0) ;
  G4VPrimitiveScorer* ModuleInteraction= new InteractionScorers::PS_Interactions("ModuleInteraction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_ModuleScorer->RegisterPrimitive(ModuleCalorimeter);
  m_ModuleScorer->RegisterPrimitive(ModuleInteraction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ModuleScorer) ;

  // Veto 
  G4VPrimitiveScorer* VetoCalorimeter= new CalorimeterScorers::PS_Calorimeter("VetoCalorimeter",level, 0) ;
  G4VPrimitiveScorer* VetoInteraction= new InteractionScorers::PS_Interactions("VetoInteraction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_VetoScorer->RegisterPrimitive(VetoCalorimeter);
  m_VetoScorer->RegisterPrimitive(VetoInteraction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_VetoScorer) ;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Nebula::Construct(){
  return  (NPS::VDetector*) new Nebula();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Nebula{
    public:
      proxy_nps_Nebula(){
        NPS::DetectorFactory::getInstance()->AddToken("Nebula","Nebula");
        NPS::DetectorFactory::getInstance()->AddDetector("Nebula",Nebula::Construct);
      }
  };

  proxy_nps_Nebula p_nps_Nebula;
}
