/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Catana simulation                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * Geometry of crystal based on official catana simulation from Samurai      *
 * Collaboration package 5.2                                                 *
 * http://be.nucl.ap.titech.ac.jp/~nebula/simulator.php                      *
 * Thanks to Togano-san                                                      *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"

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
#include "Catana.hh"
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
namespace Catana_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 4.5*ns ;
  const double ResoEnergy = 0.08*MeV ;
  double Length = 600*mm ;
  string Material = "CsI";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Catana Specific Method
Catana::Catana(){
  m_Event = new TCatanaData() ;
  m_CatanaScorer = 0;
  m_DetectorType1 = 0;
  m_DetectorType2 = 0;
  m_DetectorType3 = 0;
  m_DetectorType4 = 0;
  m_DetectorType5 = 0;

  // RGB Color + Transparency
  m_VisCrystal1 = new G4VisAttributes(G4Colour(0.8, 0.3, 0.3, 1));   
  m_VisCrystal2 = new G4VisAttributes(G4Colour(0.3, 0.8, 0.3, 1));   
  m_VisCrystal3 = new G4VisAttributes(G4Colour(0.3, 0.3, 0.8, 1));   
  m_VisCrystal4 = new G4VisAttributes(G4Colour(0.3, 0.8, 0.8, 1));   
  m_VisCrystal5 = new G4VisAttributes(G4Colour(0.8, 0.3, 0.8, 1));   
  m_VisHousing = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3, 0.2));   
  m_VisReflector = new G4VisAttributes(G4Colour(1, 1, 1, 0.1));   

}

Catana::~Catana(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Catana::AddDetector(double X,double Y, double Z, double Theta, double Phi, int ID, int Type,double Rshift){
  m_X.push_back(X);
  m_Y.push_back(Y);
  m_Z.push_back(Z);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_ID.push_back(ID);
  m_Type.push_back(Type);
  m_Rshift.push_back(Rshift);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Catana::ReadCSV(string path,double Rshift){
  std::ifstream csv(path); 
  if(!csv.is_open()){
    std::ostringstream message;
    message << "Catana csv file " << path << " not found " << std::endl;
  }


  int ID, type,layer;
  double X,Y,Z,Theta,Phi;
  string buffer;
  // ignore first line
  getline(csv,buffer);
  while(csv >> ID >> buffer >> type >> buffer >> layer >> buffer >> X >> buffer >> Y >> buffer >> Z >> buffer >> Theta >> buffer >> Phi){
      if(type<6)
      AddDetector(X,Y,Z,Theta*deg,Phi*deg,ID,type,Rshift);
      else{
        // ignore other type for which I don't have the geometry
        }
  }

  return;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Catana::BuildDetector(int Type){
  // Check if the requested detector already exist
  if(Type==1 && m_DetectorType1)
      return m_DetectorType1;
  
  else if(Type==2 && m_DetectorType2)
      return m_DetectorType2;
  
  else if(Type==3 && m_DetectorType3)
      return m_DetectorType3;

  else if(Type==4 && m_DetectorType4)
      return m_DetectorType4;

  else if(Type==5 && m_DetectorType5)
      return m_DetectorType5;


  // Define the requested geometry
    double x1,x2,x3,x4,y1,y2,z,crystalz;
    //int pmttype; // 1: H7195, 2:H11934
    double seal_dt = 12*mm, housing_dt = 0.5*mm, reflector_dt = 0.5*mm;

    if(Type == 1){ // crystal type 1
      x1 = 62.3*mm; x2 = 62.3*mm; x3 = 95.7*mm; y1 = 36.6*mm; y2 = 56.3*mm;
      z = 107*mm; crystalz = 93.*mm; //pmttype = 1;
    }
    if(Type == 2){ // crystal type 2
      x1 = 57.1*mm; x2 = 63.6*mm; x3 = 84.5*mm; y1 = 34.9*mm; y2 = 55.4*mm;
      z = 117*mm; crystalz = 103.*mm; //pmttype = 1;
    }
    if(Type == 3){ // crystal type 3
      x1 = 49.7*mm; x2 = 58.5*mm; x3 = 74.9*mm; y1 = 38.3*mm; y2 = 64.7*mm;
      z = 137*mm; crystalz = 123.*mm; //pmttype = 1;
    }
    if(Type == 4){ // crystal type 4
      x1 = 40.0*mm; x2 = 50.2*mm; x3 = 60.3*mm; y1 = 38.3*mm; y2 = 66.4*mm;
      z = 152*mm; crystalz = 138.5*mm; //pmttype = 1;
    }
    if(Type == 5){ // crystal type 5
      x1 = 28.4*mm; x2 = 39.7*mm; x3 = 41.5*mm; y1 = 38.3*mm; y2 = 69.9*mm;
      z = 162*mm; crystalz = 148.5*mm; //pmttype = 2;
    }
    x4 = x3 + (y2/y1)*(x2-x1); // planarity condition
    Double_t beta1 = 90*deg + std::atan((x2-x1)/(y1*2));// should be
    Double_t beta2 = 90*deg + std::atan((x4-x3)/(y2*2));// beta1 = beta2
    if(std::abs(beta1-beta2)>1e-4){
      std::cout << "Housing type " << Type << " is not planar" << std::endl;
    }

    // Define Material
    G4Material* CsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
    G4Material* Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* Vacuum= MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    G4Material* Teflon= MaterialManager::getInstance()->GetMaterialFromLibrary("G4_TEFLON");
    G4RotationMatrix* Rot= new G4RotationMatrix();


    // Parameters for G4Trap
    double pDz, pDy1, pDy2, pDx1, pDx2, pDx3, pDx4;
    double pTheta=0.*deg, pPhi=0.*deg, pAlp1=0.*deg, pAlp2=0*deg;

    // housing outside
    pDz  = z*0.5; 
    pDy1 = y1*0.5;
    pDy2 = y2*0.5;
    pDx1 = x1*0.5;
    pDx2 = x2*0.5;
    pDx3 = x3*0.5;
    pDx4 = x4*0.5;
    
    G4Trap *solidHousingOUT = new G4Trap("solidHousingOut", pDz, pTheta, pPhi,
					 pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					 pDx4, pAlp2);

    G4LogicalVolume* result = 0;
    if(Type==1){
      m_DetectorType1  = new G4LogicalVolume(solidHousingOUT,
        Vacuum,
        "logicDetectorType1",
        0,0,0);
      result = m_DetectorType1; 
      }
    else if(Type==2){
      m_DetectorType2  = new G4LogicalVolume(solidHousingOUT,
        Vacuum,
        "logicDetectorType2",
        0,0,0);
      result = m_DetectorType2; 
      }


    else if(Type==3){
      m_DetectorType3  = new G4LogicalVolume(solidHousingOUT,
        Vacuum,
        "logicDetectorType3",
        0,0,0);
      result = m_DetectorType3; 
      }


    else if(Type==4){
      m_DetectorType4  = new G4LogicalVolume(solidHousingOUT,
        Vacuum,
        "logicDetectorType4",
        0,0,0);
      result = m_DetectorType4; 
      }
   
   else if(Type==5){
      m_DetectorType5  = new G4LogicalVolume(solidHousingOUT,
        Vacuum,
        "logicDetectorType5",
        0,0,0);
      result = m_DetectorType5; 
      }


    result->SetVisAttributes(G4VisAttributes::Invisible);


    // -- Al housing inside
    double length = z;
    double ax1 = pDx1;
    double bx1 = pDx3;
    double ax2 = pDx2;
    double bx2 = pDx4;
    double ay1 = pDy1;
    double by1 = pDy2;

    pDz = (length - seal_dt - housing_dt)/2.;
    pDy1 = (by1-ay1)/length * housing_dt + ay1 - housing_dt;
    pDx1 = (bx1-ax1)/length * housing_dt + ax1 - housing_dt;
    pDx2 = (bx2-ax2)/length * housing_dt + ax2 - housing_dt;
    pDy2 = (by1-ay1)/length * (length - seal_dt) + ay1 - housing_dt;
    pDx3 = (bx1-ax1)/length * (length - seal_dt) + ax1 - housing_dt;
    //pDx4 = (bx2-ax2)/length * (length - seal_dt) + ax2 - housing_dt;
    pDx4 = pDx3 + (pDy2 / pDy1)*(pDx2 - pDx1); // planarity condition
    
    G4Trap* solidHousingIN = new G4Trap("solidHousingIN", pDz, pTheta, pPhi,
					pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					pDx4, pAlp2);
 

    double offset = -(length*0.5 - pDz - housing_dt);
    G4SubtractionSolid* solidHousing = 
      new G4SubtractionSolid("solidHousing", solidHousingOUT,// mother
			     solidHousingIN, Rot,
			     G4ThreeVector(0.,0.,offset));

    G4LogicalVolume* LogicHousing = new G4LogicalVolume(solidHousing, Al,"logicHousing",0,0,0);
    LogicHousing->SetVisAttributes(m_VisHousing);

    // -- Crystal --
    double space = 2.*mm; // space btw. crystal and housing
    length = pDz * 2.; // housing inner z length
    ax1 = pDx1;
    bx1 = pDx3;
    ax2 = pDx2;
    bx2 = pDx4;
    ay1 = pDy1;
    by1 = pDy2;

    pDz = crystalz*0.5;    
    pDy1 = (by1-ay1)/length * reflector_dt + ay1 - space;
    pDx1 = (bx1-ax1)/length * reflector_dt + ax1 - space;
    pDx2 = (bx2-ax2)/length * reflector_dt + ax2 - space;
    pDy2 = (by1-ay1)/length * (reflector_dt + crystalz) + ay1 - space;
    pDx3 = (bx1-ax1)/length * (reflector_dt + crystalz) + ax1 - space;
    //pDx4 = (bx2-ax2)/length * (reflector_dt + crystalz) + ax2 - space;
    pDx4 = pDx3 + (pDy2 / pDy1)*(pDx2 - pDx1); // planarity condition
    
    G4Trap* solidCrystal = new G4Trap("solidCrystal", pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);

    G4LogicalVolume* LogicCrystal = new G4LogicalVolume(solidCrystal,// solid
					       CsI,
					       "SolidCrystal",
					       0,0,0);
    if(Type==1)
      LogicCrystal->SetVisAttributes(m_VisCrystal1);
    
    else if(Type==2)
      LogicCrystal->SetVisAttributes(m_VisCrystal2);
    
    else if(Type==3)
      LogicCrystal->SetVisAttributes(m_VisCrystal3);

    else if(Type==4)
      LogicCrystal->SetVisAttributes(m_VisCrystal4);

    else if(Type==5)
      LogicCrystal->SetVisAttributes(m_VisCrystal5);


    LogicCrystal->SetSensitiveDetector(m_CatanaScorer);

    // -- Teflon reflector
    length = crystalz;
    ax1 = pDx1;
    bx1 = pDx3;
    ax2 = pDx2;
    bx2 = pDx4;
    ay1 = pDy1;
    by1 = pDy2;

    pDz = crystalz*0.5 + reflector_dt;
    pDy1 = (by1-ay1)/length * -reflector_dt + ay1 + reflector_dt;
    pDx1 = (bx1-ax1)/length * -reflector_dt + ax1 + reflector_dt;
    pDx2 = (bx2-ax2)/length * -reflector_dt + ax2 + reflector_dt;
    pDy2 = (by1-ay1)/length * (reflector_dt + crystalz) + ay1 + reflector_dt;
    pDx3 = (bx1-ax1)/length * (reflector_dt + crystalz) + ax1 + reflector_dt;
    //pDx4 = (bx2-ax2)/length * (reflector_dt + crystalz) + ax2 + reflector_dt;
    pDx4 = pDx3 + (pDy2 / pDy1)*(pDx2 - pDx1); // planarity condition

    G4Trap* solidReflector = new G4Trap("solidReflector",pDz, pTheta,
					   pPhi, pDy1, pDx1, pDx2, pAlp1,
					   pDy2, pDx3, pDx4, pAlp2);
 
    G4LogicalVolume* LogicReflector = new G4LogicalVolume(solidReflector, Teflon,
						 "logicReflector",
						 0,0,0);

    LogicReflector->SetVisAttributes(m_VisReflector);

    // Place volume in each other:
     new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
        LogicHousing,
        "CatanaHousing",result,false,0);
     

     m_Zoffset[Type] = (z - crystalz)*0.5 - housing_dt - reflector_dt;

     new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,-m_Zoffset[Type])),
        LogicReflector,
        "CatanaReflector",result,false,0);
     
     new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
        LogicCrystal,
        "CatanaCrystal",LogicReflector,false,0);
     
  return result;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Catana::ReadConfiguration(NPL::InputParser parser){
  // CSV config
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("Catana","CSV");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " CSV block found " << endl; 

  vector<string> token = {"Path","Pos","Rshift"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      string path = blocks[i]->GetString("Path");
      double Rshift = blocks[i]->GetDouble("Rshift","micrometer");
      // Reference position of the whole array
      m_Ref = NPS::ConvertVector(blocks[i]->GetTVector3("Pos","mm"));
      ReadCSV(path,Rshift);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
 
  // Type 1
  blocks = parser.GetAllBlocksWithTokenAndValue("Catana","Detector");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  token = {"X","Y","Z","Theta","Phi","ID","Type"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      double X = blocks[i]->GetDouble("X","mm");
      double Y = blocks[i]->GetDouble("Y","mm");
      double Z = blocks[i]->GetDouble("Z","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      int    ID  = blocks[i]->GetInt("ID");
      int    Type =  blocks[i]->GetInt("Type"); 
      AddDetector(X,Y,Z,Theta,Phi,ID,Type);
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
void Catana::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_X.size() ; i++) {
    if(m_Type[i]>5)
      exit(1);
    BuildDetector(m_Type[i]);
    // Reference coordinate given for center of crystal
    G4ThreeVector Det_pos = G4ThreeVector(m_X[i],m_Y[i],m_Z[i]) ;
    // But placed volume is casing which is shifted w/respect to crystal 
    G4ThreeVector Det_dir = Det_pos;
    Det_dir.unit();
    // had to add a 70micron in radius to avoid overlap when using official
    // csv simulation file
    Det_dir.setMag(m_Zoffset[m_Type[i]]+m_Rshift[i]);
    Det_pos+=Det_dir+m_Ref;
    G4RotationMatrix* Rot = new G4RotationMatrix();
    Rot->rotateX(-m_Theta[i]);
    Rot->rotateZ(m_Phi[i]);
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildDetector(m_Type[i]),
          "Catana",world,false,m_ID[i]);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Catana::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Catana")){
    pTree->Branch("Catana", "TCatanaData", &m_Event) ;
  }
  pTree->SetBranchAddress("Catana", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Catana::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_CatanaScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Catana_NS::ResoEnergy);
    if(Energy>Catana_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Catana_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Catana::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_CatanaScorer = CheckScorer("CatanaScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(2);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_CatanaScorer->RegisterPrimitive(Calorimeter);
  m_CatanaScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_CatanaScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Catana::Construct(){
  return  (NPS::VDetector*) new Catana();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Catana{
    public:
      proxy_nps_Catana(){
        NPS::DetectorFactory::getInstance()->AddToken("Catana","Catana");
        NPS::DetectorFactory::getInstance()->AddDetector("Catana",Catana::Construct);
      }
  };

  proxy_nps_Catana p_nps_Catana;
}
