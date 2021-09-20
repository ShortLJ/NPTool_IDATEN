/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: E. Tronchin                                              *
 * contact address: elidiano.tronchin@studenti.unipd.it                      *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Dali simulation                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
using namespace std;

//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
// #ifndef G4UEXTRUDEDSOLID_hh
// #define G4UEXTRUDEDSOLID_hh
// #include "G4USolid.hh"
//  #if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

//#include "G4UExtrudedSolid.hh"
#include "G4TwoVector.hh"
#include "G4TessellatedSolid.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
//#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// NPTool header
#include "Dali.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Dali_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0*MeV;
  const double ResoTime = 0.0*ns; //4.5*ns ;
  const double ResoEnergy = 0.122;  //Relative resolution DeltaE = 0.122*Sqrt(E) 
  /* const double ResoEnergy = 1.36*MeV ; // mean Resolution(FWHM) 1.7% of 80MeV from slides 20170214-SAMURAI34-setup-DALI.pdf if 1.7% of 30MeV = 0.51 MeV // 0.001*MeV ; */
  //const double Radius = 50*mm ; 
  const double Width = 49.76*mm ;
  const double Hight = 84.81*mm ;
  const double Thickness = 164.82*mm ;
  const double LengthPMT = 152.62*mm ;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Dali Specific Method
Dali::Dali(){
  m_Event = new TDaliData() ;
  m_SquareDetector = 0;
  m_DaliScorer = 0;
  Log_Dali_1Volume = 0;
  Log_Al_Cryst_can = 0;
  Log_MgO_Cryst_can =0;
  Log_Crystal = 0;
  Log_Dali_3Volume=0;
 }

Dali::~Dali(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dali::AddDetector(G4ThreeVector POS ){
  // Convert the POS value to R theta Phi as Cylindrical coordinate is easier in G4 
  m_R.push_back(POS.perp());
  m_Alpha.push_back(POS.phi());
  m_Zeta.push_back(POS.y());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dali::AddDetector(double  R, double  Theta, double  Phi){
  double m_r, m_alpha, m_zeta;
  m_r = R*cos(Phi);
  m_alpha = Theta;
  m_zeta = R*sin(Phi);
  m_R.push_back(m_r);
  m_Alpha.push_back(m_alpha);
  m_Zeta.push_back(m_zeta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dali::AddDetector(double  R, double  Alpha, double  Zeta, int Ring){
  m_R.push_back(R);
  m_Alpha.push_back(Alpha);
  m_Zeta.push_back(Zeta);
  m_Ring.push_back(Ring);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Definition Materials MgO and NaI(Tl)

 void Dali::DefinitionMaterials()
 {
  
    G4Element *Tl = new G4Element("Thallium","Tl",81.,204.383*g/mole );

    NaI_Tl = new G4Material("NaI_Tl",3.6667*g/cm3, 2);
    NaI_Tl->AddMaterial(MaterialManager::getInstance()->GetMaterialFromLibrary("NaI"), 99.6*perCent);
    NaI_Tl->AddElement(Tl, 0.4*perCent);
 
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* Dali::BuildSquareDetector(){
  if(!m_SquareDetector){
    
    std::vector<G4TwoVector> polygon;
    polygon.push_back(G4TwoVector(Dali_NS::Hight*0.5, Dali_NS::Width*0.5*3.)  ) ;
    polygon.push_back(G4TwoVector(Dali_NS::Hight*0.5, -Dali_NS::Width*0.5*3.)  ) ;
    polygon.push_back(G4TwoVector(-Dali_NS::Hight*0.5, -Dali_NS::Width*0.5*3.)  ) ;
    polygon.push_back(G4TwoVector(-Dali_NS::Hight*0.5, Dali_NS::Width*0.5*3.)  ) ;
    
    G4Material* Aria = MaterialManager::getInstance()->GetMaterialFromLibrary("Air");
    G4Material* Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* MgO = MaterialManager::getInstance()->GetMaterialFromLibrary("MgO");
    //G4Material* muMetal = MaterialManager::getInstance()->GetMaterialFromLibrary("mumetal");
    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    //G4Material* BoroSili_Glass = MaterialManager::getInstance()->GetMaterialFromLibrary("Borosillicate_Glass");
    
    ///// CONTENAIR VOLUMES 
    G4Box* Dali_3Volume = new G4Box("Dali_3Volume", 
                                 Dali_NS::Hight*0.5,
                                 Dali_NS::Width*0.5*3,
                                 Dali_NS::Thickness*0.5 + Dali_NS::LengthPMT/2.+11.5/2.*mm );        
    Log_Dali_3Volume = new G4LogicalVolume(Dali_3Volume,
                                            Aria,"log_Dali_3Volume",0,0,0);
    G4Box* Dali_1Volume = new G4Box("Dali_1Volume", 
                                      Dali_NS::Hight*0.5,
                                      Dali_NS::Width*0.5, 
                                      Dali_NS::Thickness*0.5 + Dali_NS::LengthPMT/2.+11.5/2.*mm );
    Log_Dali_1Volume = new G4LogicalVolume(Dali_1Volume, 
                                                Aria,"logic_1DaliVolume",0,0,0);
    
    G4Box* Extrudedbox_can  = new G4Box("extrude_box", Dali_NS::Hight*0.5,
                                         Dali_NS::Width*0.5, 
                                         Dali_NS::LengthPMT/2.+11.5/2.*mm);
    AriaExtrude = new G4LogicalVolume(Extrudedbox_can,Aria, "logic_Ariaextrude",0,0,0);


    /////// PMTs objets 
    G4Tubs* AlPMT = new G4Tubs("AlPMT",16.5*mm, 19.5*mm,Dali_NS::LengthPMT/2.,0*deg,360*deg);
    lAlPMT = new G4LogicalVolume(AlPMT, Vacuum ,"lAlPMT",0,0,0);
    
    G4Tubs* MuPMT = new G4Tubs("MuPMT",16.5*mm,20.*mm,Dali_NS::LengthPMT/2.,0*deg,360*deg);
    lMuPMT = new G4LogicalVolume(MuPMT, Vacuum ,"lMuPMT",0,0,0);
    
    G4Box* TopPlatePMT = new G4Box("TopPlatePMT",  Dali_NS::Hight*0.5-1*mm,
                                  Dali_NS::Width*0.5-1*mm, 
                                  11.5/2.*mm );
    lTopPlatePMT = new G4LogicalVolume(TopPlatePMT, Vacuum ,"lTopPlatePMT",0,0,0);
    
    G4Tubs* GlassPMT = new G4Tubs("GlassPMT", 0. , 16.5*mm , 11.5/2.*mm ,0*deg,360*deg);
    lGlassPMT = new G4LogicalVolume(GlassPMT , Vacuum ,"lGlassPMT",0,0,0); // Cyril


    
    //////////////// CRYSTAL part
    G4Box* Al_Cryst_can = new G4Box("Al_Cryst_can", 
                                    Dali_NS::Hight*0.5,
                                    Dali_NS::Width*0.5,  
                                    Dali_NS::Thickness*0.5);
    Log_Al_Cryst_can    = new G4LogicalVolume(Al_Cryst_can,
                                    Al,"log_Al_Cryst_can",0,0,0);

    G4Box* MgO_Cryst_can = new G4Box("MgO_Cryst_can", 
                                    Dali_NS::Hight*0.5-1*mm,
                                    Dali_NS::Width*0.5-1*mm, 
                                    Dali_NS::Thickness*0.5-1*mm);  
    Log_MgO_Cryst_can = new G4LogicalVolume(MgO_Cryst_can,
                                    MgO, "logic_Dali_CanMg0",0,0,0);

    G4Box* Crystal  = new G4Box("Crystal", 
                                Dali_NS::Hight*0.5-2.4*mm,
                                Dali_NS::Width*0.5-2.4*mm, 
                                Dali_NS::Thickness*0.5-2.4*mm);
    Log_Crystal     = new G4LogicalVolume(Crystal,
                                          NaI_Tl,"logic_Dali_Box",0,0,0);

    G4ThreeVector positionnull = G4ThreeVector(0,0,0);

    // PMT  part -
    new G4PVPlacement(0, positionnull,
                      lAlPMT ,"AlPMT",lMuPMT,false,0);
    new G4PVPlacement(0, G4ThreeVector(0,0, -11.5/2.*mm ),
                      lMuPMT ,"MuPMT",AriaExtrude,false,0);
    new G4PVPlacement(0, positionnull,lGlassPMT,"GlassPMT",
                      lTopPlatePMT,false,0);
    new G4PVPlacement(0,  G4ThreeVector(0,0, Dali_NS::LengthPMT/2. ),
                      lTopPlatePMT,"TopPlatePMT",AriaExtrude,false,0);
    new G4PVPlacement(0,  G4ThreeVector(0,0, -Dali_NS::Thickness*0.5 ),
                      AriaExtrude,"PMTVolume",Log_Dali_1Volume,false,0);
    

    // Cryst Part -
    new G4PVPlacement(0,  G4ThreeVector(0,0, 
                      Dali_NS::LengthPMT/2.+11.5/2.*mm ),
                      Log_Al_Cryst_can,
                      "DetectorVolume",
                      Log_Dali_1Volume,false,0);
                      
    new G4PVPlacement(0, positionnull,
                      Log_MgO_Cryst_can,
                      "MgO_Can",
                      Log_Al_Cryst_can,false,0); 
    new G4PVPlacement(0, positionnull,
                      Log_Crystal,
                      "CrystalNaI",
                      Log_MgO_Cryst_can,false,0); 
    new G4PVReplica("DaliArrayElement",
                    Log_Dali_1Volume,
                    Log_Dali_3Volume ,
                    kYAxis,3,Dali_NS::Width,0);
     
    G4VisAttributes* MgO_Color = new G4VisAttributes(G4Colour(1,1,1, .4));
    G4VisAttributes* Al_Color = new G4VisAttributes(G4Colour(0.5,0.5,0.5, .3));
    G4VisAttributes* Crystal_Color = new G4VisAttributes(G4Colour(0, 1, 1));   
    //G4VisAttributes* mumetal_Color = new G4VisAttributes(G4Colour(0, 0.5, 1, .3));   


    Log_Dali_3Volume->SetVisAttributes(G4VisAttributes(G4Colour(1,1,1, 0)));
    Log_Dali_1Volume->SetVisAttributes(G4VisAttributes(G4Colour(1,1,1,0)));
    AriaExtrude->SetVisAttributes(G4VisAttributes(G4Colour(1,1,1,0)));
    
    Log_MgO_Cryst_can->SetVisAttributes(MgO_Color);
    Log_Crystal->SetVisAttributes(Crystal_Color);
    Log_Al_Cryst_can->SetVisAttributes(Al_Color);
    
    lAlPMT->SetVisAttributes(Al_Color);
    lMuPMT->SetVisAttributes(Al_Color);
    lTopPlatePMT->SetVisAttributes(Al_Color);

    Log_Crystal->SetSensitiveDetector(m_DaliScorer);
  }
  return Log_Dali_3Volume;
} // end BuildSquareDetector

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Dali::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Dali");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> cyli = {"R","Alpha","Zeta"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Dali " << i+1 <<  endl;
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);
    }

    else if(blocks[i]->HasTokenList(cyli)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Dali " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Alpha = blocks[i]->GetDouble("Alpha","deg");
      double Zeta = blocks[i]->GetDouble("Zeta","mm");
      int    Ring = blocks[i]->GetInt("Ring");
      AddDetector(R,Alpha,Zeta,Ring);
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
void Dali::ConstructDetector(G4LogicalVolume* world){    

  DefinitionMaterials();
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {
    G4double wX = m_R[i] * cos(m_Alpha[i] ) ; 
    G4double wY = m_R[i] * sin(m_Alpha[i] ) ;
    G4double wZ = m_Zeta[i];
    if(m_Ring[i]==1) wZ = wZ - (Dali_NS::Thickness+11.5*mm+Dali_NS::LengthPMT)/2. + Dali_NS::Thickness/2.;
    else wZ = wZ + (Dali_NS::Thickness+11.5*mm+Dali_NS::LengthPMT)/2. - Dali_NS::Thickness/2.;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;

    G4RotationMatrix* Rot = new G4RotationMatrix();

    Rot->rotateX(180*deg);
    if(m_Ring[i]==1){
      Rot->rotateY(180*deg); Rot->rotateZ(m_Alpha[i]); 
    } else{Rot->rotateZ(m_Alpha[i]);}

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
        BuildSquareDetector(),
        "Dali",world,false,i+1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Dali::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Dali")){
    pTree->Branch("Dali", "TDaliData", &m_Event) ;
  }
  pTree->SetBranchAddress("Dali", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Dali::ReadSensitive(const G4Event* ){
  m_Event->Clear();
  ///////////

  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_DaliScorer->GetPrimitive(0);

  //cout << m_DaliScorer->GetNumberOfPrimitives()<<endl;
  unsigned int size = Scorer->GetMult();
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Dali_NS::ResoEnergy*std::sqrt(Scorer->GetEnergy(i)));
    // cout << Energy << endl;
    if(Energy>Dali_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Dali_NS::ResoTime);
      int ArrayNbr = level[1];
      int DetectinsArrayNbr = level[0]+1;
      int DetectorNbr = (ArrayNbr-1)*3+DetectinsArrayNbr;
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time);
      /* m_Event->SetParticleID(Scorer->GetParticleID(i)); */
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Dali::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;
  vector<int> NestingLevel;
  NestingLevel.push_back(3);
  NestingLevel.push_back(4);

  m_DaliScorer = CheckScorer("DaliScorer",already_exist) ;

  if(already_exist) //Necessary?
    return ;  //Necessary?
  // Otherwise the scorer is initialised
  //  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter", NestingLevel) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_DaliScorer->RegisterPrimitive(Calorimeter);
  m_DaliScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_DaliScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Dali::Construct(){
  return  (NPS::VDetector*) new Dali();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Dali{
    public:
      proxy_nps_Dali(){
        NPS::DetectorFactory::getInstance()->AddToken("Dali","Dali");
        NPS::DetectorFactory::getInstance()->AddDetector("Dali",Dali::Construct);
      }
  };

  proxy_nps_Dali p_nps_Dali;
}
