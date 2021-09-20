/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: E. Tronchin  contact address: tronchin@lpccaen.in2p3.fr  *
 * Updated by C. Lenain          contact address: lenain@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2018                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Minos simulation                                    *
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
#include "G4Trd.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh" 
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"

// G4 Field
#include "G4FieldManager.hh"
#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4MaterialPropertiesTable.hh"

// NPTool header
#include "Minos.hh"
#include "InteractionScorers.hh"
/* #include "TPCScorers.hh" */
#include "CylinderTPCScorers.hh"

#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

// ROOT
#include "TH1F.h"
#include "TF1.h"
#include "RootOutput.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Minos_NS{

  // TPC
  const G4double  ChamberInnerRadius      = 37.*mm; 
  /* const G4double  ChamberInnerRadius      = 29*mm; */ //big TPC
  //const G4double  ChamberThickness        = 2.*mm; 
  const G4double  ChamberLength           = 300.*mm;
  const G4double  KaptonThickness         = 0.125*mm; 
  const G4double  RohacellThickness  = 2.*mm;
  const G4double  TPCRadiusExt            = 91.525*mm;
  /* const G4double  TPCRadiusExt            = 150*mm; //big TPC */

  // MINOS
  const G4double  TargetRadius            =  28.*mm; 
  const G4double  WindowThickness         = 0.150*mm;

}

using namespace Minos_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Minos Specific Method
Minos::Minos(){
  m_Event = new TMinosData() ;
  m_MinosPadScorer = 0;
  m_ReactionRegion=NULL;

  // RGB Color + Transparency
  m_VisTarget= new G4VisAttributes(G4Colour(0.6,1.,1., .4));
  m_VissimpleBox= new G4VisAttributes(G4Colour(0,1,0,.6));
  m_VisTPC= new G4VisAttributes(G4Colour(1.,0.5,0.6,0.3));
  m_VisRohacell= new G4VisAttributes(G4Colour(1.,1.,1., .8));
  m_VisKapton = new G4VisAttributes(G4Colour(1.,1.,0.6,0.4));
  m_VisTargetCell = new G4VisAttributes(G4Colour(0,0,1, .4));
  m_VisTargetCell->SetForceSolid(true);
  m_VisOuterKapton = new G4VisAttributes(G4Colour(1.,1.,0.6,0.8));

  Raw_Signal = new TH1F("raw_Signal","raw_Signal",512,0,512); 
  Elec_Signal = new TH1F("Elec_Signal","Elec_Signal",512,0,512);

  /*  Raw_Signal = new TH1F;*/
  /*  Elec_Signal = new TH1F;*/


  fa1 = new TF1("fa1","abs((x>[1]&&x<512)*([0]*exp(-3.*(x-[1])/[2]) * sin((x-[1])/[2]) * pow((x-[1])/[2],3))+[3])",0,1000);
  fa1->SetNpx(512);

  solidTarget=0;   
  logicTarget=0;   
  solidChamber=0;  
  logicChamber=0;  
  solidTPC=0; 
  logicTPC=0; 
  solidWindow0=0; 
  logicWindow0=0; 
  solidRohacell=0;   
  logicRohacell=0;   
  solidKapton=0;   
  logicKapton=0;   
}

Minos::~Minos(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



/* void Minos::AddDetector(G4ThreeVector POS, double LengthOfTarget, int PresenceOfMinos){ */
void Minos::AddDetector(G4ThreeVector POS, double LengthOfTarget, G4String MaterialOfTarget,G4String MaterialOfCell, int PresenceOfMinos){
  m_POS.push_back(POS);
  m_TargetLength.push_back(LengthOfTarget);
  m_TargetMaterial.push_back(MaterialOfTarget);
  m_CellMaterial.push_back(MaterialOfCell);
  m_TPCOnly.push_back(PresenceOfMinos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// with the interactive expansion / contraction geometry system of the
// vis/OpenInventor driver :

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Minos::BuildTarget(){
  if(!logicTarget){
    //                               
    // Target
    //  
    solidTarget = new G4Tubs("Target",		//its name
        0.,TargetRadius,TargetLength/2.,0,360.);//size
    logicTarget = new G4LogicalVolume(solidTarget,	//its solid
        TargetMaterial,	//its material
        "Target");	//its name    
    logicTarget->SetVisAttributes(m_VisTarget);
  }
  return logicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Minos::BuildTPC(){
  if(!logicTPC){
    //                               
    // TPC
    //
    solidTPC = new G4Tubs("TPC",
        ChamberInnerRadius ,TPCRadiusExt,ChamberLength/2.,0,360.); 
    logicTPC = new G4LogicalVolume(solidTPC,    //its solid
        TPCMaterial, //its material
        "TPC"); //name
    logicTPC->SetVisAttributes(m_VisTPC);
  }
  return logicTPC;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/* G4LogicalVolume* Minos::BuildChamber(){ */
/*   if(!logicChamber){ */
/*     solidChamber = new G4Tubs("Chamber",			//its name */
/*         ChamberInnerRadius,ChamberInnerRadius+ChamberThickness,ChamberLength/2.,0,360.); //size */
/*     logicChamber = new G4LogicalVolume(solidChamber,	//its solid */
/*         ChamberMaterial,	//its material */
/*         "Chamber");	//its name */                               
/*     m_VissimpleBox->SetVisibility(true); */
/*     logicChamber->SetVisAttributes(m_VissimpleBox); */
/*   } */
/*   return logicChamber; */
/* } */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Minos::BuildInnerRohacell(){
  if(!logicRohacell){
    solidRohacell = new G4Tubs("InnerRohacell",			//its name
        ChamberInnerRadius ,ChamberInnerRadius + RohacellThickness ,ChamberLength/2.,0,360.); //size
    logicRohacell = new G4LogicalVolume(solidRohacell,	//its solid
        RohacellMaterial,	//its material
        "InnerRohacell");	//its name
    logicRohacell->SetVisAttributes(m_VisRohacell);
  }
  return logicRohacell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4LogicalVolume* Minos::BuildOuterRohacell(){
  if(logicRohacell){
    solidRohacell = new G4Tubs("OuterRohacell",			//its name
        TPCRadiusExt-RohacellThickness, TPCRadiusExt ,ChamberLength/2.,0,360.); //size
    logicRohacell = new G4LogicalVolume(solidRohacell,	//its solid
        RohacellMaterial,	//its material
        "OuterRohacell");	//its name
    logicRohacell->SetVisAttributes(m_VisRohacell);
  }
  return logicRohacell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Minos::BuildKapton(){
  if(!logicKapton){
    solidKapton = new G4Tubs("Kapton",			//its name
        ChamberInnerRadius+RohacellThickness ,ChamberInnerRadius+RohacellThickness+KaptonThickness,ChamberLength/2.,0,360.); //size
    logicKapton = new G4LogicalVolume(solidKapton,	//its solid
        KaptonMaterial,	//its material
        "Kapton");	//its name
    logicKapton->SetVisAttributes(m_VisKapton);
  }

  return logicKapton;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Minos::BuildOuterKapton(){
  if(logicKapton){
    solidKapton = new G4Tubs("Kapton",			//its name
        TPCRadiusExt-RohacellThickness-KaptonThickness, TPCRadiusExt-RohacellThickness,ChamberLength/2.,0,360.); //size
    logicKapton = new G4LogicalVolume(solidKapton,	//its solid
        KaptonMaterial,	//its material
        "Kapton");	//its name
    logicKapton->SetVisAttributes(m_VisOuterKapton);
  }
  return logicKapton;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//                               
// windows
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Minos::BuildWindow0(){
  if(!logicWindow0){
    solidWindow0 = new G4Tubs("WindowTube",		//its name
        0.,TargetRadius+WindowThickness,TargetLength/2.+WindowThickness,0,360.);  
    logicWindow0 = new G4LogicalVolume(solidWindow0,    //its solid
        WindowMaterial, //its material
        "WindowTube"); //name
    logicWindow0->SetVisAttributes(m_VisTargetCell);
  }
  return logicWindow0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Minos::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Minos");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> simu = {"TPCOnly"};
  vector<string> token= {"XML","Position","TargetMaterial","TargetLength","CellMaterial","TimeBin","ShapingTime","Baseline","Sampling","ZOffset"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Minos " << i+1 <<  endl;
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("Position","mm"));
      TargetLength = blocks[i]->GetDouble("TargetLength","mm");
      G4String TargetMaterialname = blocks[i]->GetString("TargetMaterial");
      G4String CellMaterial = blocks[i]->GetString("CellMaterial");
      m_ShapingTime = blocks[i]->GetDouble("ShapingTime","ns")/ns;   
      m_TimeBin = blocks[i]->GetDouble("TimeBin","ns")/ns;   
      m_Sampling= blocks[i]->GetInt("Sampling");   
      m_Baseline= blocks[i]->GetInt("BaseLine");   
      m_ZOffset = blocks[i]->GetDouble("ZOffset","mm");   
      string xmlpath = blocks[i]->GetString("XML");
      NPL::XmlParser xml;
      xml.LoadFile(xmlpath);
      ReadXML(xml);

      TPCOnly=1;
      if(blocks[i]->HasTokenList(simu))
        TPCOnly = blocks[i]->GetInt("TPCOnly");
      AddDetector(Pos,TargetLength,TargetMaterialname, CellMaterial, TPCOnly);
      
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
void Minos::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_POS.size() ; i++) {
    TargetLength = m_TargetLength[i]; 
    TPCOnly = m_TPCOnly[i];

    /* TargetMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"); */ 
    TargetMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("LH2"); 
    WindowMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");
    ChamberMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al"); 
    KaptonMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Kapton");
    RohacellMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Rohacell");
    TPCMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("mixMINOS");

    /* TargetMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"); */ 
    /* WindowMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"); */
    /* ChamberMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"); */ 
    /* KaptonMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"); */
    /* RohacellMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"); */
    /* TPCMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("mixMINOS"); */

    ///////// Drift properties
    G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();      
    MPT->AddConstProperty("DE_PAIRENERGY",30*eV);
    MPT->AddConstProperty("DE_ABSLENGTH",10*pc); 
    MPT->AddConstProperty("DE_DRIFTSPEED",3.475*cm/microsecond);
    MPT->AddConstProperty("DE_TRANSVERSALSPREAD",7e-5*mm2/ns);
    MPT->AddConstProperty("DE_LONGITUDINALSPREAD",7e-5*mm2/ns);

    /* MPT->AddConstProperty("DE_TRANSVERSALSPREAD",0*mm2/ns); */
    /* MPT->AddConstProperty("DE_LONGITUDINALSPREAD",0*mm2/ns); */

    TPCMaterial->SetMaterialPropertiesTable(MPT);

    G4ThreeVector Det_pos = m_POS[i] ;

    double MinosX = Det_pos.x();    
    double MinosY= Det_pos.y();    
    double MinosZ = Det_pos.z() + m_TargetLength[i]/2. ;    

    new G4PVPlacement(0,//its name
        G4ThreeVector(0,0,+ChamberLength/2 ),	
        /* G4ThreeVector(wX,wY, wZ + ChamberLength - m_TargetLength[i]-WindowThickness*2. - 5*mm ),	// Z positioning putting TPC beginn and Target beginning w/ difference of 5mm */ 
        BuildTPC(),	//its logical volume
        "TPC",	//its name
        world,	//its mother  volume
        false,		//no boolean operation
        0);		//copy number

    //////// ELECTRIC FIELD

    G4ElectricField* field = new G4UniformElectricField(G4ThreeVector(0.0,0.0,200*volt/cm));
    // Create an equation of motion for this field
    G4EqMagElectricField*  Equation = new G4EqMagElectricField(field); 
    G4MagIntegratorStepper* Stepper = new G4ClassicalRK4( Equation, 8 );
    // Get the global field manager 
    G4FieldManager* FieldManager= new G4FieldManager();
    // Set this field to the global field manager 
    FieldManager->SetDetectorField(field );
    BuildTPC()->SetFieldManager(FieldManager,true);

    G4MagInt_Driver* IntgrDriver = new G4MagInt_Driver(0.1*mm, 
        Stepper, 
        Stepper->GetNumberOfVariables() );

    G4ChordFinder* ChordFinder = new G4ChordFinder(IntgrDriver);
    FieldManager->SetChordFinder( ChordFinder );

    G4Material* Cu
      = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");

    ///////// CONSTRUCT PADS using parametrized volumes or replica (uncomment the correct section) 

    //// Using REPLICA

    // Construct Pads ring by ring
    for (int RingNbr = 0; RingNbr < 18; RingNbr++){ 

      int PadsPerRing[18]={144,152,156,164,172,176,184,192,196,204,212,216,224,228,236,244,248,256};  
      /* G4double InnerRadius =  (45.2+RingNbr*2.1+0.02)*mm;// 0.02 mm between each ring */
      G4double InnerRadius =  (45.75+RingNbr*2.1)*mm;
      G4double OuterRadius =  (47.75+RingNbr*2.1)*mm;

      //Volume where are placed replica pads
      G4VSolid* AnodeRing = new G4Tubs("ring",InnerRadius,OuterRadius,
          0.01*mm,0.,360*deg);
      G4LogicalVolume * AnodeRing_log = new G4LogicalVolume(AnodeRing, 
          /* mix, "ringL", 0, 0, 0); */
                      TPCMaterial, "ringL", 0, 0, 0);

      {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.8, 0.4, 0.,0.4));
        AnodeRing_log->SetVisAttributes(atb);}

      new G4PVPlacement(0,G4ThreeVector(0,0,-ChamberLength/2+0.01*mm),
          AnodeRing_log,"ring", logicTPC, false, RingNbr);

      G4double Pad_dPhi = (360./(PadsPerRing[RingNbr]+1))*deg; // longitudinal component of Pad
      G4double Pad_shift = (360./PadsPerRing[RingNbr])*deg; // dPhi between each Pads

      G4VSolid* Pad = new G4Tubs("div_ring", InnerRadius, OuterRadius,
          0.01*mm,0, Pad_dPhi);

      G4LogicalVolume* Pad_log = new G4LogicalVolume(Pad,
          Cu,"div_ringL",0,0,0);

      {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.8, 0.4, 0.,0.4));
      Pad_log->SetVisAttributes(atb);}
      Pad_log->SetSensitiveDetector(m_MinosPadScorer);

      new G4PVReplica("div_ring_phys", Pad_log, 
          AnodeRing_log, kPhi, PadsPerRing[RingNbr],Pad_shift ,0); 
    }

    //////////////////////////////////////

    new G4PVPlacement(0,		//its name
        G4ThreeVector(0,0,0),	//at (0,0,0)
        BuildInnerRohacell(),	//its logical volume
        "InnerRohacell",	//its name
        logicTPC,	//its mother  volume
        false,		//no boolean operation
        0);		//copy number

    new G4PVPlacement(0,		//its name
        G4ThreeVector(0,0,0),	//at (0,0,0)
        BuildOuterRohacell(),	//its logical volume
        "OuterRohacell",	//its name
        logicTPC,	//its mother  volume
        /* logicTPC,	//its mother  volume */
        false,		//no boolean operation
        0);		//copy number

    new G4PVPlacement(0, G4ThreeVector(0,0,0),	//at (0,0,0)
        BuildKapton(),	//its logical volume
        "InnerKapton",	//its name
        logicTPC,	//its mother  volume
        false,		//no boolean operation
        0);		//copy number

    new G4PVPlacement(0,		//its name
        G4ThreeVector(0,0,0),	//at (0,0,0)
        BuildOuterKapton(),	//its logical volume
        "OuterKapton",	//its name
        logicTPC,	//its mother  volume
        false,		//no boolean operation
        0);		//copy number

    if(TPCOnly!=1){   
      new G4PVPlacement(0,		//its name
          G4ThreeVector(MinosX,MinosY,MinosZ),	//at (0,0,0)
          BuildWindow0(),	//its logical volume
          "WindowTube",	//its name
          world,	//its mother  volume
          false,		//no boolean operation
          0);

      new G4PVPlacement(0,//no rotation
          G4ThreeVector(0,0,0/*m_TargetLength[i]*/),	//at (0,0,0)
          BuildTarget(),	//its logical volume
          "Target",	//its name
          logicWindow0,	//its mother  volume
          false,		//no boolean operation
          0);		//copy number

      if(!m_ReactionRegion){
        m_ReactionRegion= new G4Region("NPSimulationProcess");
        m_ReactionRegion -> AddRootLogicalVolume(logicTarget);
        m_ReactionRegion->SetUserLimits(new G4UserLimits(1.2*mm)); //???
      }

      G4FastSimulationManager* mng = m_ReactionRegion->GetFastSimulationManager();
      unsigned int size = m_ReactionModel.size();
      for(unsigned int o = 0 ; o < size ; o++){
        mng->RemoveFastSimulationModel(m_ReactionModel[o]);
      }
      m_ReactionModel.clear();
      G4VFastSimulationModel* fsm;
      fsm = new NPS::BeamReaction("BeamReaction",m_ReactionRegion);
      ((NPS::BeamReaction*) fsm)->SetStepSize(1.2*mm);
      m_ReactionModel.push_back(fsm);
      fsm = new NPS::Decay("Decay",m_ReactionRegion);
      m_ReactionModel.push_back(fsm);
    }
  }

  // Visualization attributes
  world->SetVisAttributes (G4VisAttributes::Invisible);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Minos::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Minos")){
    pTree->Branch("Minos", "TMinosData", &m_Event) ;
  }
  pTree->SetBranchAddress("Minos", &m_Event) ;
}

//....oooOO0OOooo.;
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAction
void Minos::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // DriftElectron scorer
  CylinderTPCScorers::PS_TPCAnode* Scorer2= (CylinderTPCScorers::PS_TPCAnode*) m_MinosPadScorer->GetPrimitive(0);
  unsigned int size2 = Scorer2->GetMult();
  for(unsigned int i = 0 ; i < size2 ; i++){
    int Pad = FindPadID(Scorer2->GetPad(i),Scorer2->GetX(i),Scorer2->GetY(i));
    SimulateGainAndDigitizer(Scorer2->GetT(i), Scorer2->GetQ(i),T,Q);

    m_Event->SetPad(Pad, Scorer2->GetQ(i)->size(),&T,&Q);
  }

}

void Minos::SimulateGainAndDigitizer(vector<double>* rawT, vector<double>* rawQ,vector<int>& T, vector<int>& Q){
  T.clear();
  Q.clear();

  Raw_Signal->Reset();
  Elec_Signal->Reset();

  // Add white noise; 
  /* for( int bin = 0 ; bin < Elec_Signal->GetNbinsX() ; bin++) */
  /*   Elec_Signal->SetBinContent(bin,20+(10-gRandom->Uniform(20)*20)); */   

  // Reallocate electrons from each pack
  /* if (rawQ.sze()==1) { */
  /*     for ( int j = 0; j < rawQ[0]; j++){ */
  /*       Raw_Signal->Fill(rawT[0]/30.); */
  /*     } */
  /* } */
  /* else{ */
  /*   for ( int j = 1; j < rawQ.size()-1; j++){ */
  /*     time = (rawT[j+1]-rawT[j-1])/rawQ[j]; */    
  /*     for ( int k = -rawQ[0]/2; k < (rawQ[0]/2); k++){ */
  /*       Raw_Signal->Fill((rawT[j]+k*time)/30.); */
  /*     } */
  /*   } */
  /* } */

  for ( unsigned int j = 0; j < rawQ->size(); j++){
    for ( int i = 0 ; i < (*rawQ)[j]; i++ ){ 
      Raw_Signal->Fill(((*rawT)[j])/m_TimeBin);
    }
  }

  // Application of the electronic reponse function
  for ( int x=0; x <  Raw_Signal->GetNbinsX(); x++){
    if(Raw_Signal->GetBinContent(x) > 0.5){
      start = Raw_Signal->GetBinCenter(x);
      end = Raw_Signal->GetBinCenter(x)+512;  
      // DriftTime = Raw_Signal->GetBinCenter(x);
      fa1->SetRange(start, end);
      fa1->SetParameter(0,1500);
      fa1->SetParameter(1,start);
      fa1->SetParameter(2,m_ShapingTime/(log(2)*m_TimeBin));
      fa1->SetParameter(3,0);

      for (int p=0; p < Raw_Signal->GetBinContent(x)*1500; p++)
        Elec_Signal->Fill(fa1->GetRandom(start,end));
    }  
  }

  for ( int bin=0; bin < Elec_Signal->GetNbinsX(); bin++){
      if(Elec_Signal->GetBinContent(bin)){
        Q.push_back(Elec_Signal->GetBinContent(bin)+m_Baseline);
        T.push_back(Elec_Signal->GetBinCenter(bin));
      }
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Minos::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;
  bool already_exist2 = false;
  bool already_exist3 = false;

  m_MinosPadScorer = CheckScorer("MinosPadScorer",already_exist3) ;

  if(already_exist && already_exist2 && already_exist3 ) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);

  G4VPrimitiveScorer* DriftElectronMinosTPCScorer= new CylinderTPCScorers::PS_TPCAnode("DriftElectronsScore",level, 0) ;

  //and register it to the multifunctionnal detector
  m_MinosPadScorer->RegisterPrimitive(DriftElectronMinosTPCScorer);

  G4VPrimitiveScorer* PadScorer= new CylinderTPCScorers::PS_TPCAnode("MinosTPCAnode",level, 0);
  m_MinosPadScorer->RegisterPrimitive(PadScorer);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_MinosPadScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Minos::ReadXML(NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("MINOS");  
  unsigned int size = b.size();
  for(unsigned int i = 0 ; i < size ; i++){
    unsigned short ID = b[i]->AsInt("ID"); 
    m_XY[ID] = std::make_pair(b[i]->AsDouble("X"),b[i]->AsDouble("Y"));  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
unsigned int Minos::FindPadID(unsigned int G4ID, double X,double Y){
 // if no XML is provided, do nothing
 if(m_XY.size()==0)
   return G4ID;
 // The pad is already identified
 if(m_ID.find(G4ID)!=m_ID.end())
   return m_ID[G4ID];

 // look for the closest pad
 else{
  double d = 1e6;
  double id=0;
  for(auto it = m_XY.begin();it!=m_XY.end();it++){
    double dd = sqrt((it->second.first-X)*(it->second.first-X)+(it->second.second-Y)*(it->second.second-Y));
    if(dd<d){
      d=dd;
      id=it->first;
    } 
  } 
  //cout << G4ID << " " << id << endl;
  m_ID[G4ID]=id;
  return id;
 }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Minos::Construct(){
  return  (NPS::VDetector*) new Minos();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Minos{
    public:
      proxy_nps_Minos(){
        NPS::DetectorFactory::getInstance()->AddToken("Minos","Minos");
        NPS::DetectorFactory::getInstance()->AddDetector("Minos",Minos::Construct);
      }
  };

  proxy_nps_Minos p_nps_Minos;
}
