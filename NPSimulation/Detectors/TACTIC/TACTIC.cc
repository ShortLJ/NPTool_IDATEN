/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Warren Lynch  contact address: Warren.Lynch@york.ac.uk                        *
 *                                                                           *
 * Creation Date  : June 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  TACTIC simulation                             *
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
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolumeStore.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4THitsMap.hh"
#include "G4SDParticleFilter.hh"

// NPTool header
#include "TACTIC.hh"
#include "TACTICScorer.hh"

#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

// For reaction file
#include "BeamReaction.hh"
#include "NPFunction.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace TACTIC_NS{
  const double drift_radius = 50.*mm;
  const double active_length = 251.9*mm;
  const double window_pos = 104.*mm; //from centre of TACTIC from https://elog.triumf.ca/Tactic/Documentation/18                                     
  const double window_radius = 12.*mm; //guess
  const double window_width = 1.5e-03*mm;
  const int NumberOfStrips = 60;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// TACTIC Specific Method
TACTIC::TACTIC(){
  m_Event = new TTACTICData() ;
  m_Scorer = 0;
  m_CylindricalDetector = 0;

  m_ReactionRegion = NULL; 
  // RGB Color + Transparency

  m_VisChamber        = new G4VisAttributes(G4Colour(1., 0.0, 0.0, 0.0));
  m_VisWindows        = new G4VisAttributes(G4Colour(1, 1, 1, 1.0));
  m_VisGas            = new G4VisAttributes(G4Colour(1, 1, 1, 0.1));
  m_VisVacuum         = new G4VisAttributes(G4Colour(1,1,1,1.));

  m_VisChamber->SetForceWireframe(true);
  m_VisVacuum->SetForceWireframe(true);
}

TACTIC::~TACTIC(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TACTIC::AddDetector(G4ThreeVector POS, string  Shape){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TACTIC::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* TACTIC::BuildCylindricalDetector(){
  //  if(!m_CylindricalDetector){
  
    // G4Material* Cu = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
    G4Material* Mylar = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");
    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_Galactic");

    //    std::cout << "Vacuum: " << Vacuum << std::endl;
    
    unsigned const int NumberOfGasMix = m_GasMaterial.size();

    double density=0;
    vector<G4Material*> GasComponent;
    
    for(unsigned int i=0; i<NumberOfGasMix; i++){
      //      if(m_GasMaterial[i] == "CO2") GasComponent.push_back(MaterialManager::getInstance()->GetGasFromLibrary(m_GasMaterial[i], 1.0/bar, m_Temperature));
      if(m_GasMaterial[i] == "P10_gas") {
	//G4Material *P10_gas = new G4Material("P10_gas", 0.00156 * g /cm3, 3, kStateGas, m_Temperature, 1.0/bar); //density for SRIM (1 atm);
	G4Material *P10_gas = new G4Material("P10_gas", 0.00156*g/cm3, 3);
	G4Element *elAr = new G4Element("Argon","Ar",18.,39.948*g/mole);
	G4Element *elC = new G4Element("Carbon","C",6.,12.0107*g/mole);
	G4Element *elH = new G4Element("Hydrogen","H",1.,1.00784*g/mole);
	P10_gas->AddElement(elAr,90);
	P10_gas->AddElement(elC,2);
	P10_gas->AddElement(elH,8);
	GasComponent.push_back(P10_gas);
      }
      else GasComponent.push_back(MaterialManager::getInstance()->GetMaterialFromLibrary(m_GasMaterial[i]));
    }
    
    for(unsigned int i=0; i<NumberOfGasMix; i++){
      density += ((double)m_GasFraction[i]/100)*GasComponent[i]->GetDensity()*(m_Pressure/bar)/(GasComponent[i]->GetPressure()/bar);
      // p2 = p1*(P2/P1) //e.g.  p2 = p1*(0.5/1) p scales with P, T=const 
    }

    G4Material* TACTIC_gas = new G4Material("TACTIC_gas", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);

    for(unsigned int i=0; i<NumberOfGasMix; i++) TACTIC_gas->AddMaterial(GasComponent[i], (double)m_GasFraction[i]/100);

    cout << TACTIC_gas << endl;

    
    if(m_Shape[0] == "Cylindrical") {
    
      G4Tubs* world_volume = new G4Tubs("world_volume",0,TACTIC_NS::drift_radius+1*mm,TACTIC_NS::active_length*0.5+1*mm,0,360*deg);
      G4Tubs* window = new G4Tubs("window",0,TACTIC_NS::window_radius,TACTIC_NS::window_width*0.5,0,360*deg);

      //G4Tubs* anode = new G4Tubs("anode",TACTIC_NS::drift_radius,TACTIC_NS::anode_radius,TACTIC_NS::active_length*0.5,0,360*deg);
      G4Tubs* gas_cathode = new G4Tubs("gas_cathode",0,TACTIC_NS::window_radius,TACTIC_NS::window_pos-TACTIC_NS::window_width*0.5,0,360*deg); //window pos doesn't need halving
      G4Tubs* gas_drift = new G4Tubs("gas_drift",TACTIC_NS::window_radius,TACTIC_NS::drift_radius,TACTIC_NS::active_length*0.5,0,360*deg);
      G4UnionSolid* gas_volume = new G4UnionSolid("gas_volume",gas_cathode,gas_drift);
      
      //G4Tubs* gas_volume = new G4Tubs("gas_volume",0,TACTIC_NS::drift_radius,TACTIC_NS::active_length*0.5,0,360*deg); 
      
      m_CylindricalDetector = new G4LogicalVolume(world_volume, Vacuum, "m_CylindricalDetector_log",0,0,0);
      //anode_log = new G4LogicalVolume(anode, Cu, "anode_log", 0,0,0);
      gas_volume_log = new G4LogicalVolume(gas_volume, TACTIC_gas, "gas_volume_log",0,0,0);
      window_log = new G4LogicalVolume(window, Mylar, "window_log",0,0,0);
      
      new G4PVPlacement(0,G4ThreeVector(0,0,0),gas_volume_log,"gas_volume_phys",m_CylindricalDetector,false,0);
      //new G4PVPlacement(0,G4ThreeVector(0,0,0),anode_log,"anode_phys",m_CylindricalDetector,false,0,true);
      new G4PVPlacement(0,G4ThreeVector(0,0,TACTIC_NS::window_pos),window_log,"window_phys1",m_CylindricalDetector,false,0,true);
      new G4PVPlacement(0,G4ThreeVector(0,0,-TACTIC_NS::window_pos),window_log,"window_phys2",m_CylindricalDetector,false,0,true);

    }
    
    if(m_Shape[0] == "Long_Chamber") {

      G4Box* world_volume = new G4Box("world_volume",34.7,25.,TACTIC_NS::active_length/2.); //Pad width = 34.7 mm (2 rows of pads), LC drift = 50 mm
      m_CylindricalDetector = new G4LogicalVolume(world_volume, Vacuum, "m_CylindricalDetector_log",0,0,0);
      G4Box* gas_volume = new G4Box("gas_volume",34.7,25.,TACTIC_NS::active_length/2.);
      gas_volume_log = new G4LogicalVolume(gas_volume, TACTIC_gas, "gas_volume_log",0,0,0);
      new G4PVPlacement(0,G4ThreeVector(0,0,0),gas_volume_log,"gas_volume_phys",m_CylindricalDetector,false,0);
      
    }

    m_CylindricalDetector->SetVisAttributes(m_VisGas);
    gas_volume_log->SetVisAttributes(m_VisGas);

    G4UserLimits *gas_volume_step = new G4UserLimits();
    G4double maxStep = 0.1*mm;
    gas_volume_step->SetMaxAllowedStep(maxStep);
    gas_volume_log->SetUserLimits(gas_volume_step);    
    gas_volume_log->SetSensitiveDetector(m_Scorer);
        
  return m_CylindricalDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void TACTIC::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("TACTIC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 
  
  vector<string> cart = {"POS","Shape","GasMaterial_1","GasMaterial_2","GasFraction_1","GasFraction_2","Temperature","Pressure","Active"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TACTIC " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      Shape = blocks[i]->GetString("Shape");
      m_GasMaterial.push_back(blocks[i]->GetString("GasMaterial_1"));
      m_GasMaterial.push_back(blocks[i]->GetString("GasMaterial_2"));
      m_GasFraction.push_back(blocks[i]->GetInt("GasFraction_1"));
      m_GasFraction.push_back(blocks[i]->GetInt("GasFraction_2"));
      m_Temperature = blocks[i]->GetDouble("Temperature","kelvin");
      m_Pressure = blocks[i]->GetDouble("Pressure","bar");
      m_Active = blocks[i]->GetString("Active");
      m_p0 = blocks[i]->GetDouble("p0","0");
      m_p1 = blocks[i]->GetDouble("p1","0");
      m_p2 = blocks[i]->GetDouble("p2","0");
      m_p3 = blocks[i]->GetDouble("p3","0");
      //cout << "ACTIVE: " << m_Active << endl;
      AddDetector(Pos,Shape);
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
void TACTIC::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*TACTIC_NS::active_length*0.5;
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
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),BuildCylindricalDetector(),"TACTIC",world,false,i+1);

    if(!m_ReactionRegion){
      G4ProductionCuts* ecut = new G4ProductionCuts();
      ecut->SetProductionCut(1e06,"e-"); //I think lowest is 900 eV for delta electron production, this is 1e06 MeV to cut all electrons produced this way
      m_ReactionRegion= new G4Region("NPSimulationProcess");
      m_ReactionRegion->SetProductionCuts(ecut);
      //m_ReactionRegion->AddRootLogicalVolume(gas_volume_log);
      //m_ReactionRegion->AddRootLogicalVolume(m_CylindricalDetector);

      if(m_Active == "windows") m_ReactionRegion->AddRootLogicalVolume(window_log);
      if(m_Active == "gas") m_ReactionRegion->AddRootLogicalVolume(gas_volume_log);
    }

    G4FastSimulationManager* mng = m_ReactionRegion->GetFastSimulationManager();
    unsigned int size = m_ReactionModel.size();
    for(unsigned int j = 0 ; j < size ; j++){
      mng->RemoveFastSimulationModel(m_ReactionModel[j]);
    }
    m_ReactionModel.clear();
    
    G4VFastSimulationModel* fsm;
    fsm = new NPS::BeamReaction("BeamReaction",m_ReactionRegion);
    m_ReactionModel.push_back(fsm);
        
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void TACTIC::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("TACTIC")){
    pTree->Branch("TACTIC", "TTACTICData", &m_Event) ;
  }
  pTree->SetBranchAddress("TACTIC", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void TACTIC::ReadSensitive(const G4Event* event ){
  m_Event->Clear();

  ofstream file;
  /*
  file.open("signal.dat", std::ios::app);
  file << "Event" << endl;
  file.close();
  */
  file.open("out.dat",std::ios::app);
  
  G4THitsMap<G4double*>* LightHitMap;
  std::map<G4int, G4double**>::iterator Light_itr;
  G4int LightCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TACTICScorer/LightScorer");
  LightHitMap = (G4THitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(LightCollectionID));

  for (Light_itr = LightHitMap->GetMap()->begin(); Light_itr != LightHitMap->GetMap()->end(); Light_itr++) {
    G4double* Info = *(Light_itr->second);
    //file << event->GetEventID() << "\t";
    for(int s = 0; s<13; s++) {
      file << Info[s] << "\t";
    }
    file << Info[13] << endl;
  }
  
  LightHitMap->clear();

  G4THitsMap<G4double*>* HeavyHitMap;
  std::map<G4int, G4double**>::iterator Heavy_itr;
  G4int HeavyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TACTICScorer/HeavyScorer");
  HeavyHitMap = (G4THitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(HeavyCollectionID));

  for (Heavy_itr = HeavyHitMap->GetMap()->begin(); Heavy_itr != HeavyHitMap->GetMap()->end(); Heavy_itr++) {
    G4double* Info = *(Heavy_itr->second);
    //file << event->GetEventID() << "\t";
    for(int s = 0; s<13; s++) {
      file << Info[s] << "\t";
    }
    file << Info[13] << endl;
  }
  
  HeavyHitMap->clear();
  
  G4THitsMap<G4double*>* BeamHitMap;
  std::map<G4int, G4double**>::iterator Beam_itr;
  G4int BeamCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TACTICScorer/BeamScorer");
  BeamHitMap = (G4THitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(BeamCollectionID));
  
  for (Beam_itr = BeamHitMap->GetMap()->begin(); Beam_itr != BeamHitMap->GetMap()->end(); Beam_itr++) {
    G4double* Info = *(Beam_itr->second);
    //file << event->GetEventID() << "\t";
    for(int s = 0; s<13; s++) {
      file << Info[s] << "\t";
    }
    file << Info[13] << endl;
  }

  BeamHitMap->clear();
 
  file.close();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void TACTIC::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_Scorer = CheckScorer("TACTICScorer",already_exist) ;
  if(already_exist) 
    return ;
  
  // Otherwise the scorer is initialised

  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());
  NPL::Reaction m_Reaction;
  m_Reaction.ReadConfigurationFile(input);
  
  if(input.GetAllBlocksWithToken("TwoBodyReaction").size()>0) {
    // NEED SEPERATE SCORERS FOR EACH PARTICLE! OTHERWISE SCORER JUST RETURNS OUTPUT FOR ONE PARTICLE WHEN THERE IS OVERLAP
    G4VPrimitiveScorer* LightScorer = new TACTICScorer::Gas_Scorer("LightScorer",1,TACTIC_NS::active_length,TACTIC_NS::NumberOfStrips,0,m_p0,m_p1,m_p2,m_p3,Shape);
    G4VPrimitiveScorer* HeavyScorer = new TACTICScorer::Gas_Scorer("HeavyScorer",1,TACTIC_NS::active_length,TACTIC_NS::NumberOfStrips,0,m_p0,m_p1,m_p2,m_p3,Shape);
    G4VPrimitiveScorer* BeamScorer = new TACTICScorer::Gas_Scorer("BeamScorer",1,TACTIC_NS::active_length,TACTIC_NS::NumberOfStrips,0,m_p0,m_p1,m_p2,m_p3,Shape);
    G4SDParticleFilter* LightFilter = new G4SDParticleFilter("LightFilter");
    G4SDParticleFilter* HeavyFilter = new G4SDParticleFilter("HeavyFilter");
    G4SDParticleFilter* BeamFilter = new G4SDParticleFilter("BeamFilter");
  
  
    cout << m_Reaction.GetParticle1()->GetName() << "\t" << m_Reaction.GetParticle1()->GetA() << "\t" << m_Reaction.GetParticle1()->GetZ() << endl;
    cout << m_Reaction.GetParticle2()->GetName() << "\t" << m_Reaction.GetParticle2()->GetA() << "\t" << m_Reaction.GetParticle2()->GetZ() << endl;
    cout << m_Reaction.GetParticle3()->GetName() << "\t" << m_Reaction.GetParticle3()->GetA() << "\t" << m_Reaction.GetParticle3()->GetZ() << endl;
    cout << m_Reaction.GetParticle4()->GetName() << "\t" << m_Reaction.GetParticle4()->GetA() << "\t" << m_Reaction.GetParticle4()->GetZ() << endl;
  
    LightFilter->addIon(m_Reaction.GetParticle3()->GetZ(),m_Reaction.GetParticle3()->GetA());
    HeavyFilter->addIon(m_Reaction.GetParticle4()->GetZ(),m_Reaction.GetParticle4()->GetA());
    if(m_Reaction.GetParticle1()->GetZ() == m_Reaction.GetParticle4()->GetZ()) BeamFilter->add("geantino");
    else BeamFilter->addIon(m_Reaction.GetParticle1()->GetZ(),m_Reaction.GetParticle1()->GetA());
    
    LightScorer->SetFilter(LightFilter);
    HeavyScorer->SetFilter(HeavyFilter);
    BeamScorer->SetFilter(BeamFilter);

    m_Scorer->RegisterPrimitive(LightScorer);
    m_Scorer->RegisterPrimitive(HeavyScorer);
    m_Scorer->RegisterPrimitive(BeamScorer);

  }

  if(input.GetAllBlocksWithToken("Isotropic").size()>0) { //For alpha, proton or neutron source (obviously neutron wont do anything so this is spare).

    G4VPrimitiveScorer* LightScorer = new TACTICScorer::Gas_Scorer("LightScorer",1,TACTIC_NS::active_length,TACTIC_NS::NumberOfStrips,0,m_p0,m_p1,m_p2,m_p3,Shape);
    G4VPrimitiveScorer* HeavyScorer = new TACTICScorer::Gas_Scorer("HeavyScorer",1,TACTIC_NS::active_length,TACTIC_NS::NumberOfStrips,0,m_p0,m_p1,m_p2,m_p3,Shape);
    G4VPrimitiveScorer* BeamScorer = new TACTICScorer::Gas_Scorer("BeamScorer",1,TACTIC_NS::active_length,TACTIC_NS::NumberOfStrips,0,m_p0,m_p1,m_p2,m_p3,Shape);
    G4SDParticleFilter* LightFilter = new G4SDParticleFilter("LightFilter","alpha");
    G4SDParticleFilter* HeavyFilter = new G4SDParticleFilter("HeavyFilter","proton");
    G4SDParticleFilter* BeamFilter = new G4SDParticleFilter("BeamFilter","neutron");
    LightScorer->SetFilter(LightFilter);
    HeavyScorer->SetFilter(HeavyFilter);
    BeamScorer->SetFilter(BeamFilter);
    
    m_Scorer->RegisterPrimitive(LightScorer);
    m_Scorer->RegisterPrimitive(HeavyScorer);
    m_Scorer->RegisterPrimitive(BeamScorer);

  }
    
    
  G4SDManager::GetSDMpointer()->AddNewDetector(m_Scorer);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* TACTIC::Construct(){
  return  (NPS::VDetector*) new TACTIC();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_TACTIC{
    public:
      proxy_nps_TACTIC(){
        NPS::DetectorFactory::getInstance()->AddToken("TACTIC","TACTIC");
        NPS::DetectorFactory::getInstance()->AddDetector("TACTIC",TACTIC::Construct);
      }
  };

  proxy_nps_TACTIC p_nps_TACTIC;
}
