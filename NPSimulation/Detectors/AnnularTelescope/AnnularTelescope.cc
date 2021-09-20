/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Greg Christian  contact address: gchristian@tamu.edu     *
 *                                                                           *
 * Creation Date  : March 2018                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  AnnularTelescope simulation                         *
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
#include "AnnularTelescope.hh"
#include "CalorimeterScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace AnnularTelescope_Utils;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace AnnularTelescope_NS{
// Energy and time Resolution
// const double EnergyThreshold = 0.*MeV;
const double CsIResoTime = 10*ns ; // ??? (doesn't really matter)
const double CsIResoEnergy = 0.01; // percent, from Wilton
const double SiResoTime = 1*ns ;  // ??? (doesn't really matter)
const double SiResoEnergy = 0.0149; // absolute, from AnnularS1.hh
// const double Radius = 50*mm ; // not used 
// const double Width = 100*mm ; // not used
// const double Thickness = 10*mm ; // from Wilton
// const string Material = "CsI";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// AnnularTelescope Specific Method
AnnularTelescope::AnnularTelescope(){
  m_Event = new TAnnularTelescopeData() ;
  m_AnnularTelescopeScorer_CsI= 0;
  m_AnnularTelescopeScorer_Si = 0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));   
  m_VisCylinder = new G4VisAttributes(G4Colour(0, 0, 1, 0.5));   
}

AnnularTelescope::~AnnularTelescope(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnnularTelescope::AddDetector(AnnularTelescope_Utils::Geometry& geo){
	m_Geo.push_back(geo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::pair<G4LogicalVolume*, G4LogicalVolume*>
AnnularTelescope::BuildDetector(unsigned short i){
																														
	if(!(i<m_Geo.size())){
		G4cerr << "ERROR (AnnularTelescope::BuildDetector):"
					 << "Invalid index: " << i << endl;
		exit(1);
	}

	std::pair<G4LogicalVolume*, G4LogicalVolume*> logic; // output: <Si, CsI>
	//
	// Si Detector
	{
		std::stringstream name;
		name << "AnnularTelescope_Si_" << i+1;
		G4Tubs* tub = new G4Tubs(name.str().c_str(),
														 m_Geo[i].R_min, m_Geo[i].R_max,
														 m_Geo[i].SiThickness*0.5,
														 0, 360);

		G4Material* DetectorMaterial = 
			MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
		logic.first = new G4LogicalVolume(
			tub,DetectorMaterial,"logic_AnnularTelescope_Si_tub",0,0,0);
		logic.first->SetVisAttributes(m_VisSquare);
		logic.first->SetSensitiveDetector(m_AnnularTelescopeScorer_Si);
	}
	//
	// CsI Detector
	{
		std::stringstream name;
		name << "AnnularTelescope_CsI_" << i+1;
		G4Tubs* tub = new G4Tubs(name.str().c_str(),
														 m_Geo[i].R_min, m_Geo[i].R_max,
														 m_Geo[i].CsIThickness*0.5,
														 0, 360);

		G4Material* DetectorMaterial = 
			MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
		logic.second = new G4LogicalVolume(
			tub,DetectorMaterial,"logic_AnnularTelescope_CsI_tub",0,0,0);
		logic.second->SetVisAttributes(m_VisSquare);
		logic.second->SetSensitiveDetector(m_AnnularTelescopeScorer_CsI);
	}
	
  return logic; // <Si, CsI>
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void AnnularTelescope::ReadConfiguration(NPL::InputParser parser){
	vector<Geometry> geometries =
		AnnularTelescope_Utils::ReadConfiguration(parser);
	for(auto& geo : geometries){ 	AddDetector(geo);  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void AnnularTelescope::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_Geo.size() ; i++) {

		auto logic = BuildDetector(i); // <Si, CsI>
		//
		// Si
		// Put FACE of detector at m_Geo[i].Z
		double Z_Si = m_Geo[i].Z + m_Geo[i].SiThickness*0.5;
		G4ThreeVector Si_pos (0, 0, Z_Si);		
		new G4PVPlacement(
			0, Si_pos, logic.first,"AnnularTelescope",world,false,i+1 );
		//
		// CsI
		// Put FACE of detector at m_Geo[i].Z + Si thickness
		double Z_CsI = 
			m_Geo[i].Z + m_Geo[i].SiThickness + m_Geo[i].CsIThickness*0.5;
		G4ThreeVector CsI_pos (0, 0, Z_CsI);
		new G4PVPlacement(
			0, CsI_pos, logic.second,"AnnularTelescope",world,false,i+1 );		
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void AnnularTelescope::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("AnnularTelescope")){
    pTree->Branch("AnnularTelescope", "TAnnularTelescopeData", &m_Event) ;
  }
  pTree->SetBranchAddress("AnnularTelescope", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read scorer
// Helper for ReadSensitive()
void AnnularTelescope::ReadScorer(
	const G4Event* , const string& scorerName,
	std::vector<HitInfo_t>& hits) {
  CalorimeterScorers::PS_Calorimeter* Scorer=0;
   if(scorerName =="AnnularTelescopeScorer_CsI/Calorimeter_CsI" )
    Scorer= (CalorimeterScorers::PS_Calorimeter*) m_AnnularTelescopeScorer_CsI->GetPrimitive(0);
   else if (scorerName == "AnnularTelescopeScorer_Si/Calorimeter_Si")
    Scorer= (CalorimeterScorers::PS_Calorimeter*) m_AnnularTelescopeScorer_Si->GetPrimitive(0);

   unsigned int size = Scorer->GetMult(); 
   for(unsigned int i = 0 ; i < size ; i++){
        vector<unsigned int> level = Scorer->GetLevel(i); 
		hits.push_back(HitInfo_t());
		hits.back().detector = level[0];
		hits.back().energy = Scorer->GetEnergy(i);
		hits.back().time = Scorer->GetTime(i);
	}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void AnnularTelescope::ReadSensitive(const G4Event* event){
  m_Event->Clear();

	// CSI
	std::vector<HitInfo_t> csiHits;
	ReadScorer(event, "AnnularTelescopeScorer_CsI/Calorimeter_CsI", csiHits);
	for(const auto& hit : csiHits) {
		FillCsIData(hit.detector, hit.energy, hit.time, G4ThreeVector(hit.x,hit.y,hit.z));
	}

	// SI
	std::vector<HitInfo_t> siHits;
	ReadScorer(event, "AnnularTelescopeScorer_Si/Calorimeter_Si", siHits);
	for(const auto& hit : siHits) {
		FillSiData(hit.detector, hit.energy, hit.time, G4ThreeVector(hit.x,hit.y,hit.z));
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void AnnularTelescope::FillCsIData(
	int detector_number, double energy, double time, const G4ThreeVector& pos){
	//
	// Figure out wedge
	int wedge_number = 0;
	if(size_t(detector_number - 1) >= m_Geo.size()) {
		std::cerr << "ERROR: AnnularTelescope::FillCsIData: "
							<< "invalid detector number: " << detector_number
							<< ", tried to access index: " << detector_number-1
							<< " in length " << m_Geo.size() << " array\n";
		exit(1);
	}
	const Geometry& geo = m_Geo.at(detector_number - 1);
	double pitch = geo.CsI_Wedge_Angle_Pitch;
	for(const double& phi : geo.CsI_Wedge_Phi_Angle){
		if(pos.phi() >= phi - pitch/2. && pos.phi() < phi + pitch/2.) {
			// match!
			break;
		}
		++wedge_number;
	}
	// start counting at 1, not 0
	++wedge_number;
	//
	// Add resolutions
	double eres = AnnularTelescope_NS::CsIResoEnergy*energy; // percent
	double energy_res = -1000;
	while(energy_res < 0) { // make sure it's a positive number
		energy_res = RandGauss::shoot(energy, eres);
	}
	double tres = AnnularTelescope_NS::CsIResoTime; // absolute
	double time_res = RandGauss::shoot(time, tres);
	//
	// Set CsI energy and time
	m_Event->SetCsIEnergy(detector_number, wedge_number, energy_res);
	m_Event->SetCsITime(detector_number, wedge_number, time_res);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void AnnularTelescope::FillSiData(
	int detector_number, double energy, double time, const G4ThreeVector& pos){
	//
	// Figure out strips
	// Phi
	if(size_t(detector_number - 1) >= m_Geo.size()) {
		std::cerr << "ERROR: AnnularTelescope::FillSiData: "
							<< "invalid detector number: " << detector_number
							<< ", tried to access index: " << detector_number-1
							<< " in length " << m_Geo.size() << " array\n";
		exit(1);
	}
	const Geometry& geo = m_Geo.at(detector_number - 1);
	size_t phi_strip = 1;
	{
		double pitch = geo.Si_Phi_Angle_Pitch;
		for(const double& phi : geo.Si_Strip_Phi_Angle){
			if(pos.phi() >= phi - pitch/2. && pos.phi() < phi + pitch/2.) {
				// match!
				break;
			}
			++phi_strip;
		}
		if(phi_strip == geo.Si_Strip_Phi_Angle.size() + 1) {
			std::cerr << "\nWARNING: Phi strip number not found: "
								<< "Angle of hit [deg]: " << pos.phi()/deg << "\n";
			std::cerr << "Strip angles:\n";
			for(const double& R : geo.Si_Strip_Phi_Angle){
				std::cerr << "   " << (R-pitch/2.)/deg << " - " << (R+pitch/2.)/deg << "\n";
			}
			std::cerr << "---------\n";
		}
	}
	// Theta
	size_t theta_strip = 1;
	{
		double pitch = geo.Si_Theta_Radius_Pitch;
		for(const double& R : geo.Si_Strip_Theta_Radius){
			if(pos.perp() >= R - pitch/2. && pos.perp() < R + pitch/2.) {
				// match!
				break;
			}
			++theta_strip;
		}
		if(theta_strip == geo.Si_Strip_Theta_Radius.size() + 1) {
			// Strip not found			
			// Check for error in floating point
			if(
				pos.perp() < geo.Si_Strip_Theta_Radius.front() - pitch/2. &&
				fabs(pos.perp() - (geo.Si_Strip_Theta_Radius.front() - pitch/2.)) < 1e-3) {
				// okay!
				theta_strip = 1; // inner strip
			}
			else if (
				pos.perp() >= geo.Si_Strip_Theta_Radius.back() + pitch/2. &&
				fabs(pos.perp() - (geo.Si_Strip_Theta_Radius.back() + pitch/2.)) < 1e-3) {
				// okay!
				--theta_strip; // outer strip
			}
			else {
				//
				// Genuine problem
				std::cerr << "\nWARNING: Theta strip number not found: "
									<< "Radius of hit [mm]: " << pos.perp()/mm << "\n";
				std::cerr << "Strip radii:\n";
				for(const double& R : geo.Si_Strip_Theta_Radius){
					std::cout << "   " << (R-pitch/2.)/mm << " - " << (R+pitch/2.)/mm << "\n";
				}
				std::cerr << "---------\n";
			}
		}
	}
	//
	// Add resolutions
	double eres = AnnularTelescope_NS::SiResoEnergy; // absolute
	double energy_theta = -1000;
	while(energy_theta < 0) { // make sure it's a positive number
		energy_theta = RandGauss::shoot(energy, eres);
	}
	double energy_phi = -1000;
	while(energy_phi < 0) { // make sure it's a positive number
		energy_phi = RandGauss::shoot(energy, eres);
	}
	double tres = AnnularTelescope_NS::SiResoTime; // absolute
	double time_theta = RandGauss::shoot(time, tres);
	double time_phi = RandGauss::shoot(time, tres);
	//
	// Set Si energy and time
	m_Event->SetSiThetaEnergy(detector_number, theta_strip, energy_theta);
	m_Event->SetSiPhiEnergy(detector_number, phi_strip, energy_phi);
	m_Event->SetSiThetaTime(detector_number, theta_strip, time_theta);
	m_Event->SetSiPhiTime(detector_number, phi_strip, time_phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void AnnularTelescope::InitializeScorers() {
	// CSI
	//
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_AnnularTelescopeScorer_CsI = 
		CheckScorer("AnnularTelescopeScorer_CsI",already_exist) ;

  if(!already_exist) {
		// Otherwise the scorer is initialised
		vector<int> level; level.push_back(0);
		G4VPrimitiveScorer* Calorimeter_CsI =
			new CalorimeterScorers::PS_Calorimeter(
				"Calorimeter_CsI",level, 0) ;
		//and register it to the multifunctionnal detector
		m_AnnularTelescopeScorer_CsI->RegisterPrimitive(Calorimeter_CsI);
		G4SDManager::GetSDMpointer()->AddNewDetector(m_AnnularTelescopeScorer_CsI) ;
	}
	// SI
	//
  // This check is necessary in case the geometry is reloaded
	already_exist = false; 
  m_AnnularTelescopeScorer_Si = 
		CheckScorer("AnnularTelescopeScorer_Si",already_exist) ;

  if(!already_exist) {
		// Otherwise the scorer is initialised
		vector<int> level; level.push_back(0);
		G4VPrimitiveScorer* Calorimeter_Si =
			new CalorimeterScorers::PS_Calorimeter(
				"Calorimeter_Si",level, 0) ;
		//and register it to the multifunctionnal detector
		m_AnnularTelescopeScorer_Si->RegisterPrimitive(Calorimeter_Si);
		G4SDManager::GetSDMpointer()->AddNewDetector(m_AnnularTelescopeScorer_Si) ;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* AnnularTelescope::Construct(){
  return  (NPS::VDetector*) new AnnularTelescope();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
class proxy_nps_AnnularTelescope{
public:
	proxy_nps_AnnularTelescope(){
		NPS::DetectorFactory::getInstance()->AddToken("AnnularTelescope","AnnularTelescope");
		NPS::DetectorFactory::getInstance()->AddDetector("AnnularTelescope",AnnularTelescope::Construct);
	}
};

proxy_nps_AnnularTelescope p_nps_AnnularTelescope;
}
