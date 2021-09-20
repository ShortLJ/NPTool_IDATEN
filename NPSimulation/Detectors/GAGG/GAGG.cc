/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Valerian Alcindor  contact address: *
 *                                                                           *
 * Creation Date  : October 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  GAGG simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>
#include <limits>
#include <sstream>
// G4 Geometry object
#include "G4Box.hh"
#include "G4Tubs.hh"

// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// G4 various object
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// NPTool header
#include "CalorimeterScorers.hh"
#include "GAGG.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "RootOutput.h"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace GAGG_NS {
// Energy and time Resolution
const double EnergyThreshold = 0.1 * MeV;
const double ResoTime        = 4.5 * ns;
// const double ResoEnergy = 20*keV ;
const double ResoEnergy = 4.2 / 100.;
// const double ResoEnergy = 0.00001*keV ;
const double Radius    = 20 * mm;
const double Width     = 40 * mm;
const double Thickness = 80 * mm;
const string Material = "GAGG";
// const string Material = "CsI";
} // namespace GAGG_NS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// GAGG Specific Method
GAGG::GAGG() {
  m_Event               = new TGAGGData();
  m_GAGGScorer          = 0;
  m_SquareDetector      = 0;
  m_CylindricalDetector = 0;

  // RGB Color + Transparency
  m_VisSquare   = new G4VisAttributes(G4Colour(1, 1, 0, 0.5));
  m_VisCylinder = new G4VisAttributes(G4Colour(1, 1, 0, 0.5));
}

GAGG::~GAGG() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GAGG::AddDetector(G4ThreeVector POS, string Shape) {
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in
  // G4
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GAGG::AddDetector(double R, double Theta, double Phi, string Shape) {
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* GAGG::BuildSquareDetector() {
  if (!m_SquareDetector) {
    G4Box*      box = new G4Box("GAGG_Box", GAGG_NS::Width * 0.5,
                           GAGG_NS::Width * 0.5, GAGG_NS::Thickness * 0.5);
    G4Material* DetectorMaterial
        = MaterialManager::getInstance()->GetMaterialFromLibrary(
            GAGG_NS::Material);
    m_SquareDetector
        = new G4LogicalVolume(box, DetectorMaterial, "logic_GAGG_Box", 0, 0, 0);
    m_SquareDetector->SetVisAttributes(m_VisSquare);
    m_SquareDetector->SetSensitiveDetector(m_GAGGScorer);
  }
  return m_SquareDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* GAGG::BuildCylindricalDetector() {
  if (!m_CylindricalDetector) {
    G4Tubs* tub = new G4Tubs("GAGG_Cyl", 0, GAGG_NS::Radius,
                             GAGG_NS::Thickness * 0.5, 0, 360 * deg);

    G4Material* DetectorMaterial
        = MaterialManager::getInstance()->GetMaterialFromLibrary(
            GAGG_NS::Material);
    m_CylindricalDetector
        = new G4LogicalVolume(tub, DetectorMaterial, "logic_GAGG_tub", 0, 0, 0);
    m_CylindricalDetector->SetVisAttributes(m_VisSquare);
    m_CylindricalDetector->SetSensitiveDetector(m_GAGGScorer);
  }
  return m_CylindricalDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void GAGG::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("GAGG");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> cart = {"POS", "Shape"};
  vector<string> sphe = {"R", "Theta", "Phi", "Shape"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  GAGG " << i + 1 << endl;

      G4ThreeVector Pos
          = NPS::ConvertVector(blocks[i]->GetTVector3("POS", "mm"));
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos, Shape);
    } else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  GAGG " << i + 1 << endl;
      double R     = blocks[i]->GetDouble("R", "mm");
      double Theta = blocks[i]->GetDouble("Theta", "deg");
      double Phi   = blocks[i]->GetDouble("Phi", "deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R, Theta, Phi, Shape);
    } else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void GAGG::ConstructDetector(G4LogicalVolume* world) {
  for (unsigned short i = 0; i < m_R.size(); i++) {

    G4double      wX      = m_R[i] * sin(m_Theta[i]) * cos(m_Phi[i]);
    G4double      wY      = m_R[i] * sin(m_Theta[i]) * sin(m_Phi[i]);
    G4double      wZ      = m_R[i] * cos(m_Theta[i]);
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ);
    // So the face of the detector is at R instead of the middle
    Det_pos += Det_pos.unit() * GAGG_NS::Thickness * 0.5;
    // Building Detector reference frame
    G4double      ii = cos(m_Theta[i]) * cos(m_Phi[i]);
    G4double      jj = cos(m_Theta[i]) * sin(m_Phi[i]);
    G4double      kk = -sin(m_Theta[i]);
    G4ThreeVector Y(ii, jj, kk);
    G4ThreeVector w = Det_pos.unit();
    G4ThreeVector u = w.cross(Y);
    G4ThreeVector v = w.cross(u);
    v               = v.unit();
    u               = u.unit();

    G4RotationMatrix* Rot = new G4RotationMatrix(u, v, w);

    if (m_Shape[i] == "Cylindrical") {
      new G4PVPlacement(G4Transform3D(*Rot, Det_pos),
                        BuildCylindricalDetector(), "GAGG", world, false,
                        i + 1);
    }

    else if (m_Shape[i] == "Square") {
      new G4PVPlacement(G4Transform3D(*Rot, Det_pos), BuildSquareDetector(),
                        "GAGG", world, false, i + 1);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void GAGG::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree*      pTree     = pAnalysis->GetTree();
  if (!pTree->FindBranch("GAGG")) {
    pTree->Branch("GAGG", "TGAGGData", &m_Event);
  }
  pTree->SetBranchAddress("GAGG", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void GAGG::ReadSensitive(const G4Event*) {
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer
      = (CalorimeterScorers::PS_Calorimeter*)m_GAGGScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult();
  for (unsigned int i = 0; i < size; i++) {
    vector<unsigned int> level = Scorer->GetLevel(i);
    double               Energy
        = RandGauss::shoot(Scorer->GetEnergy(i),
                           GAGG_NS::ResoEnergy / 2.35 * Scorer->GetEnergy(i));
    if (Energy > GAGG_NS::EnergyThreshold) {
      double Time = RandGauss::shoot(Scorer->GetTime(i), GAGG_NS::ResoTime);
      int    DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr, Energy);
      m_Event->SetTime(DetectorNbr, Time);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void GAGG::InitializeScorers() {
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;
  m_GAGGScorer       = CheckScorer("GAGGScorer", already_exist);

  if (already_exist)
    return;

  // Otherwise the scorer is initialised
  vector<int> level;
  level.push_back(0);
  G4VPrimitiveScorer* Calorimeter
      = new CalorimeterScorers::PS_Calorimeter("Calorimeter", level, 0);
  G4VPrimitiveScorer* Interaction = new InteractionScorers::PS_Interactions(
      "Interaction", ms_InterCoord, 0);
  // and register it to the multifunctionnal detector
  m_GAGGScorer->RegisterPrimitive(Calorimeter);
  m_GAGGScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_GAGGScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* GAGG::Construct() { return (NPS::VDetector*)new GAGG(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_GAGG {
public:
  proxy_nps_GAGG() {
    NPS::DetectorFactory::getInstance()->AddToken("GAGG", "GAGG");
    NPS::DetectorFactory::getInstance()->AddDetector("GAGG", GAGG::Construct);
  }
};

proxy_nps_GAGG p_nps_GAGG;
}
