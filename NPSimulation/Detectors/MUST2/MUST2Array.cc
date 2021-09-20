/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    : October 2009                                             *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This file describe the MUST2 charge particle Detector                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * MUST2 is a modular array made of Telescope (1 to 8 telescope). Each       *
 *  Telescope is made of Three Stage:                                        *
 *  - A 300um Silicium, double-sided strip                                   *
 *  - 16 Si(Li) pad                                                          *
 *  - 16 CsI scintillator Crystal                                            *
 *****************************************************************************/
#include "MUST2Array.hh"

// Geant4
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4MaterialTable.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4VisAttributes.hh"
#include "Randomize.hh"
// NPS
#include "CalorimeterScorers.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
// NPL
#include "NPCore.h"
// ROOT
#include "RootOutput.h"

// CLHEP
#include "CLHEP/Random/RandGauss.h"

// STL
#include <cmath>
#include <set>
#include <sstream>
#include <string>
using namespace std;
using namespace CLHEP;
using namespace MUST2;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// MUST2Array Specific Method
MUST2Array::MUST2Array() {
  m_Event = new TMust2Data();
  InitializeMaterial();
  m_StripScorer = 0;
  m_SiLiScorer  = 0;
  m_CsIScorer   = 0;
}

MUST2Array::~MUST2Array() {}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MUST2Array::AddTelescope(G4ThreeVector X1_Y1, G4ThreeVector X128_Y1,
    G4ThreeVector X1_Y128, G4ThreeVector X128_Y128,
    bool wSi, bool wSiLi, bool wCsI) {
  m_DefinitionType.push_back(true);

  m_X1_Y1.push_back(X1_Y1);
  m_X128_Y1.push_back(X128_Y1);
  m_X1_Y128.push_back(X1_Y128);
  m_X128_Y128.push_back(X128_Y128);
  m_wSi.push_back(wSi);
  m_wSiLi.push_back(wSiLi);
  m_wCsI.push_back(wCsI);

  m_R.push_back(0);
  m_Theta.push_back(0);
  m_Phi.push_back(0);
  m_beta_u.push_back(0);
  m_beta_v.push_back(0);
  m_beta_w.push_back(0);
}

void MUST2Array::AddTelescope(G4double R, G4double Theta, G4double Phi,
    G4double beta_u, G4double beta_v, G4double beta_w,
    bool wSi, bool wSiLi, bool wCsI) {
  G4ThreeVector empty = G4ThreeVector(0, 0, 0);

  m_DefinitionType.push_back(false);

  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_beta_u.push_back(beta_u);
  m_beta_v.push_back(beta_v);
  m_beta_w.push_back(beta_w);
  m_wSi.push_back(wSi);
  m_wSiLi.push_back(wSiLi);
  m_wCsI.push_back(wCsI);

  m_X1_Y1.push_back(empty);
  m_X128_Y1.push_back(empty);
  m_X1_Y128.push_back(empty);
  m_X128_Y128.push_back(empty);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MUST2Array::VolumeMaker(G4int TelescopeNumber, G4ThreeVector MMpos,
    G4RotationMatrix* MMrot, bool wSi, bool wSiLi,
    bool wCsI, G4LogicalVolume* world) {
  G4double           NbrTelescopes = TelescopeNumber;
  G4String           DetectorNumber;
  std::ostringstream Number;
  Number << NbrTelescopes;
  DetectorNumber = Number.str();

  ////////////////////////////////////////////////////////////////
  ////////////// Starting Volume Definition //////////////////////
  ////////////////////////////////////////////////////////////////
  G4Trd* solidMM = new G4Trd("MUST2Telescope" + DetectorNumber, 0.5 * FaceFront,
      0.5 * FaceBack, 0.5 * FaceFront, 0.5 * FaceBack,
      0.5 * Length);
  G4LogicalVolume* logicMM = new G4LogicalVolume(
      solidMM, m_MaterialIron, "MUST2Telescope" + DetectorNumber, 0, 0, 0);
  G4String Name = "MUST2Telescope" + DetectorNumber;

  new G4PVPlacement(G4Transform3D(*MMrot, MMpos), logicMM, Name, world, false,
      TelescopeNumber);

  if (m_non_sensitive_part_visiualisation) {
    G4VisAttributes* FrameVisAtt
      = new G4VisAttributes(G4Colour(0.80, 0.80, 0.80));
    FrameVisAtt->SetForceWireframe(true);
    logicMM->SetVisAttributes(FrameVisAtt);
  } else
    logicMM->SetVisAttributes(G4VisAttributes::Invisible);

  G4ThreeVector positionVacBox = G4ThreeVector(0, 0, VacBox_PosZ);

  G4Trd* solidVacBox
    = new G4Trd("solidVacBox", 0.5 * SiliconFace, 0.5 * CsIFaceFront,
        0.5 * SiliconFace, 0.5 * CsIFaceFront, 0.5 * VacBoxThickness);
  G4LogicalVolume* logicVacBox = new G4LogicalVolume(
      solidVacBox, m_MaterialVacuum, "logicVacBox", 0, 0, 0);

  new G4PVPlacement(0, positionVacBox, logicVacBox, Name + "_VacBox", logicMM,
      false, TelescopeNumber);

  logicVacBox->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* SiliconVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));

  ////////////////////////////////////////////////////////////////
  /////////////////Si Strip Construction//////////////////////////
  ////////////////////////////////////////////////////////////////

  if (wSi) {
    G4ThreeVector positionAluStripFront
      = G4ThreeVector(0, 0, AluStripFront_PosZ);
    G4ThreeVector positionAluStripBack = G4ThreeVector(0, 0, AluStripBack_PosZ);

    G4Box* solidAluStrip
      = new G4Box("AluBox", 0.5 * SiliconFace, 0.5 * SiliconFace,
          0.5 * AluStripThickness);
    G4LogicalVolume* logicAluStrip = new G4LogicalVolume(
        solidAluStrip, m_MaterialAluminium, "logicAluStrip", 0, 0, 0);

    new G4PVPlacement(0, positionAluStripFront, logicAluStrip,
        Name + "_AluStripFront", logicMM, false, TelescopeNumber);
    new G4PVPlacement(0, positionAluStripBack, logicAluStrip,
        Name + "_AluStripBack", logicMM, false, TelescopeNumber);

    logicAluStrip->SetVisAttributes(G4VisAttributes::Invisible);

    G4ThreeVector positionSilicon = G4ThreeVector(0, 0, Silicon_PosZ);

    G4Box*           solidSilicon = new G4Box("solidSilicon", 0.5 * SiliconFace,
        0.5 * SiliconFace, 0.5 * SiliconThickness);
    G4LogicalVolume* logicSilicon = new G4LogicalVolume(
        solidSilicon, m_MaterialSilicon, "logicSilicon", 0, 0, 0);

    new G4PVPlacement(0, positionSilicon, logicSilicon, Name + "_Silicon",
        logicMM, false, TelescopeNumber);

    /// Set Silicon strip sensible
    logicSilicon->SetSensitiveDetector(m_StripScorer);

    /// Visualisation of Silicon Strip
    logicSilicon->SetVisAttributes(SiliconVisAtt);
  }

  ////////////////////////////////////////////////////////////////
  //////////////////// SiLi  Construction ////////////////////////
  ////////////////////////////////////////////////////////////////

  if (wSiLi) {
    G4double          SiLiSpace = 8 * mm;
    G4RotationMatrix* rotSiLi   = new G4RotationMatrix(0, 0, 0);
    G4Box* solidSiLi = new G4Box("SiLi", 0.5 * SiliconFace + 0.5 * SiLiSpace,
        0.5 * SiliconFace, 0.5 * SiLiThickness);
    G4LogicalVolume* logicSiLi = new G4LogicalVolume(
        solidSiLi, m_MaterialAluminium, Name + "_SiLi", 0, 0, 0);

    logicSiLi->SetVisAttributes(G4VisAttributes::Invisible);

    new G4PVPlacement(G4Transform3D(*rotSiLi, G4ThreeVector(0, 0, 0)),
        logicSiLi, Name + "_SiLi", logicVacBox, false, 0);

    // SiLi are placed inside of the VacBox...
    // Left/Right define when looking to detector from Si to CsI
    G4double SiLi_HighY_Upper  = 19.86 * mm;
    G4double SiLi_HighY_Center = 25.39 * mm;
    G4double SiLi_WidthX_Left  = 22.85 * mm;
    G4double SiLi_WidthX_Right = 24.9 * mm;
    G4double SiLi_ShiftX       = 0.775 * mm;

    //	SiLi are organized by two group of 8 Up(9 to 15) and Down(1 to 8).
    G4ThreeVector ShiftSiLiUp
      = G4ThreeVector(-0.25 * SiliconFace - 0.5 * SiLiSpace, 0, 0);
    G4ThreeVector ShiftSiLiDown
      = G4ThreeVector(0.25 * SiliconFace + 0.5 * SiLiSpace, 0, 0);

    // SiLi : left side of SiLi detector
    G4Box* solidSiLi_LT
      = new G4Box("SiLi_LT", 0.5 * SiLi_WidthX_Left, 0.5 * SiLi_HighY_Upper,
          0.5 * SiLiThickness);
    G4Box* solidSiLi_RT
      = new G4Box("SiLi_RT", 0.5 * SiLi_WidthX_Right, 0.5 * SiLi_HighY_Upper,
          0.5 * SiLiThickness);
    G4Box* solidSiLi_LC1
      = new G4Box("SiLi_LC1", 0.5 * SiLi_WidthX_Left, 0.5 * SiLi_HighY_Center,
          0.5 * SiLiThickness);
    G4Box* solidSiLi_RC1
      = new G4Box("SiLi_RC1", 0.5 * SiLi_WidthX_Right,
          0.5 * SiLi_HighY_Center, 0.5 * SiLiThickness);
    G4Box* solidSiLi_LB
      = new G4Box("SiLi_LB", 0.5 * SiLi_WidthX_Left, 0.5 * SiLi_HighY_Upper,
          0.5 * SiLiThickness);
    G4Box* solidSiLi_RB
      = new G4Box("SiLi_RB", 0.5 * SiLi_WidthX_Right, 0.5 * SiLi_HighY_Upper,
          0.5 * SiLiThickness);
    G4Box* solidSiLi_LC2
      = new G4Box("SiLi_LC2", 0.5 * SiLi_WidthX_Left, 0.5 * SiLi_HighY_Center,
          0.5 * SiLiThickness);
    G4Box* solidSiLi_RC2
      = new G4Box("SiLi_RC2", 0.5 * SiLi_WidthX_Right,
          0.5 * SiLi_HighY_Center, 0.5 * SiLiThickness);

    G4LogicalVolume* logicSiLi_LT = new G4LogicalVolume(
        solidSiLi_LT, m_MaterialSilicon, "SiLi_LT", 0, 0, 0);
    G4LogicalVolume* logicSiLi_RT = new G4LogicalVolume(
        solidSiLi_RT, m_MaterialSilicon, "SiLi_RT", 0, 0, 0);
    G4LogicalVolume* logicSiLi_LC1 = new G4LogicalVolume(
        solidSiLi_LC1, m_MaterialSilicon, "SiLi_LC1", 0, 0, 0);
    G4LogicalVolume* logicSiLi_RC1 = new G4LogicalVolume(
        solidSiLi_RC1, m_MaterialSilicon, "SiLi_RC1", 0, 0, 0);
    G4LogicalVolume* logicSiLi_LB = new G4LogicalVolume(
        solidSiLi_LB, m_MaterialSilicon, "SiLi_LB", 0, 0, 0);
    G4LogicalVolume* logicSiLi_RB = new G4LogicalVolume(
        solidSiLi_RB, m_MaterialSilicon, "SiLi_RB", 0, 0, 0);
    G4LogicalVolume* logicSiLi_LC2 = new G4LogicalVolume(
        solidSiLi_LC2, m_MaterialSilicon, "SiLi_LC2", 0, 0, 0);
    G4LogicalVolume* logicSiLi_RC2 = new G4LogicalVolume(
        solidSiLi_RC2, m_MaterialSilicon, "SiLi_RC2", 0, 0, 0);

    G4double interSiLi = 0.5 * mm;

    // Top
    G4ThreeVector positionSiLi_LT_up = G4ThreeVector(
        -0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
        0.5 * SiLi_HighY_Upper + SiLi_HighY_Center + 1.5 * interSiLi, 0);

    positionSiLi_LT_up += ShiftSiLiUp;

    G4ThreeVector positionSiLi_RT_up = G4ThreeVector(
        0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
        0.5 * SiLi_HighY_Upper + SiLi_HighY_Center + 1.5 * interSiLi, 0);

    positionSiLi_RT_up += ShiftSiLiUp;

    G4ThreeVector positionSiLi_LC1_up
      = G4ThreeVector(-0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
          0.5 * SiLi_HighY_Center + 0.5 * interSiLi, 0);

    positionSiLi_LC1_up += ShiftSiLiUp;

    G4ThreeVector positionSiLi_RC1_up
      = G4ThreeVector(0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
          0.5 * SiLi_HighY_Center + 0.5 * interSiLi, 0);

    positionSiLi_RC1_up += ShiftSiLiUp;

    G4ThreeVector positionSiLi_LB_up = G4ThreeVector(
        -0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
        -0.5 * SiLi_HighY_Upper - SiLi_HighY_Center - 1.5 * interSiLi, 0);

    positionSiLi_LB_up += ShiftSiLiUp;

    G4ThreeVector positionSiLi_RB_up = G4ThreeVector(
        0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
        -0.5 * SiLi_HighY_Upper - SiLi_HighY_Center - 1.5 * interSiLi, 0);

    positionSiLi_RB_up += ShiftSiLiUp;

    G4ThreeVector positionSiLi_LC2_up
      = G4ThreeVector(-0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
          -0.5 * SiLi_HighY_Center - 0.5 * interSiLi, 0);

    positionSiLi_LC2_up += ShiftSiLiUp;

    G4ThreeVector positionSiLi_RC2_up
      = G4ThreeVector(0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
          -0.5 * SiLi_HighY_Center - 0.5 * interSiLi, 0);

    positionSiLi_RC2_up += ShiftSiLiUp;

    // Down
    G4ThreeVector positionSiLi_LT_down = G4ThreeVector(
        -0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
        0.5 * SiLi_HighY_Upper + SiLi_HighY_Center + 1.5 * interSiLi, 0);

    positionSiLi_LT_down += ShiftSiLiDown;

    G4ThreeVector positionSiLi_RT_down = G4ThreeVector(
        0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
        0.5 * SiLi_HighY_Upper + SiLi_HighY_Center + 1.5 * interSiLi, 0);

    positionSiLi_RT_down += ShiftSiLiDown;

    G4ThreeVector positionSiLi_LC1_down
      = G4ThreeVector(-0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
          0.5 * SiLi_HighY_Center + 0.5 * interSiLi, 0);

    positionSiLi_LC1_down += ShiftSiLiDown;

    G4ThreeVector positionSiLi_RC1_down
      = G4ThreeVector(0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
          0.5 * SiLi_HighY_Center + 0.5 * interSiLi, 0);

    positionSiLi_RC1_down += ShiftSiLiDown;

    G4ThreeVector positionSiLi_LB_down = G4ThreeVector(
        -0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
        -0.5 * SiLi_HighY_Upper - SiLi_HighY_Center - 1.5 * interSiLi, 0);

    positionSiLi_LB_down += ShiftSiLiDown;

    G4ThreeVector positionSiLi_RB_down = G4ThreeVector(
        0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
        -0.5 * SiLi_HighY_Upper - SiLi_HighY_Center - 1.5 * interSiLi, 0);

    positionSiLi_RB_down += ShiftSiLiDown;

    G4ThreeVector positionSiLi_LC2_down
      = G4ThreeVector(-0.5 * SiLi_WidthX_Left - interSiLi - SiLi_ShiftX,
          -0.5 * SiLi_HighY_Center - 0.5 * interSiLi, 0);

    positionSiLi_LC2_down += ShiftSiLiDown;

    G4ThreeVector positionSiLi_RC2_down
      = G4ThreeVector(0.5 * SiLi_WidthX_Right - SiLi_ShiftX,
          -0.5 * SiLi_HighY_Center - 0.5 * interSiLi, 0);

    positionSiLi_RC2_down += ShiftSiLiDown;

    new G4PVPlacement(0, positionSiLi_RT_down, logicSiLi_RT,
        Name + "_SiLi_Pad1", logicSiLi, false, 1);
    new G4PVPlacement(0, positionSiLi_LT_down, logicSiLi_LT,
        Name + "_SiLi_Pad2", logicSiLi, false, 2);
    new G4PVPlacement(0, positionSiLi_RC1_down, logicSiLi_RC1,
        Name + "_SiLi_Pad3", logicSiLi, false, 3);
    new G4PVPlacement(0, positionSiLi_LC1_down, logicSiLi_LC1,
        Name + "_SiLi_Pad4", logicSiLi, false, 4);
    new G4PVPlacement(0, positionSiLi_LC2_down, logicSiLi_LC2,
        Name + "_SiLi_Pad5", logicSiLi, false, 5);
    new G4PVPlacement(0, positionSiLi_RC2_down, logicSiLi_RC2,
        Name + "_SiLi_Pad6", logicSiLi, false, 6);
    new G4PVPlacement(0, positionSiLi_LB_down, logicSiLi_LB,
        Name + "_SiLi_Pad7", logicSiLi, false, 7);
    new G4PVPlacement(0, positionSiLi_RB_down, logicSiLi_RB,
        Name + "_SiLi_Pad8", logicSiLi, false, 8);
    new G4PVPlacement(0, positionSiLi_LT_up, logicSiLi_LT, Name + "_SiLi_Pad9",
        logicSiLi, false, 9);
    new G4PVPlacement(0, positionSiLi_RT_up, logicSiLi_RT, Name + "_SiLi_Pad10",
        logicSiLi, false, 10);
    new G4PVPlacement(0, positionSiLi_LC1_up, logicSiLi_LC1,
        Name + "_SiLi_Pad11", logicSiLi, false, 11);
    new G4PVPlacement(0, positionSiLi_RC1_up, logicSiLi_RC1,
        Name + "_SiLi_Pad12", logicSiLi, false, 12);
    new G4PVPlacement(0, positionSiLi_RC2_up, logicSiLi_RC2,
        Name + "_SiLi_Pad13", logicSiLi, false, 13);
    new G4PVPlacement(0, positionSiLi_LC2_up, logicSiLi_LC2,
        Name + "_SiLi_Pad14", logicSiLi, false, 14);
    new G4PVPlacement(0, positionSiLi_RB_up, logicSiLi_RB, Name + "_SiLi_Pad15",
        logicSiLi, false, 15);
    new G4PVPlacement(0, positionSiLi_LB_up, logicSiLi_LB, Name + "_SiLi_Pad16",
        logicSiLi, false, 16);

    // Set SiLi sensible
    logicSiLi_LT->SetSensitiveDetector(m_SiLiScorer);
    logicSiLi_RT->SetSensitiveDetector(m_SiLiScorer);
    logicSiLi_LC1->SetSensitiveDetector(m_SiLiScorer);
    logicSiLi_RC1->SetSensitiveDetector(m_SiLiScorer);

    logicSiLi_LB->SetSensitiveDetector(m_SiLiScorer);
    logicSiLi_RB->SetSensitiveDetector(m_SiLiScorer);
    logicSiLi_LC2->SetSensitiveDetector(m_SiLiScorer);
    logicSiLi_RC2->SetSensitiveDetector(m_SiLiScorer);

    // Mark blue a SiLi to see telescope orientation
    G4VisAttributes* SiLiVisAtt = new G4VisAttributes(G4Colour(0.3, 1, 0.3));

    logicSiLi_LT->SetVisAttributes(SiLiVisAtt);
    logicSiLi_RT->SetVisAttributes(SiLiVisAtt);
    logicSiLi_LC1->SetVisAttributes(SiLiVisAtt);
    logicSiLi_RC1->SetVisAttributes(SiLiVisAtt);

    logicSiLi_LB->SetVisAttributes(SiLiVisAtt);
    logicSiLi_RB->SetVisAttributes(SiLiVisAtt);
    logicSiLi_LC2->SetVisAttributes(SiLiVisAtt);
    logicSiLi_RC2->SetVisAttributes(SiLiVisAtt);

    delete rotSiLi;
  }

  ////////////////////////////////////////////////////////////////
  //////////////////// CsI  Construction//////////////////////////
  ////////////////////////////////////////////////////////////////

  if (wCsI) {
    m_MaterialMyl
      = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");
    m_MaterialCsI
      = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");

    G4ThreeVector positionCsI = G4ThreeVector(0, 0, CsI_PosZ);
    G4Trd*        solidCsI
      = new G4Trd("csI", 0.5 * CsIFaceFront, 0.5 * CsIFaceBack,
          0.5 * CsIFaceFront, 0.5 * CsIFaceBack, 0.5 * CsIThickness);

    G4LogicalVolume* logicCsI = new G4LogicalVolume(
        solidCsI, m_MaterialAluminium, Name + "_CsI_Mylar", 0, 0, 0);
    new G4PVPlacement(0, positionCsI, logicCsI, Name + "_CsI_Mylar", logicMM,
        false, 0);

    G4ThreeVector positionMylarCsI
      = G4ThreeVector(0, 0, MylarCsIThickness * 0.5 - CsIThickness * 0.5);

    G4Box* solidMylarCsI
      = new G4Box("MylarCsIBox", 0.5 * CsIFaceFront, 0.5 * CsIFaceFront,
          0.5 * MylarCsIThickness);
    G4LogicalVolume* logicMylarCsI = new G4LogicalVolume(
        solidMylarCsI, m_MaterialMyl, Name + "_CsI_Mylar", 0, 0, 0);

    new G4PVPlacement(0, positionMylarCsI, logicMylarCsI, Name + "_CsI_Mylar",
        logicCsI, false, 0);

    logicCsI->SetVisAttributes(G4VisAttributes::Invisible);
    logicMylarCsI->SetVisAttributes(G4VisAttributes::Invisible);

    // Cristal1
    G4Trap* solidCristal1 = new G4Trap(
        "Cristal1", 40. * mm / 2., 6.693896 * deg, 41.97814 * deg,
        33.1 * mm / 2., 37.39 * mm / 2., 37.39 * mm / 2., 0. * deg,
        26.9 * mm / 2., 30.41 * mm / 2., 30.41 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal1 = new G4LogicalVolume(
        solidCristal1, m_MaterialCsI, Name + "_CsI_Cristal1", 0, 0, 0);

    // Cristal2
    G4Trap* solidCristal2 = new G4Trap(
        "Cristal2", 40. * mm / 2., 17.8836 * deg, (74.3122 + 180) * deg,
        43.49 * mm / 2., 37.39 * mm / 2., 37.39 * mm / 2., 0. * deg,
        31.0377 * mm / 2., 30.41 * mm / 2., 30.41 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal2 = new G4LogicalVolume(
        solidCristal2, m_MaterialCsI, Name + "_CsI_Cristal2", 0, 0, 0);

    // Cristal3
    G4Trap* solidCristal3 = new G4Trap(
        "Cristal3", 40. * mm / 2., 18.243 * deg, 13.5988 * deg, 33.11 * mm / 2.,
        39.25 * mm / 2., 39.25 * mm / 2., 0. * deg, 26.91 * mm / 2.,
        27.58 * mm / 2., 27.58 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal3 = new G4LogicalVolume(
        solidCristal3, m_MaterialCsI, Name + "_CsI_Cristal3", 0, 0, 0);

    // Cristal4

    G4Trap* solidCristal4 = new G4Trap(
        "Cristal4", 40. * mm / 2., 24.0482 * deg, 44.1148 * deg,
        43.49 * mm / 2., 39.19 * mm / 2., 39.19 * mm / 2., 0. * deg,
        31.04 * mm / 2., 27.52 * mm / 2., 27.52 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal4 = new G4LogicalVolume(
        solidCristal4, m_MaterialCsI, Name + "_CsI_Cristal4", 0, 0, 0);

    // Cristal1s

    G4Trap* solidCristal1s = new G4Trap(
        "Cristal1s", 40. * mm / 2., 6.693896 * deg, -41.97814 * deg,
        33.1 * mm / 2., 37.39 * mm / 2., 37.39 * mm / 2., 0. * deg,
        26.9 * mm / 2., 30.41 * mm / 2., 30.41 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal1s = new G4LogicalVolume(
        solidCristal1s, m_MaterialCsI, Name + "_CsI_Cristal1s", 0, 0, 0);

    // Cristal2s

    G4Trap* solidCristal2s = new G4Trap(
        "Cristal2s", 40. * mm / 2., 17.8836 * deg, -(74.3122 + 180) * deg,
        43.49 * mm / 2., 37.39 * mm / 2., 37.39 * mm / 2., 0. * deg,
        31.0377 * mm / 2., 30.41 * mm / 2., 30.41 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal2s = new G4LogicalVolume(
        solidCristal2s, m_MaterialCsI, Name + "_CsI_Cristal2s", 0, 0, 0);

    // Cristal3s

    G4Trap* solidCristal3s = new G4Trap(
        "Cristal3s", 40. * mm / 2., 18.243 * deg, -13.5988 * deg,
        33.11 * mm / 2., 39.25 * mm / 2., 39.25 * mm / 2., 0. * deg,
        26.91 * mm / 2., 27.58 * mm / 2., 27.58 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal3s = new G4LogicalVolume(
        solidCristal3s, m_MaterialCsI, Name + "_CsI_Cristal3s", 0, 0, 0);

    // Cristal4s

    G4Trap* solidCristal4s = new G4Trap(
        "Cristal4s", 40. * mm / 2., 24.0482 * deg, -44.1148 * deg,
        43.49 * mm / 2., 39.19 * mm / 2., 39.19 * mm / 2., 0. * deg,
        31.04 * mm / 2., 27.52 * mm / 2., 27.52 * mm / 2., 0. * deg);
    G4LogicalVolume* logicCristal4s = new G4LogicalVolume(
        solidCristal4s, m_MaterialCsI, Name + "_CsI_Cristal4s", 0, 0, 0);

    G4double XEdge1 = 16.96 * mm + DistInterCsI * 0.5;
    G4double YEdge1 = 15.01 * mm + DistInterCsI * 0.5;
    G4double XEdge2 = 50.63 * mm + DistInterCsI * 1.5;
    G4double YEdge2 = 48.67 * mm + DistInterCsI * 1.5;

    G4ThreeVector positionCristal1 = G4ThreeVector(-XEdge1, YEdge1, 0 * mm);
    G4ThreeVector positionCristal2 = G4ThreeVector(-XEdge1, YEdge2, 0);
    G4ThreeVector positionCristal3 = G4ThreeVector(-XEdge2, YEdge1, 0);
    G4ThreeVector positionCristal4 = G4ThreeVector(-XEdge2, YEdge2, 0);

    new G4PVPlacement(Rotation(180., 0., 0.), positionCristal1, logicCristal1,
        Name + "_CsI_Cristal1", logicCsI, false, 6); // 1
    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal2, logicCristal2,
        Name + "_CsI_Cristal2", logicCsI, false, 3); // 2
    new G4PVPlacement(Rotation(180., 0., 0.), positionCristal3, logicCristal3,
        Name + "_CsI_Cristal3", logicCsI, false, 5); // 3
    new G4PVPlacement(Rotation(180., 0., 0.), positionCristal4, logicCristal4,
        Name + "_CsI_Cristal4", logicCsI, false, 4); // 4

    G4ThreeVector positionCristal1b = G4ThreeVector(XEdge1, -YEdge1, 0 * mm);
    G4ThreeVector positionCristal2b = G4ThreeVector(XEdge1, -YEdge2, 0);
    G4ThreeVector positionCristal3b = G4ThreeVector(XEdge2, -YEdge1, 0);
    G4ThreeVector positionCristal4b = G4ThreeVector(XEdge2, -YEdge2, 0);

    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal1b,
        logicCristal1, Name + "_CsI_Cristal5", logicCsI, false,
        11); // 5
    new G4PVPlacement(Rotation(180., 0., 0.), positionCristal2b, logicCristal2,
        Name + "_CsI_Cristal6", logicCsI, false, 15); // 6
    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal3b,
        logicCristal3, Name + "_CsI_Cristal7", logicCsI, false,
        12); // 7
    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal4b,
        logicCristal4, Name + "_CsI_Cristal8", logicCsI, false,
        16); // 8

    G4ThreeVector positionCristal1s = G4ThreeVector(-XEdge1, -YEdge1, 0 * mm);
    G4ThreeVector positionCristal2s = G4ThreeVector(-XEdge1, -YEdge2, 0);
    G4ThreeVector positionCristal3s = G4ThreeVector(-XEdge2, -YEdge1, 0);
    G4ThreeVector positionCristal4s = G4ThreeVector(-XEdge2, -YEdge2, 0);

    new G4PVPlacement(Rotation(180., 0., 0.), positionCristal1s, logicCristal1s,
        Name + "_CsI_Cristal9", logicCsI, false, 10); // 9
    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal2s,
        logicCristal2s, Name + "_CsI_Cristal10", logicCsI, false,
        14); // 10
    new G4PVPlacement(Rotation(180., 0., 0.), positionCristal3s, logicCristal3s,
        Name + "_CsI_Cristal11", logicCsI, false, 9); // 11
    new G4PVPlacement(Rotation(180., 0., 0.), positionCristal4s, logicCristal4s,
        Name + "_CsI_Cristal12", logicCsI, false, 13); // 12

    G4ThreeVector positionCristal1sb = G4ThreeVector(XEdge1, YEdge1, 0 * mm);
    G4ThreeVector positionCristal2sb = G4ThreeVector(XEdge1, YEdge2, 0);
    G4ThreeVector positionCristal3sb = G4ThreeVector(XEdge2, YEdge1, 0);
    G4ThreeVector positionCristal4sb = G4ThreeVector(XEdge2, YEdge2, 0);

    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal1sb,
        logicCristal1s, Name + "_CsI_Cristal13", logicCsI, false,
        7); // 13
    new G4PVPlacement(Rotation(180, 0, 0), positionCristal2sb, logicCristal2s,
        Name + "_CsI_Cristal14", logicCsI, false, 2); // 14
    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal3sb,
        logicCristal3s, Name + "_CsI_Cristal15", logicCsI, false,
        8); // 15
    new G4PVPlacement(Rotation(180., 0., 180.), positionCristal4sb,
        logicCristal4s, Name + "_CsI_Cristal16", logicCsI, false,
        1); // 16

    /// Set CsI sensible
    logicCristal1->SetSensitiveDetector(m_CsIScorer);
    logicCristal2->SetSensitiveDetector(m_CsIScorer);
    logicCristal3->SetSensitiveDetector(m_CsIScorer);
    logicCristal4->SetSensitiveDetector(m_CsIScorer);

    logicCristal1s->SetSensitiveDetector(m_CsIScorer);
    logicCristal2s->SetSensitiveDetector(m_CsIScorer);
    logicCristal3s->SetSensitiveDetector(m_CsIScorer);
    logicCristal4s->SetSensitiveDetector(m_CsIScorer);

    // Mark blue a CsI corner crystal to see telescope orientation
    G4VisAttributes* CsIVisAtt = new G4VisAttributes(G4Colour(1, 0.5, 0));
    logicCristal1->SetVisAttributes(CsIVisAtt);
    logicCristal2->SetVisAttributes(CsIVisAtt);
    logicCristal3->SetVisAttributes(CsIVisAtt);
    logicCristal4->SetVisAttributes(CsIVisAtt);
    logicCristal1s->SetVisAttributes(CsIVisAtt);
    logicCristal2s->SetVisAttributes(CsIVisAtt);
    logicCristal3s->SetVisAttributes(CsIVisAtt);
    logicCristal4s->SetVisAttributes(CsIVisAtt);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void MUST2Array::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("M2Telescope");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " telescope found" << endl;
  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (NPOptionManager::getInstance()->GetVerboseLevel())
      cout << endl << "//// Must 2 Telecope " << i + 1 << endl;
    // Cartesian Case
    vector<string> cart
      = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128", "SI", "SILI", "CSI"};
    // Spherical Case
    vector<string> sphe = {"R", "THETA", "PHI", "BETA", "SI", "SILI", "CSI"};

    if (blocks[i]->HasTokenList(cart)) {
      G4ThreeVector A
        = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1", "mm"));
      G4ThreeVector B
        = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y1", "mm"));
      G4ThreeVector C
        = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y128", "mm"));
      G4ThreeVector D
        = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y128", "mm"));
      int SI   = blocks[i]->GetInt("SI");
      int SILI = blocks[i]->GetInt("SILI");
      int CSI  = blocks[i]->GetInt("CSI");
      AddTelescope(A, B, C, D, SI == 1, SILI == 1, CSI == 1);
    }

    else if (blocks[i]->HasTokenList(sphe)) {

      double         Theta = blocks[i]->GetDouble("THETA", "deg");
      double         Phi   = blocks[i]->GetDouble("PHI", "deg");
      double         R     = blocks[i]->GetDouble("R", "mm");
      vector<double> beta  = blocks[i]->GetVectorDouble("BETA", "deg");
      int            SI    = blocks[i]->GetInt("SI");
      int            SILI  = blocks[i]->GetInt("SILI");
      int            CSI   = blocks[i]->GetInt("CSI");
      AddTelescope(R, Theta, Phi, beta[0], beta[1], beta[2], SI == 1, SILI == 1,
          CSI == 1);
    }

    else {
      cout << "WARNING: Missing token for M2Telescope blocks, check your input "
        "file"
        << endl;
      exit(1);
    }

    if (blocks[i]->GetString("VIS") == "all")
      m_non_sensitive_part_visiualisation = true;
  }

  ////////////////////
  //Read the thresholds from the analysis config 
  ////////////////////  
  
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigMust2.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigMust2.dat found: Default parameters loaded for "
            "Analysis "
         << FileName << endl;
    return;
  }
  cout << " Loading user parameters for Analysis from ConfigMust2.dat " << endl;

  // read analysis config file
  string LineBuffer, DataBuffer, whatToDo;

  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);
    // search for "header"
    if (LineBuffer.compare(0, 11, "ConfigMust2") == 0)
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {

      whatToDo = "";
      AnalysisConfigFile >> whatToDo;
      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }
      //Resolutions
      else if (whatToDo == "SI_E_RESOLUTION") {
        AnalysisConfigFile >> DataBuffer;
        ResoStrip = atof(DataBuffer.c_str());
        ResoStrip = ResoStrip*keV/2.35;
        cout << whatToDo << " " << ResoStrip  << " MeV/2.35 "<< endl;
      }
      else if (whatToDo == "SILI_E_RESOLUTION") {
        AnalysisConfigFile >> DataBuffer;
        ResoSiLi = atof(DataBuffer.c_str());
        ResoSiLi = ResoSiLi*keV/2.35;
        cout << whatToDo << " " << ResoSiLi  << " MeV/2.35 "<< endl;
      }
      else if (whatToDo == "CSI_E_RESOLUTION") {
        AnalysisConfigFile >> DataBuffer;
        ResoCsI = atof(DataBuffer.c_str());
        ResoCsI = ResoCsI*keV/2.35;
        cout << whatToDo << " " << ResoCsI  << " MeV/2.35 "<< endl;
      }
      //Time
      else if (whatToDo == "MUST_T_RESOLUTION") {
        AnalysisConfigFile >> DataBuffer;
        ResoTimeMust = atof(DataBuffer.c_str());
        ResoTimeMust = ResoTimeMust*ns/2.35;
        cout << whatToDo << " " << ResoTimeMust  << " ns/2.35 "<< endl;
      }
      else if (whatToDo == "SI_T_OFFSET") {
        AnalysisConfigFile >> DataBuffer;
        TimeOffset = atof(DataBuffer.c_str());
        TimeOffset = TimeOffset*ns;
        cout << whatToDo << " " << TimeOffset  << " ns "<< endl;
      }
      //Thresholds
      else if (whatToDo == "SI_X_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        ThresholdSiX = atof(DataBuffer.c_str());
        ThresholdSiX = ThresholdSiX*keV;
        cout << whatToDo << " " << ThresholdSiX  << " MeV "<< endl;
      }
      else if (whatToDo == "SI_Y_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        ThresholdSiY = atof(DataBuffer.c_str());
        ThresholdSiY = ThresholdSiY*keV;
        cout << whatToDo << " " << ThresholdSiY  << " MeV "<< endl;
      }
      else if (whatToDo == "SILI_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        ThresholdSiLi = atof(DataBuffer.c_str());
        ThresholdSiLi = ThresholdSiLi*keV;
        cout << whatToDo << " " << ThresholdSiLi  << " MeV "<< endl;
      }
      else if (whatToDo == "CSI_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        ThresholdCsI = atof(DataBuffer.c_str());
        ThresholdCsI = ThresholdCsI*keV;
        cout << whatToDo << " " << ThresholdCsI  << " MeV "<< endl;
      } 
      else if (AnalysisConfigFile.eof()) ReadingStatus = false;
    }
  }
}


// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void MUST2Array::ConstructDetector(G4LogicalVolume* world) {
  G4RotationMatrix* MMrot    = NULL;
  G4ThreeVector     MMpos    = G4ThreeVector(0, 0, 0);
  G4ThreeVector     MMu      = G4ThreeVector(0, 0, 0);
  G4ThreeVector     MMv      = G4ThreeVector(0, 0, 0);
  G4ThreeVector     MMw      = G4ThreeVector(0, 0, 0);
  G4ThreeVector     MMCenter = G4ThreeVector(0, 0, 0);
  bool              Si       = true;
  bool              SiLi     = true;
  bool              CsI      = true;

  G4int NumberOfTelescope = m_DefinitionType.size();

  for (G4int i = 0; i < NumberOfTelescope; i++) {
    // By Point
    if (m_DefinitionType[i]) {
      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing CsI
      MMu = m_X128_Y1[i] - m_X1_Y1[i];
      MMu = MMu.unit();

      MMv = m_X1_Y128[i] - m_X1_Y1[i];
      MMv = MMv.unit();

      MMw = MMv.cross(MMu);
      // if (MMw.z() > 0)MMw = MMv.cross(MMu)  ;
      MMw = MMw.unit();

      MMCenter
        = (m_X1_Y1[i] + m_X1_Y128[i] + m_X128_Y1[i] + m_X128_Y128[i]) / 4;

      // Passage Matrix from Lab Referential to Telescope Referential
      MMrot = new G4RotationMatrix(MMv, MMu, MMw);
      MMpos = MMw * Length * 0.5 + MMCenter;
    }

    // By Angle
    else {
      G4double Theta = m_Theta[i];
      G4double Phi   = m_Phi[i];

      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      // Phi is angle between X axis and projection in (X,Y) plan
      // Theta is angle between  position vector and z axis
      G4double wX = m_R[i] * sin(Theta / rad) * cos(Phi / rad);
      G4double wY = m_R[i] * sin(Theta / rad) * sin(Phi / rad);
      G4double wZ = m_R[i] * cos(Theta / rad);
      MMw         = G4ThreeVector(wX, wY, wZ);

      // vector corresponding to the center of the module
      G4ThreeVector CT = MMw;

      // vector parallel to one axis of silicon plane
      G4double      ii = cos(Theta / rad) * cos(Phi / rad);
      G4double      jj = cos(Theta / rad) * sin(Phi / rad);
      G4double      kk = -sin(Theta / rad);
      G4ThreeVector Y  = G4ThreeVector(ii, jj, kk);

      MMw = MMw.unit();
      MMu = MMw.cross(Y);
      MMv = MMw.cross(MMu);
      MMv = MMv.unit();
      MMu = MMu.unit();

      // Passage Matrix from Lab Referential to Telescope Referential
      // MUST2
      MMrot = new G4RotationMatrix(MMu, MMv, MMw);
      // Telescope is rotate of Beta angle around MMv axis.
      MMrot->rotate(m_beta_u[i], MMu);
      MMrot->rotate(m_beta_v[i], MMv);
      MMrot->rotate(m_beta_w[i], MMw);
      // translation to place Telescope
      MMpos = MMw * Length * 0.5 + CT;
    }

    Si   = m_wSi[i];
    SiLi = m_wSiLi[i];
    CsI  = m_wCsI[i];

    VolumeMaker(i + 1, MMpos, MMrot, Si, SiLi, CsI, world);
  }

  delete MMrot;
}

// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method

void MUST2Array::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree*      pTree     = pAnalysis->GetTree();
  if (!pTree->FindBranch("MUST2")) {
    pTree->Branch("MUST2", "TMust2Data", &m_Event);
  }
  pTree->SetBranchAddress("MUST2", &m_Event);
}

// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void MUST2Array::ReadSensitive(const G4Event*) {
  m_Event->Clear();

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////// Used to Read Event Map of detector
  /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////

  /////////////////////
  // Read the Scorer associate to the Silicon Strip
  DSSDScorers::PS_Images* SiScorer
    = (DSSDScorers::PS_Images*)m_StripScorer->GetPrimitive(0);

  bool     SiScoredHit; // flag true if first stage scores a hit above threshold
  set<int> trig; // list of telescope that got a Si trigger
  unsigned int sizeFront = SiScorer->GetFrontMult();
  unsigned int sizeBack  = SiScorer->GetBackMult();

  // Check for double match Strip : 
  // rare case where a particle hit a strip and then an interstrip
  // since the map idex is build on pixel value, we end up with the same strip
  // fired twice, which is impossible in reality.
  std::map< unsigned int, std::pair<double,double> > mapFront;
  std::map< unsigned int, std::pair<double,double> >::iterator it;

  for (unsigned int i = 0; i < sizeFront; i++) {
    double energy      = SiScorer->GetEnergyFront(i);
    int    detectorNbr = SiScorer->GetDetectorFront(i);
    double time        = SiScorer->GetTimeFront(i);
    // Pixel value at interaction point
    unsigned int a, r, g, b;
    //  pixel
    SiScorer->GetARGBFront(i, a, r, g, b);
    if (r == 0) {
      mapFront[b+detectorNbr*1e6].first+=energy;
      mapFront[b+detectorNbr*1e6].second=time;
    } 

    else { // Interstrip X, keep maximum shared energy
      double rand = G4UniformRand();
      if (rand > 0.5) {
        double energyX = rand * energy;
          mapFront[b+detectorNbr*1e6].first+=energyX;
          mapFront[b+detectorNbr*1e6].second=time;
        }

      else {
          double energyX = (1 - rand) * energy;
          mapFront[g+detectorNbr*1e6].first+=energyX;
          mapFront[g+detectorNbr*1e6].second=time;
        }
      }
    }

  for(it=mapFront.begin();it!=mapFront.end();it++){
    double energyX = RandGauss::shoot(it->second.first, ResoStrip);
    double timeX = TimeOffset - RandGauss::shoot(it->second.second, ResoTimeMust);
    unsigned int strip = it->first-1000000*(it->first/1000000);
    unsigned int det   = it->first/1000000;
    if (energyX > ThresholdSiX) {
      trig.insert(det);
      SiScoredHit = true;
      m_Event->SetStripXE(det, strip ,
          NPL::EnergyToADC(energyX, 0, 63, 8192, 16384)); 
      m_Event->SetStripXT(det, strip ,
          NPL::EnergyToADC(timeX, 0, 1000, 8192, 16384));
    }
  }

  // Check for double match Strip : 
  // rare case where a particle hit a strip and then an interstrip
  // since the map idex is build on pixel value, we end up with the same strip
  // fired twice, which is impossible in reality.
  std::map< unsigned int, std::pair<double,double> > mapBack;

  for (unsigned int i = 0; i < sizeBack; i++) {
    double energy      = SiScorer->GetEnergyBack(i);
    int    detectorNbr = SiScorer->GetDetectorBack(i);
    double time        = SiScorer->GetTimeBack(i);

      // Pixel value at interaction point
      unsigned int a, r, g, b;
      //  pixel
      SiScorer->GetARGBBack(i, a, r, g, b);
      if (r == 0) {
          mapBack[b+detectorNbr*1e6].first+=energy;
          mapBack[b+detectorNbr*1e6].second=time;
      }
      else { // Interstrip Y, keep both strip with shared energy
        double rand     = G4UniformRand();
        double energyY1 = rand * energy;
          mapBack[b+detectorNbr*1e6].first+=energyY1;
          mapBack[b+detectorNbr*1e6].second=time;

        double energyY2 = (1 - rand) * energy;
        mapBack[g+detectorNbr*1e6].first+=energyY2;
        mapBack[g+detectorNbr*1e6].second=time;
        }
      }

  for(it=mapBack.begin();it!=mapBack.end();it++){
    double energyY = RandGauss::shoot(it->second.first, ResoStrip);
    double timeY = TimeOffset - RandGauss::shoot(it->second.second, ResoTimeMust);
    unsigned int strip = it->first-1000000*(it->first/1000000);
    unsigned int det   = it->first/1000000;
    if (energyY > ThresholdSiY) {
      trig.insert(det);
      SiScoredHit = true;
      m_Event->SetStripYE(det, strip ,
          NPL::EnergyToADC(energyY, 0, 63, 8192, 0)); 
      m_Event->SetStripYT(det, strip ,
          NPL::EnergyToADC(timeY, 0, 1000, 8192, 16384));
    }
  }


  // Look for 2nd and 3rd stage only if 1st stage is hit
  if (SiScoredHit) {
    // SiLi //
    CalorimeterScorers::PS_Calorimeter* SiLiScorer
      = (CalorimeterScorers::PS_Calorimeter*)m_SiLiScorer->GetPrimitive(0);

    unsigned int sizeSiLi = SiLiScorer->GetMult();
    for (unsigned int i = 0; i < sizeSiLi; i++) {
      double ESiLi = RandGauss::shoot(SiLiScorer->GetEnergy(i), ResoSiLi);
      vector<unsigned int> level = SiLiScorer->GetLevel(i);
      if(ESiLi>ThresholdSiLi){
        m_Event->SetSiLiE(level[0], level[1],
            NPL::EnergyToADC(ESiLi, 0, 250, 8192, 16384));
        double timeSiLi = RandGauss::shoot(SiLiScorer->GetTime(i), ResoTimeMust);
        m_Event->SetSiLiT(level[0], level[1],
            NPL::EnergyToADC(timeSiLi, 0, 1000, 16384, 8192));
      }
    }

    // CsI //
    CalorimeterScorers::PS_Calorimeter* CsIScorer
      = (CalorimeterScorers::PS_Calorimeter*)m_CsIScorer->GetPrimitive(0);

    unsigned int sizeCsI = CsIScorer->GetMult();
    for (unsigned int i = 0; i < sizeCsI; i++) {
      double ECsI = RandGauss::shoot(CsIScorer->GetEnergy(i), ResoCsI);
      vector<unsigned int> level = CsIScorer->GetLevel(i);
      if(ECsI>ThresholdCsI){
        m_Event->SetCsIE(level[0], level[1],
            NPL::EnergyToADC(ECsI, 0, 250, 8192, 16384));
        double timeCsI = RandGauss::shoot(CsIScorer->GetTime(i), ResoTimeMust);
        m_Event->SetCsIT(level[0], level[1],
            NPL::EnergyToADC(timeCsI, 0, 1000, 16384, 8192));
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MUST2Array::InitializeScorers() {
  //	Silicon Associate Scorer

  bool already_exist = false;
  m_StripScorer      = CheckScorer("MUST2_StripScorer", already_exist);
  m_SiLiScorer       = CheckScorer("MUST2_SiLiScorer", already_exist);
  m_CsIScorer        = CheckScorer("MUST2_CsIScorer", already_exist);

  // if the scorer were created previously nothing else need to be made
  if (already_exist)
    return;

  string              nptool   = getenv("NPTOOL");
  G4VPrimitiveScorer* SiScorer = new DSSDScorers::PS_Images(
      "SiScorer", nptool + "/NPLib/Detectors/MUST2/ressources/maskFront.png",
      nptool + "/NPLib/Detectors/MUST2/ressources/maskBack.png", 97.22/12800, 97.22/12800, 0,
      0, 0xffff0000, 0);

  G4VPrimitiveScorer* InterScorer
    = new InteractionScorers::PS_Interactions("SiScorer", ms_InterCoord, 0);

  // and register it to the multifunctionnal detector
  m_StripScorer->RegisterPrimitive(SiScorer);
  m_StripScorer->RegisterPrimitive(InterScorer);

  //	SiLi Associate Scorer
  vector<int>         SiLi_nesting = {3, 0};
  G4VPrimitiveScorer* SiLiScorer
    = new CalorimeterScorers::PS_Calorimeter("SiLiScorer", SiLi_nesting);
  m_SiLiScorer->RegisterPrimitive(SiLiScorer);

  //	CsI Associate Scorer
  vector<int>         CsI_nesting = {2, 0};
  G4VPrimitiveScorer* CsIScorer
    = new CalorimeterScorers::PS_Calorimeter("CsIScorer", CsI_nesting, 0);
  m_CsIScorer->RegisterPrimitive(CsIScorer);

  //	Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_StripScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SiLiScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_CsIScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MUST2Array::InitializeMaterial() {

  m_MaterialSilicon
    = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialAluminium
    = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialIron = MaterialManager::getInstance()->GetMaterialFromLibrary("Fe");
  m_MaterialVacuum
    = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4RotationMatrix* Rotation(double tetaX, double tetaY, double tetaZ) {
  double PI   = 3.141592653589793238;
  double radX = tetaX * PI / 180.;
  double radY = tetaY * PI / 180.;
  double radZ = tetaZ * PI / 180.;

  G4ThreeVector col1 = G4ThreeVector(cos(radZ) * cos(radY),
      -sin(radZ) * cos(radY), -sin(radY));
  G4ThreeVector col2
    = G4ThreeVector(sin(radZ) * cos(radX) - cos(radZ) * sin(radY) * sin(radX),
        cos(radZ) * cos(radX) + sin(radZ) * sin(radY) * sin(radX),
        -cos(radY) * sin(radX));
  G4ThreeVector col3
    = G4ThreeVector(sin(radZ) * sin(radX) + cos(radZ) * sin(radY) * sin(radX),
        cos(radZ) * sin(radX) - sin(radZ) * sin(radY) * cos(radX),
        cos(radY) * cos(radX));

  return (new G4RotationMatrix(col1, col2, col3));
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* MUST2Array::Construct() {
  return (NPS::VDetector*)new MUST2Array();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
  class proxy_nps_must2 {
    public:
      proxy_nps_must2() {
        NPS::DetectorFactory::getInstance()->AddToken("M2Telescope", "MUST2");
        NPS::DetectorFactory::getInstance()->AddDetector("M2Telescope",
            MUST2Array::Construct);
      }
  };

  proxy_nps_must2 p_nps_must2;
}
