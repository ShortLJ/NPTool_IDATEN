/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Goigoux  contact address: thomas.goigoux@cea.fr   *
 *                                                                           *
 * Creation Date  : july 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Exogam simulation                                   *
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
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Para.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

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
#include "Exogam.hh"
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Exogam_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 10*keV;
  //const double ResoTime = 4.5*ns ;  //not used
  const double ResoEnergy = 2.*keV ;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Exogam Specific Method
Exogam::Exogam(){
  m_Event = new TExogamData() ;
  m_ExogamScorer = 0;
  
  InitializeMaterials();

  HalfLengthCan = 7.35*cm;
  TaperLengthCan = 4.325*cm;
  distCollimatorToBGOSShield = 2.95*cm;

  rm90.rotateZ(90.*deg);
  rm90m.rotateZ(-90.*deg);
  rm180.rotateZ(180.*deg);
  rm270.rotateZ(270.*deg);
}

Exogam::~Exogam(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int Exogam::InitializeMaterials()
{
  m_Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  m_Aluminum = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_Copper = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
  m_Germanium = MaterialManager::getInstance()->GetMaterialFromLibrary("Ge");

  m_BGO = new G4Material("BGO", 7.13*g/cm3, 3, kStateSolid);  //BGO does not exist in nptool !!
  m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("Bi"),4);
  m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("Ge"),3);
  m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("O"),12);

  m_CsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Exogam::BuildClover(int i_clo, G4LogicalVolume* world)
{
  // enveloppe of the whole Clover (i.e. suppressed Clover) called 'SupClover' including: 
  //  the cryostat, dewar, side shield, back catcher, collimator
  G4double dzEnv = 40.472*cm;
  G4double dx1Env = 3.17*cm;
  G4double dy1Env = 3.17*cm;
  G4double dx2Env = 2.*dzEnv*tan(22.5*deg)+dx1Env;
  G4double dy2Env = 2.*dzEnv*tan(22.5*deg)+dy1Env;

  G4Trd* solidSupClover = new G4Trd("SupClover",dx1Env,dx2Env,dy1Env,dy2Env,dzEnv);
  G4LogicalVolume * logicSupClover = new G4LogicalVolume(solidSupClover,m_Vacuum,"SupClover"); 

  Offset=dzEnv;//-distCollimatorToGeCan;

  G4RotationMatrix rm;
  rm.rotateX(m_ThetaX[i_clo]/rad).rotateY(m_ThetaY[i_clo]/rad).rotateZ(m_ThetaZ[i_clo]/rad);

  new G4PVPlacement(G4Transform3D(rm,
		     G4ThreeVector(m_X[i_clo]*mm, m_Y[i_clo]*mm, m_Z[i_clo]*mm+Offset)),
	       logicSupClover,"Clover",world,false,i_clo+1,false);  //this void overlaps the whole setup

  // The Cryostat
  ////////////////
  // The Aluminum Clover can ( "CloverCan" )...
  //
  G4double PhiStartCan = 45.*deg;
  G4double PhiTotCan = 360.*deg;

  G4double zPlaneCan[3];
  G4double rInnerCan[3];

  G4double zPlaneVac[3];
  G4double rInnerVac[3];
  G4double rOuterVac[3];

  zPlaneCan[0] = -HalfLengthCan;
  zPlaneCan[1] = -HalfLengthCan+TaperLengthCan;
  zPlaneCan[2] =  HalfLengthCan;

  G4double rOuterCan[3];	// used to build the shield
  rOuterCan[0] = 4.4085*cm;
  rOuterCan[1] = 6.2*cm;
  rOuterCan[2] = 6.2*cm;

  rInnerCan[0] = rInnerCan[1] = rInnerCan[2] = 0.*cm;

  G4Polyhedra* solidCloverCan = new G4Polyhedra("CloverCan",PhiStartCan,PhiTotCan,4,3,
						zPlaneCan,rInnerCan,rOuterCan);

  G4LogicalVolume* logicCloverCan = new G4LogicalVolume(solidCloverCan,m_Aluminum,"CloverCan");

  // The position of the Clover can in the SupClover:
G4ThreeVector posClover(0.*cm,0.*cm,-Offset+HalfLengthCan+0.001*mm); //+0.001mm to avoid roundoff errors

new G4PVPlacement(0,posClover, logicCloverCan,"CloverCan",logicSupClover,false,i_clo+1,true); //There is an overlap with vacuum SupClover

// The vacuum clover ( "Vac" ) ...
//
G4double HalfLengthVac = 7.175*cm;
G4double TaperLengthVac = 4.0842*cm;

zPlaneVac[0] = -HalfLengthVac;
zPlaneVac[1] = -HalfLengthVac+TaperLengthVac;
zPlaneVac[2] =  HalfLengthVac;
rOuterVac[0] = 4.3083*cm;
rOuterVac[1] = 6.0*cm;
rOuterVac[2] = 6.0*cm;

rInnerVac[0] = rInnerVac[1] = rInnerVac[2] = 0.*cm;

G4Polyhedra* solidVac = new G4Polyhedra("Vac",PhiStartCan,PhiTotCan,4,3,
          zPlaneVac,rInnerVac,rOuterVac);
G4LogicalVolume * logicVac = new G4LogicalVolume(solidVac,m_Vacuum,"Vac");

G4ThreeVector positionVac = G4ThreeVector(0.*cm,0.*cm,-0.25*mm);
new G4PVPlacement(0,positionVac, logicVac,"Vac",logicCloverCan,false,i_clo+1,true);


//
// The enveloppe of the cold finger from the back side of the can to the Dewar
//

G4double zPlaneEnvColdFinger[6];
G4double rInnerEnvColdFinger[6];
G4double rOuterEnvColdFinger[6];

G4double PhiStart = 0.*deg;
G4double PhiTot = 360.*deg;
G4double EnvColdFingerHalfLength = 7.24*cm;

zPlaneEnvColdFinger[0] = -EnvColdFingerHalfLength;
zPlaneEnvColdFinger[1] = -EnvColdFingerHalfLength+4.1*cm;
zPlaneEnvColdFinger[2] = -EnvColdFingerHalfLength+4.1*cm;
zPlaneEnvColdFinger[3] = -EnvColdFingerHalfLength+4.9*cm;
zPlaneEnvColdFinger[4] = -EnvColdFingerHalfLength+4.9*cm;
zPlaneEnvColdFinger[5] =  EnvColdFingerHalfLength;

rInnerEnvColdFinger[0]=rInnerEnvColdFinger[1]=rInnerEnvColdFinger[2]=0.*cm;
rInnerEnvColdFinger[3]=rInnerEnvColdFinger[4]=rInnerEnvColdFinger[5]=0.*cm;

rOuterEnvColdFinger[0]=2.225*cm;
rOuterEnvColdFinger[1]=2.225*cm;
rOuterEnvColdFinger[2]=3.1*cm;
rOuterEnvColdFinger[3]=3.1*cm;
rOuterEnvColdFinger[4]=2.225*cm;
rOuterEnvColdFinger[5]=2.225*cm;

G4Polycone* solidEnvColdFinger = new G4Polycone("EnvColdFinger",PhiStart,PhiTot,6,
            zPlaneEnvColdFinger,rInnerEnvColdFinger,rOuterEnvColdFinger);

G4LogicalVolume* logicEnvColdFinger = new G4LogicalVolume(solidEnvColdFinger,m_Aluminum,"EnvColdFinger");

G4ThreeVector posEnvColdFinger = 
  G4ThreeVector(0.*cm,0.*cm,-Offset+2.*HalfLengthCan+EnvColdFingerHalfLength+0.005*mm);

  new G4PVPlacement(0,posEnvColdFinger,logicEnvColdFinger,"EnvColdFinger",logicSupClover,false,i_clo+1,true);

  // Its internal vacuum...
  G4double minRadiusIntEnvColdFinger = 0.*cm;
  G4double maxRadiusIntEnvColdFinger = 2.025*cm;
  G4double HalfLengthIntEnvColdFinger = 7.24*cm;
  G4double startPhiIntEnvColdFinger = 0.*deg;
  G4double deltaPhiIntEnvColdFinger = 360.*deg;

  G4Tubs* solidIntEnvColdFinger = new G4Tubs("IntDewar",minRadiusIntEnvColdFinger,maxRadiusIntEnvColdFinger,
					     HalfLengthIntEnvColdFinger,startPhiIntEnvColdFinger,deltaPhiIntEnvColdFinger);

  G4LogicalVolume* logicIntEnvColdFinger = new G4LogicalVolume(solidIntEnvColdFinger,m_Vacuum,"IntEnvColdFinger");

  // and its position in the cold finger enveloppe.
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicIntEnvColdFinger,"IntEnvColdFinger",logicEnvColdFinger,false,i_clo+1,true);

  
  // The cold finger and the associated plate
  //
  G4double xHalfLengthCFPlate = 5.04*cm;
  G4double yHalfLengthCFPlate = 5.04*cm;
  G4double zHalfLengthCFPlate = 1.*mm;

  G4Box* solidCFPlate = new G4Box("CFPlate",xHalfLengthCFPlate,yHalfLengthCFPlate,
				  zHalfLengthCFPlate);

  G4LogicalVolume* logicCFPlate = new G4LogicalVolume(solidCFPlate,m_Copper,"CFPlate");

  G4ThreeVector posCFPlate(0.*cm,0.*cm,-HalfLengthVac+9.65*cm); // 0.55(d(IntCan-Ge)
                                                                // +9.(Ge length)+0.1(half length plate)
  new G4PVPlacement(0,posCFPlate,logicCFPlate,"CFPlate",logicVac,false,i_clo+1,true);

  // The cold finger (internal part)
  //
  G4double minRadiusIntCF = 0.*cm;
  G4double maxRadiusIntCF = 1.5*cm;
  G4double HalfLengthIntCF = 2.30*cm;
  G4double startPhiIntCF = 0.*deg;
  G4double deltaPhiIntCF = 360.*deg;

  G4Tubs* solidIntCF = new G4Tubs("IntCF",minRadiusIntCF,maxRadiusIntCF,
				  HalfLengthIntCF,startPhiIntCF,deltaPhiIntCF);

  G4LogicalVolume* logicIntCF = 
    new G4LogicalVolume(solidIntCF,m_Copper,"IntCF");

  // its position vs CloverCan...
  G4ThreeVector posIntCF(0.*cm,0.*cm,4.875*cm); // -7.175 (halflengthcan internal)
						// +0.55 (ext Can - Ge)
						// +9.0 (Ge length)
						// +0.2 (CF plate)
						// +2.3 (IntCF length)

  new G4PVPlacement(0,posIntCF,logicIntCF,"IntCF",logicVac,false,i_clo+1,true);

  // The cold finger (external part)
  //
  G4double minRadiusExtCF = 0.*cm;
  G4double maxRadiusExtCF = 2.0*cm;
  G4double HalfLengthExtCF = 7.2*cm;
  G4double startPhiExtCF = 0.*deg;
  G4double deltaPhiExtCF = 360.*deg;

  G4Tubs* solidExtCF = new G4Tubs("IntCF",minRadiusExtCF,maxRadiusExtCF,
				  HalfLengthExtCF,startPhiExtCF,deltaPhiExtCF);

  G4LogicalVolume* logicExtCF = 
    new G4LogicalVolume(solidExtCF,m_Copper,"ExtCF");

  // its position vs EnvColdFinger...
  G4ThreeVector posExtCF(0.*cm,0.*cm,0.*cm); 
  new G4PVPlacement(0,posExtCF,logicExtCF,"ExtCF",logicIntEnvColdFinger,false,i_clo+1,true);

  // The Dewar
  //
  G4double minRadiusDewar = 0.*cm;
  G4double maxRadiusDewar = 10.9*cm;
  G4double HalfLengthDewar = 15.2*cm;
  G4double startPhiDewar = 0.*deg;
  G4double deltaPhiDewar = 360.*deg;

  G4Tubs* solidDewar = new G4Tubs("Dewar",minRadiusDewar,maxRadiusDewar,
				  HalfLengthDewar,startPhiDewar,deltaPhiDewar);

  G4LogicalVolume* logicDewar = new G4LogicalVolume(solidDewar,m_Aluminum,"Dewar");

  G4double distFrontToMidDewar = 
    -Offset+2.*(HalfLengthCan+EnvColdFingerHalfLength)+HalfLengthDewar+0.01*mm; 
  //+0.01mm to avoid roundoff errors

  G4ThreeVector posDewar = G4ThreeVector(0.*cm,0.*cm,distFrontToMidDewar);
  new G4PVPlacement(0,posDewar,logicDewar,"Dewar",logicSupClover,false,i_clo+1,true);

  /////////////////////////////////////////
  //  Construction of the active Ge volume:
  /////////////////////////////////////////
  //  A: Ge diode built from cuts subtracted from a cylinder (the "GeDiode")
  //
  G4double minRadiusGeDiode = 0.*cm;
  G4double maxRadiusGeDiode = 3.0*cm;
  G4double HalfLengthGeDiode = 4.5*cm;
  G4double startPhiGeDiode = 0.*deg;
  G4double deltaPhiGeDiode = 360.*deg;

  G4Tubs* solidGeDiode = new G4Tubs("GeDiode",minRadiusGeDiode,maxRadiusGeDiode,
				    HalfLengthGeDiode,startPhiGeDiode,deltaPhiGeDiode);
  //
  // External Tapered volume all along the diode ( "Cut1&2" )
  //
  //
  // Cut 1 :
  //
  G4double dummy = acos(2.9/3.0);
  G4double xHalfLengthCut1 = 0.5*mm;
  G4double yHalfLengthCut1 = 2.9*tan(dummy)*cm;
  G4double zHalfLengthCut1 = 4.55*cm;

  G4Box* solidCut1 = new G4Box("Cut1",xHalfLengthCut1,yHalfLengthCut1,
			       zHalfLengthCut1);

  //
  //... and its position vs GeDiode
  //

  G4ThreeVector transCut1(2.95*cm,0.*cm,0.*cm);
  G4SubtractionSolid* solidGeMinusCut1 = 
    new G4SubtractionSolid("GeMinusCut1",solidGeDiode,solidCut1,0,transCut1);

  G4ThreeVector transCut2(0.,2.95*cm,0.);
  G4Transform3D positionCut2(rm90,transCut2);

  G4SubtractionSolid* solidGeMinusCut12 = 
    new G4SubtractionSolid("GeMinusCut12",solidGeMinusCut1,solidCut1,positionCut2);
  //
  // External Tapered volume at the front face ( "Cut3&4" )

  G4double cosTap = cos(22.5*deg);  
  G4double sinTap = sin(22.5*deg);
  G4double tanTap = tan(22.5*deg);

  G4double xHalfLengthCut3 = 3.0*cm;
  G4double yHalfLengthCut3 = 1.5*cm*sinTap;
  G4double zHalfLengthCut3 = 1.5*cm/cosTap;

  G4Box* solidCut3 = new G4Box("Cut3",xHalfLengthCut3,yHalfLengthCut3,
			       zHalfLengthCut3+0.5*cm);
  
  G4double yCut3 = 2.9*cm-1.5*cm*tanTap+yHalfLengthCut3*cosTap;

  G4double temp = zHalfLengthCut3*cosTap-yHalfLengthCut3*sinTap;
  G4double zCut3 = -HalfLengthGeDiode+temp;

  G4RotationMatrix rmCut3;
  rmCut3.rotateX(-22.5*deg);

  G4ThreeVector transCut3(0.,yCut3,zCut3);
  G4Transform3D positionCut3(rmCut3,transCut3);

  G4SubtractionSolid* solidGeMinusCut123 = 
    new G4SubtractionSolid("GeMinusCut123",solidGeMinusCut12,solidCut3,positionCut3);

  G4Box* solidCut4 = new G4Box("Cut4",yHalfLengthCut3,xHalfLengthCut3,
			       zHalfLengthCut3);

  G4RotationMatrix rmCut4;
  rmCut4.rotateY(22.5*deg);

  G4ThreeVector transCut4(yCut3,0.,zCut3);
  G4Transform3D positionCut4(rmCut4,transCut4);

  G4SubtractionSolid* solidGeMinusCut1234 = 
    new G4SubtractionSolid("GeMinusCut1234",solidGeMinusCut123,solidCut4,positionCut4);

  dummy = acos(2.45/3.0);
  G4double xHalfLengthCut5 = 5.5*mm;
  G4double yHalfLengthCut5 = 2.45*tan(dummy)*cm;
  G4double zHalfLengthCut5 = 4.55*cm;

  G4Box* solidCut5 = new G4Box("Cut5",xHalfLengthCut5,yHalfLengthCut5,
			       zHalfLengthCut5);

  G4ThreeVector transCut5(-3.0*cm,0.*cm,0.*cm);

  G4SubtractionSolid* solidGeMinusCut12345 = 
    new G4SubtractionSolid("GeMinusCut12345",solidGeMinusCut1234,solidCut5,0,transCut5);

  G4ThreeVector transCut6(0.,-3.0*cm,0.);
  G4Transform3D positionCut6(rm90,transCut6);

  G4SubtractionSolid* solidGe = 
    new G4SubtractionSolid("Ge",solidGeMinusCut12345,solidCut5,positionCut6);

  // Now the individual diode is built; create logical volumes for each of
  // the four individual diodes A, B, C and D:

  G4LogicalVolume * logicGeA = new G4LogicalVolume(solidGe,m_Germanium,"GeA");
  G4LogicalVolume * logicGeB = new G4LogicalVolume(solidGe,m_Germanium,"GeB");
  G4LogicalVolume * logicGeC = new G4LogicalVolume(solidGe,m_Germanium,"GeC");
  G4LogicalVolume * logicGeD = new G4LogicalVolume(solidGe,m_Germanium,"GeD");

  logicGeA -> SetSensitiveDetector(m_ExogamScorer);
  logicGeB -> SetSensitiveDetector(m_ExogamScorer);
  logicGeC -> SetSensitiveDetector(m_ExogamScorer);
  logicGeD -> SetSensitiveDetector(m_ExogamScorer);

  // positioning the tapered partial diodes (A to D)
  // into the real vacuum of the can
  G4double HalfDistanceBetweenDiodes = 0.5*mm;

  G4double xDumVac = 2.45*cm+HalfDistanceBetweenDiodes; 
  G4double yDumVac = 2.45*cm+HalfDistanceBetweenDiodes;
  G4double zDumVac = -HalfLengthVac+5.05*cm; 	// 5.05 = 0.55 d(int can to Ge) +4.5(half length Ge) 

  G4ThreeVector positionVacA(xDumVac,yDumVac,zDumVac);

  G4ThreeVector posDumVacB(xDumVac,-yDumVac,zDumVac);
  G4Transform3D positionVacB(rm270,posDumVacB);

  G4ThreeVector posDumVacC(-xDumVac,-yDumVac,zDumVac);
  G4Transform3D positionVacC(rm180,posDumVacC);

  G4ThreeVector posDumVacD(-xDumVac,yDumVac,zDumVac);
  G4Transform3D positionVacD(rm90,posDumVacD);

  new G4PVPlacement(0,positionVacA,logicGeA,"GeA",logicVac,false,1,true);  //There is an overlap with vacumm Vac
  new G4PVPlacement(positionVacB,logicGeB,"GeB",logicVac,false,2,true);   
  new G4PVPlacement(positionVacC,logicGeC,"GeC",logicVac,false,3,true);   
  new G4PVPlacement(positionVacD,logicGeD,"GeD",logicVac, false,4, true);

  //
  // some material between the diodes to reproduce the experimental addback factor ...
  //

  G4double xAbsorb1 = 4.16*cm;
  G4double yAbsorb1 = 200.*um; // max = HalfDistanceBetweenDiodes = 0.5*mm;
  G4double zAbsorb1 = 4.5*cm;

  G4Box* solidAbsorb1 = new G4Box("Absorb1",xAbsorb1,yAbsorb1,zAbsorb1);

  G4double xAbsorb2 = 200*um; // max = HalfDistanceBetweenDiodes = 0.5*mm;
  G4double yAbsorb2 = 4.16*cm;
  G4double zAbsorb2 = 4.5*cm;

  G4Box* solidAbsorb2 = new G4Box("Absorb2",xAbsorb2,yAbsorb2,zAbsorb2);

  //G4UnionSolid* solidAbsorb = 
  //new G4UnionSolid("Absorb",solidAbsorb1,solidAbsorb2,0,0);
  G4UnionSolid* solidAbsorb = 
    new G4UnionSolid("Absorb",solidAbsorb1,solidAbsorb2);

  G4LogicalVolume* logicAbsorb = new G4LogicalVolume(solidAbsorb,m_Copper,"Absorb");

  G4ThreeVector positionAbsorb(0.,0.,zDumVac);

  new G4PVPlacement(0,positionAbsorb,logicAbsorb,"Absorb",logicVac,false,i_clo+1,true);

  //
  // Now: takes care of the holes and amorphous Ge in each diode:
  // Central hole with amorphous Ge for each diode. 
  //

  G4double minRadiusAGe1 = 0.*cm;
  G4double maxRadiusAGe1 = 0.52*cm;
  G4double HalfLengthAGe1 = 3.75*cm;
  G4double startPhiAGe1 = 0.*deg;
  G4double deltaPhiAGe1 = 360.*deg;

  //G4Tubs* solidAGe1 = new G4Tubs("AGe1",minRadiusAGe1,maxRadiusAGe1,
	//			 HalfLengthAGe1,startPhiAGe1,deltaPhiAGe1);

  //G4LogicalVolume* logicAGe1 = new G4LogicalVolume(solidAGe1,m_Germanium,"AGe1");
  
  // ... and second the hole in it:

  G4Tubs* solidHole1 = new G4Tubs("Hole1",minRadiusAGe1,maxRadiusAGe1-2.*mm,
				  HalfLengthAGe1,startPhiAGe1,deltaPhiAGe1);

  G4LogicalVolume* logicHole1 = new G4LogicalVolume(solidHole1,m_Vacuum,"Hole1");

  // Visu
  G4VisAttributes* CanVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5, 0.7));  // Grey
  G4VisAttributes* DewarVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));  // Grey
  
  logicCloverCan ->SetVisAttributes(CanVisAtt);
  logicEnvColdFinger ->SetVisAttributes(CanVisAtt);
  logicDewar ->SetVisAttributes(DewarVisAtt);
  logicSupClover->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* HoleVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,0.0));  // Black
  G4VisAttributes* AbsorbVisAtt= new G4VisAttributes(G4Colour(0.5,0.0,0.5,1));  // purple
  G4VisAttributes* GeAVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0, 0.6)); //Red
  G4VisAttributes* GeBVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0, 0.6)); //Green
  G4VisAttributes* GeCVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0, 0.6)); //Blue
  G4VisAttributes* GeDVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0, 0.6)); //White 

  logicHole1 ->SetVisAttributes(HoleVisAtt);
  logicAbsorb ->SetVisAttributes(AbsorbVisAtt);
  logicGeA ->SetVisAttributes(GeAVisAtt);
  logicGeB ->SetVisAttributes(GeBVisAtt);
  logicGeC ->SetVisAttributes(GeCVisAtt);
  logicGeD ->SetVisAttributes(GeDVisAtt);
  logicVac->SetVisAttributes (G4VisAttributes::Invisible);

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Exogam::AddDetector(double  X, double Y, double Z, double  ThetaX, double ThetaY, double ThetaZ){
  m_X.push_back(X);
  m_Y.push_back(Y);
  m_Z.push_back(Z);
  m_ThetaX.push_back(ThetaX);
  m_ThetaY.push_back(ThetaY);
  m_ThetaZ.push_back(ThetaZ);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Exogam::ReadConfiguration(NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Exogam");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> coord = {"X", "Y", "Z","ThetaX","ThetaY","ThetaZ"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel()) 
        cout << endl << "////  Exogam " << i+1 <<  endl;
      double X = blocks[i]->GetDouble("X","mm");
      double Y = blocks[i]->GetDouble("Y","mm");
      double Z = blocks[i]->GetDouble("Z","mm");
      double ThetaX = blocks[i]->GetDouble("ThetaX","deg");
      double ThetaY = blocks[i]->GetDouble("ThetaY","deg");
      double ThetaZ = blocks[i]->GetDouble("ThetaZ","deg");
      AddDetector(X,Y,Z,ThetaX, ThetaY, ThetaZ);
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
void Exogam::ConstructDetector(G4LogicalVolume* world){
  //G4double distBGOSShieldToGeCan = 3.2*cm;	 	// distance from the front face of the 
						// BGO Side Shield to the front face of the 
						// Ge can (theory: 3.2*cm)
  //G4double distCollimatorToGeCan=6.15*cm;		// distance from front face of the collimator
						// to the front face of the Ge can
  
  for ( unsigned i = 0; i < m_X.size(); ++i )
  {
    // Build and place Clover and its enveloppe
    BuildClover(i, world);
    
    //BuildSideCatcher();
    //BuildBackCatcher();

    //BuildSideShield();
    //BuildCollimator();
  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Exogam::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Exogam")){
    pTree->Branch("Exogam", "TExogamData", &m_Event) ;
  }
  pTree->SetBranchAddress("Exogam", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Exogam::ReadSensitive(const G4Event*){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_ExogamScorer->GetPrimitive(0);

    unsigned int size = Scorer->GetMult(); 
    for(unsigned int i = 0 ; i < size ; i++){
      double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Exogam_NS::ResoEnergy);
    if(Energy>Exogam_NS::EnergyThreshold){
      double Time = Scorer->GetTime(i);
      int CristalNbr = Scorer->GetLevel(i)[0];
      int CloverNbr = Scorer->GetLevel(i)[1];
      m_Event->SetCristal(CristalNbr);
      m_Event->SetClover(CloverNbr);
      m_Event->SetEnergy(Energy);
      m_Event->SetTime(Time);
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Exogam::InitializeScorers() { 
  bool already_exist = false; 
  m_ExogamScorer = CheckScorer("ExogamScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level({0, 1}); 
 
  m_ExogamScorer->RegisterPrimitive( 
          new CalorimeterScorers::PS_Calorimeter("Cristal",level, 0)); 
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ExogamScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Exogam::Construct(){
  return  (NPS::VDetector*) new Exogam();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Exogam{
    public:
      proxy_nps_Exogam(){
        NPS::DetectorFactory::getInstance()->AddToken("Exogam","Exogam");
        NPS::DetectorFactory::getInstance()->AddDetector("Exogam",Exogam::Construct);
      }
  };

  proxy_nps_Exogam p_nps_Exogam;
}
