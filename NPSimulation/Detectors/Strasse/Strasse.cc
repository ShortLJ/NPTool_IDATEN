/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta  contact address: matta@lpccaen.in2p3.fr        *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Strasse simulation                                  *
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
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4TwoVector.hh"

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
#include "Strasse.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

// CAD Mesh
#include "CADMesh.hh"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Strasse_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 10*keV;
  const double ResoEnergy = 0.015*MeV ;

  ////////////////////
  // Inner Detector //
  ////////////////////
  // Wafer parameter
  double Inner_Wafer_Length=100*mm;
  double Inner_Wafer_Width=50*mm;
  double Inner_Wafer_Thickness=300*micrometer;
  double Inner_Wafer_AlThickness=0.4*micrometer;
  double Inner_Wafer_PADExternal=1*cm;
  double Inner_Wafer_PADInternal=1*mm;
  double Inner_Wafer_GuardRing=0.5*mm;

  // PCB parameter
  double Inner_PCB_PortWidth=1*cm;
  double Inner_PCB_StarboardWidth=2*mm;
  double Inner_PCB_BevelAngle= 60*deg;
  double Inner_PCB_UpstreamWidth=1*cm;
  double Inner_PCB_DownstreamWidth=2*mm;
  double Inner_PCB_MidWidth=2*mm;
  double Inner_PCB_Thickness=3*mm;
  double Inner_PCB_Ledge = 1*mm ;
  double Inner_PCB_Step = 2*mm ;
  double Inner_Wafer_TransverseStrips= 128;
  double Inner_Wafer_LongitudinalStrips= 128;

  ////////////////////
  // Outer Detector //
  ////////////////////
  // Wafer parameter
  double Outer_Wafer_Length=150*mm;
  double Outer_Wafer_Width=75*mm;
  double Outer_Wafer_Thickness=300*micrometer;
  double Outer_Wafer_AlThickness=0.4*micrometer;
  double Outer_Wafer_PADExternal=1*cm;
  double Outer_Wafer_PADInternal=1*mm;
  double Outer_Wafer_GuardRing=0.5*mm;

  // PCB parameter
  double Outer_PCB_PortWidth=1*cm;
  double Outer_PCB_StarboardWidth=2*mm;
  double Outer_PCB_BevelAngle= 60*deg;
  double Outer_PCB_UpstreamWidth=1*cm;
  double Outer_PCB_DownstreamWidth=2*mm;
  double Outer_PCB_MidWidth=2*mm;
  double Outer_PCB_Thickness=3*mm;
  double Outer_PCB_Ledge = 1*mm ;
  double Outer_PCB_Step = 2*mm ;
  double Outer_Wafer_TransverseStrips= 128;
  double Outer_Wafer_LongitudinalStrips= 128;

  ////////////////////
  // Vacuum Chamber //
  ////////////////////
  double Chamber_Thickness= 3*mm;
  double Chamber_Cylinder_Length= 400*mm;
  double Chamber_Radius= 180*mm;
  double Chamber_ExitTube_Radius= 79.5*mm ;
  double Chamber_ExitTube_Length= 100*mm;
  double Chamber_Flange_Inner_Radius= 150*mm;
  double Chamber_Sphere_Radius= 220*mm ;
  double Chamber_Sphere_Shift= 60*mm;

}

using namespace Strasse_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Strasse Specific Method
Strasse::Strasse(){
  InitializeMaterial();
  m_Event = new TStrasseData() ;
  m_InnerScorer1 = 0;
  m_OuterScorer1 = 0;
  m_InnerScorer2 = 0;
  m_OuterScorer2 = 0;
  m_InnerDetector=0;
  m_OuterDetector=0;
  m_Chamber=0;
  m_Blades=0;
  m_Stars=0;
  m_Base=0;
  m_Electronic=0;
  found_chamber = false;
  found_stars = false;
  found_blades = false;
  found_base = false;
  // Dark Grey
  SiliconVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
  // Green
  PCBVisAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2)) ;
  // Gold Yellow
  PADVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.2)) ;
  // Light Grey
  StarsVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
  // Space transparent
  ChamberVisAtt = new G4VisAttributes(G4Colour(0.3, 0.4, 0.5,0.2)) ;
  // Light Blue
  GuardRingVisAtt = new G4VisAttributes(G4Colour(0.85, 0.85, 0.85,0.5)) ;
  // Light Blue
  BladeVisAtt = new G4VisAttributes(G4Colour(1, 0.65, 0.0,0.7)) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Strasse::~Strasse(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::AddInnerDetector(double  R, double  Z, double  Phi, double Shift, G4ThreeVector Ref){
  m_Inner_R.push_back(R);
  m_Inner_Z.push_back(Z);
  m_Inner_Phi.push_back(Phi);
  m_Inner_Shift.push_back(Shift);
  m_Inner_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::AddOuterDetector(double  R, double  Z, double  Phi, double Shift, G4ThreeVector Ref){
  m_Outer_R.push_back(R);
  m_Outer_Z.push_back(Z);
  m_Outer_Phi.push_back(Phi);
  m_Outer_Shift.push_back(Shift);
  m_Outer_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::AddChamber(double  Z){
  m_Chamber_Z.push_back(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildInnerDetector(){
  if(!m_InnerDetector){
    // Compute the needed full length of the PCB
    // along beam axis
    double Inner_PCB_Length= 2*Inner_Wafer_Length
      +Inner_PCB_UpstreamWidth
      +Inner_PCB_MidWidth
      +Inner_PCB_DownstreamWidth;

    // perpendicular to beam axis
    double Inner_PCB_Width= Inner_Wafer_Width
      +Inner_PCB_StarboardWidth
      +Inner_PCB_PortWidth;

    vector<G4TwoVector> PCBCrossSection;
    
    double l1;
    if(Inner_PCB_BevelAngle==90) l1 = 0;
    else l1 = Inner_PCB_Thickness*0.5/tan(Inner_PCB_BevelAngle);

    PCBCrossSection.push_back(G4TwoVector(Inner_PCB_Width*0.5-l1,-Inner_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(Inner_PCB_Width*0.5,Inner_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Inner_PCB_Width*0.5,Inner_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Inner_PCB_Width*0.5+l1,-Inner_PCB_Thickness*0.5));

    G4ExtrudedSolid* PCBFull =
      new G4ExtrudedSolid("PCBFull",
          PCBCrossSection,
          Inner_PCB_Length*0.5,// half length
          G4TwoVector(0,0),1,// offset, scale
          G4TwoVector(0,0),1);// offset, scale

    // Master Volume that encompass everything else
    m_InnerDetector =
      new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector", 0, 0, 0);
    m_InnerDetector->SetVisAttributes(G4VisAttributes::Invisible);

    ///////////////////////////////////////////////////////////////////////////
    // Build the external PCB frame
    // Calculate the hole shift within the PCB
    double Width_Shift= -0.5*Inner_PCB_Width + 0.5*Inner_Wafer_Width // Flush to border
      +Inner_PCB_PortWidth; // add the port side shift

    double Length_Shift1 = -0.5*Inner_PCB_Length + 0.5*Inner_Wafer_Length // Flush to border
      + Inner_PCB_UpstreamWidth;// add Upstream side shift

    double Length_Shift2 = Length_Shift1 // overlap detector 1
      + Inner_Wafer_Length // at opposing edge
      + Inner_PCB_MidWidth; // after mid width

    G4ThreeVector HoleShift1 = G4ThreeVector(Width_Shift, 0, Length_Shift1);
    G4ThreeVector HoleShift2 = G4ThreeVector(Width_Shift, 0, Length_Shift2);

    G4Box*  HoleShape = new G4Box("HoleShape",
        Inner_Wafer_Width*0.5,
        Inner_PCB_Thickness*0.5+0.1*mm,
        Inner_Wafer_Length*0.5);

    // Substracting the hole Shape from the Stock PCB
    G4SubtractionSolid* PCB_1 = new G4SubtractionSolid("PCB_1", PCBFull, HoleShape,
        new G4RotationMatrix,HoleShift1);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_1, HoleShape,
        new G4RotationMatrix,HoleShift2);

    // Sub Volume PCB
    G4LogicalVolume* logicPCB =
      new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(PCBVisAtt);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0),
        logicPCB,"Strasse_Inner_PCB",m_InnerDetector,
        false,0);

    ///////////////////////////////////////////////////////////////////////////
    // Build the PCB Ledge on wich Silicon is glued
    double Inner_PCB2_Thickness = Inner_PCB_Step; //Step size 
    double offsetPCB2 = Inner_PCB2_Thickness - Inner_PCB_Thickness; 

    double Inner_PCB2_MidWidth = Inner_PCB_MidWidth; 

    // perpendicular to beam axis
    double Inner_PCB2_Width= Inner_Wafer_Width;

    vector<G4TwoVector> PCB2CrossSection;
    PCB2CrossSection.push_back(G4TwoVector(Inner_PCB2_Width*0.5,-Inner_PCB2_Thickness*0.5));
    PCB2CrossSection.push_back(G4TwoVector(Inner_PCB2_Width*0.5,Inner_PCB2_Thickness*0.5));
    PCB2CrossSection.push_back(G4TwoVector(-Inner_PCB2_Width*0.5,Inner_PCB2_Thickness*0.5));
    PCB2CrossSection.push_back(G4TwoVector(-Inner_PCB2_Width*0.5,-Inner_PCB2_Thickness*0.5));

    //double Inner_PCB2_Length= Inner_PCB_Length;
    double Inner_PCB2_Length= 2*Inner_Wafer_Length + Inner_PCB_MidWidth;

    G4ExtrudedSolid* PCB2Full =
      new G4ExtrudedSolid("PCB2Full",
          PCB2CrossSection,
          Inner_PCB2_Length*0.5,// half length
          G4TwoVector(0,0),1,// offset, scale
          G4TwoVector(0,0),1);// offset, scale


    double Length_Shift21 = -0.5*Inner_PCB_Length  // Flush to border
                           + 0.5*(Inner_PCB_UpstreamWidth+Inner_PCB_DownstreamWidth) // add Upstream side shift
                           + 0.5*Inner_Wafer_Length;
    double Length_Shift22 = Length_Shift21 // overlap detector 1
      + Inner_Wafer_Length // at opposing edge
      + Inner_PCB_MidWidth; // after mid width

    G4ThreeVector HoleShift21 = G4ThreeVector(0, 0, Length_Shift21);
    G4ThreeVector HoleShift22 = G4ThreeVector(0, 0, Length_Shift22);

    G4Box* HoleShape2 = new G4Box("HoleShape2",
        Inner_Wafer_Width*0.5 - Inner_PCB_Ledge,
        Inner_PCB2_Thickness,
        Inner_Wafer_Length*0.5 - Inner_PCB_Ledge);

    // Substracting the hole Shape from the Stock PCB
    G4SubtractionSolid* PCB2_1 = new G4SubtractionSolid("PCB2_1", PCB2Full, HoleShape2,
        new G4RotationMatrix,HoleShift21);
    G4SubtractionSolid* PCB2_2 = new G4SubtractionSolid("PCB2_2", PCB2_1, HoleShape2,
        new G4RotationMatrix,HoleShift22);

    G4ThreeVector HoleCenterBar = G4ThreeVector(0, 0, 0);
    G4Box* HoleShapeCenterBar = new G4Box("HoleShapeCenterBar",
        Inner_PCB2_Width*0.5+0.1,
        Inner_PCB2_Thickness,
        Inner_PCB2_MidWidth*0.5);

    G4SubtractionSolid* PCB2_3 = new G4SubtractionSolid("PCB2_3", PCB2_2, HoleShapeCenterBar,
        new G4RotationMatrix,HoleCenterBar);

    // Sub Volume PCB
    G4LogicalVolume* logicPCB2 =
      new G4LogicalVolume(PCB2_3,m_MaterialPCB,"logicPCB2", 0, 0, 0);
    logicPCB2->SetVisAttributes(PADVisAtt);

    // Offset along beam axis between PCB middle and (2*Wafer+MiddleBar) volume center
    double CentralZOffset = - Inner_PCB_Length*0.5
                    + Inner_PCB_UpstreamWidth
                    + Inner_Wafer_Length
                    + Inner_PCB_MidWidth*0.5;

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,-0.5*offsetPCB2,CentralZOffset),
        logicPCB2,"Strasse_Inner_PCB2",m_InnerDetector,
        false,0);


    /////////////////////////////////////////////////////////////////////////// 
    // Build the Wafer
    // Sub volume Wafer
    G4Box*  WaferShape = new G4Box("WaferShape",
        Inner_Wafer_Width*0.5,
        Inner_Wafer_Thickness*0.5+Inner_Wafer_AlThickness,
        Inner_Wafer_Length*0.5);

    G4LogicalVolume* logicWafer1 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer1", 0, 0, 0);
    logicWafer1->SetVisAttributes(GuardRingVisAtt);

    G4LogicalVolume* logicWafer2 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer2", 0, 0, 0);
    logicWafer2->SetVisAttributes(GuardRingVisAtt);

    // Shift along Y to flush the wafer to the pcb ledge on one side
    G4ThreeVector WaferShiftY = G4ThreeVector(0,-0.5*Inner_Wafer_Thickness
          -Inner_Wafer_AlThickness
          -0.5*Inner_PCB_Thickness
          -offsetPCB2-0.05,0);

    //G4ThreeVector WaferShiftZ = G4ThreeVector(0,0,-Inner_PCB_UpstreamWidth-Inner_PCB_MidWidth);
    G4ThreeVector WaferShiftZ = G4ThreeVector(0,0,CentralZOffset);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        WaferShiftY+WaferShiftZ // Shift along Y
        +HoleShift21, // Shift along Z to putwafer in the 1st hole 
        logicWafer1,"Strasse_Inner_Wafer1",m_InnerDetector,
        false,0);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        WaferShiftY+WaferShiftZ// Shift along Y
        +HoleShift22, // Shift along Z to put wafer in the 2nd hole 
        logicWafer2,"Strasse_Inner_Wafer2",m_InnerDetector,
        false,0);

    // Sub volume Active Wafer
    G4Box*  ActiveWaferShape = new G4Box("InnerActiveWaferShape",
        0.5*m_Active_InnerWafer_Width,
        0.5*Inner_Wafer_Thickness,
        0.5*m_Active_InnerWafer_Length);

    G4LogicalVolume* logicActiveWafer1 =
      new G4LogicalVolume(ActiveWaferShape,m_MaterialSilicon,"logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(SiliconVisAtt);
    logicActiveWafer1->SetSensitiveDetector(m_InnerScorer1);

    G4LogicalVolume* logicActiveWafer2 =
      new G4LogicalVolume(ActiveWaferShape,m_MaterialSilicon,"logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(SiliconVisAtt);
    logicActiveWafer2->SetSensitiveDetector(m_InnerScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0.5*(Inner_Wafer_PADExternal-Inner_Wafer_PADInternal)), // assymetric pading for bounding
        logicActiveWafer1,"Strasse_Inner_ActiveWafer1",logicWafer1,
        false,1);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,-0.5*(Inner_Wafer_PADExternal-Inner_Wafer_PADInternal)), // assymetric pading for bounding
        logicActiveWafer2,"Strasse_Inner_ActiveWafer2",logicWafer2,
        false,2);

  }

  return m_InnerDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildOuterDetector(){
  if(!m_OuterDetector){
    // Compute the needed full length of the PCB
    // along beam axis
    double Outer_PCB_Length= 2*Outer_Wafer_Length
      +Outer_PCB_UpstreamWidth
      +Outer_PCB_MidWidth
      +Outer_PCB_DownstreamWidth;

    // perpendicular to beam axis
    double Outer_PCB_Width= Outer_Wafer_Width
      +Outer_PCB_StarboardWidth
      +Outer_PCB_PortWidth;


    vector<G4TwoVector> PCBCrossSection;
    double l1;
    if(Outer_PCB_BevelAngle==90) l1 = 0;
    else l1 = Outer_PCB_Thickness*0.5/tan(Outer_PCB_BevelAngle);

    PCBCrossSection.push_back(G4TwoVector(Outer_PCB_Width*0.5-l1,-Outer_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(Outer_PCB_Width*0.5,Outer_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Outer_PCB_Width*0.5,Outer_PCB_Thickness*0.5));
    PCBCrossSection.push_back(G4TwoVector(-Outer_PCB_Width*0.5+l1,-Outer_PCB_Thickness*0.5));

    G4ExtrudedSolid* PCBFull =
      new G4ExtrudedSolid("PCBFull",
          PCBCrossSection,
          Outer_PCB_Length*0.5,// half length
          G4TwoVector(0,0),1,// offset, scale
          G4TwoVector(0,0),1);// offset, scale

    // Master Volume that encompass everything else
    m_OuterDetector =
      new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector", 0, 0, 0);
    m_OuterDetector->SetVisAttributes(G4VisAttributes::Invisible);

    // Build the PCB
    // Calculate the hole shift within the PCB
    double Width_Shift= -0.5*Outer_PCB_Width + 0.5*Outer_Wafer_Width // Flush to border
      +Outer_PCB_PortWidth; // add the port side shift

    double Length_Shift1 = -0.5*Outer_PCB_Length + 0.5*Outer_Wafer_Length // Flush to border
      + Outer_PCB_UpstreamWidth;// add Upstream side shift

    double Length_Shift2 = Length_Shift1 // overlap detector 1
      + Outer_Wafer_Length // at opposing edge
      + Outer_PCB_MidWidth; // after mid width

    G4ThreeVector HoleShift1 = G4ThreeVector(Width_Shift, 0, Length_Shift1);
    G4ThreeVector HoleShift2 = G4ThreeVector(Width_Shift, 0, Length_Shift2);

    G4Box*  HoleShape = new G4Box("HoleShape",
        Outer_Wafer_Width*0.5,
        Outer_PCB_Thickness*0.5+0.1*mm,
        Outer_Wafer_Length*0.5);

    // Substracting the hole Shape from the Stock PCB
    G4SubtractionSolid* PCB_1 = new G4SubtractionSolid("PCB_1", PCBFull, HoleShape,
        new G4RotationMatrix,HoleShift1);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_1, HoleShape,
        new G4RotationMatrix,HoleShift2);

    // Sub Volume PCB
    G4LogicalVolume* logicPCB =
      new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(PCBVisAtt);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0),
        logicPCB,"Strasse_Outer_PCB",m_OuterDetector,
        false,0);

    ///////////////////////////////////////////////////////////////////////////
    // Build the internal PCB layer
    double Outer_PCB2_Thickness = Outer_PCB_Step;//Step size
    double offsetPCB2 = Outer_PCB2_Thickness - Outer_PCB_Thickness;

    double Outer_PCB2_MidWidth = Outer_PCB_MidWidth; 

    // perpendicular to beam axis
    double Outer_PCB2_Width= Outer_Wafer_Width;

    vector<G4TwoVector> PCB2CrossSection;
    PCB2CrossSection.push_back(G4TwoVector(Outer_PCB2_Width*0.5,-Outer_PCB2_Thickness*0.5));
    PCB2CrossSection.push_back(G4TwoVector(Outer_PCB2_Width*0.5,Outer_PCB2_Thickness*0.5));
    PCB2CrossSection.push_back(G4TwoVector(-Outer_PCB2_Width*0.5,Outer_PCB2_Thickness*0.5));
    PCB2CrossSection.push_back(G4TwoVector(-Outer_PCB2_Width*0.5,-Outer_PCB2_Thickness*0.5));


    double Outer_PCB2_Length= 2*Outer_Wafer_Length + Outer_PCB_MidWidth;

    G4ExtrudedSolid* PCB2Full =
      new G4ExtrudedSolid("PCB2Full",
          PCB2CrossSection,
          Outer_PCB2_Length*0.5,// half length
          G4TwoVector(0,0),1,// offset, scale
          G4TwoVector(0,0),1);// offset, scale


    double Length_Shift21 = -0.5*Outer_PCB_Length  // Flush to border
                           + 0.5*(Outer_PCB_UpstreamWidth+Outer_PCB_DownstreamWidth) 
                           + 0.5*Outer_Wafer_Length;
    double Length_Shift22 = Length_Shift21 // overlap detector 1
      + Outer_Wafer_Length // at opposing edge
      + Outer_PCB_MidWidth; // after mid width

    G4ThreeVector HoleShift21 = G4ThreeVector(0, 0, Length_Shift21);
    G4ThreeVector HoleShift22 = G4ThreeVector(0, 0, Length_Shift22);

    G4Box* HoleShape2 = new G4Box("HoleShape2",
        Outer_Wafer_Width*0.5 - Outer_PCB_Ledge,
        Outer_PCB2_Thickness,
        Outer_Wafer_Length*0.5 - Outer_PCB_Ledge);

    // Substracting the hole Shape from the Stock PCB
    G4SubtractionSolid* PCB2_1 = new G4SubtractionSolid("PCB2_1", PCB2Full, HoleShape2,
        new G4RotationMatrix,HoleShift21);
    G4SubtractionSolid* PCB2_2 = new G4SubtractionSolid("PCB2_2", PCB2_1, HoleShape2,
        new G4RotationMatrix,HoleShift22);
      
    G4ThreeVector HoleCenterBar = G4ThreeVector(0, 0, 0);
    G4Box* HoleShapeCenterBar = new G4Box("HoleShapeCenterBar",
        Outer_PCB2_Width*0.5+0.1,
        Outer_PCB2_Thickness,
        Outer_PCB2_MidWidth*0.5);

    G4SubtractionSolid* PCB2_3 = new G4SubtractionSolid("PCB2_3", PCB2_2, HoleShapeCenterBar,
        new G4RotationMatrix,HoleCenterBar);

    // Sub Volume PCB
    G4LogicalVolume* logicPCB2 =
      new G4LogicalVolume(PCB2_3,m_MaterialPCB,"logicPCB2", 0, 0, 0);
    logicPCB2->SetVisAttributes(PADVisAtt);

    // Offset along beam axis between PCB middle and (2*Wafer+MiddleBar) volume center
    double CentralZOffset = - Outer_PCB_Length*0.5
                    + Outer_PCB_UpstreamWidth
                    + Outer_Wafer_Length
                    + Outer_PCB_MidWidth*0.5;

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,-0.5*offsetPCB2,CentralZOffset),
        logicPCB2,"Strasse_Outer_PCB2",m_OuterDetector,
        false,0);

    //////////////////////////////////////////////////////////////////
    // Build the Wafer
    // Sub volume Wafer
    G4Box*  WaferShape = new G4Box("WaferShape",
        Outer_Wafer_Width*0.5,
        Outer_Wafer_Thickness*0.5+Outer_Wafer_AlThickness,
        Outer_Wafer_Length*0.5);

    G4LogicalVolume* logicWafer1 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer1", 0, 0, 0);
    logicWafer1->SetVisAttributes(GuardRingVisAtt);

    G4LogicalVolume* logicWafer2 =
      new G4LogicalVolume(WaferShape,m_MaterialSilicon,"logicWafer2", 0, 0, 0);
    logicWafer2->SetVisAttributes(GuardRingVisAtt);

    // Shift along Y to flush the wafer to the pcb ledge on one side
    G4ThreeVector WaferShiftY = G4ThreeVector(0,-0.5*Outer_Wafer_Thickness
                    -Outer_Wafer_AlThickness
                    -0.5*Outer_PCB_Thickness
                    -offsetPCB2-0.05,0);

    //G4ThreeVector WaferShiftZ = G4ThreeVector(0,0,-Outer_PCB_UpstreamWidth-Outer_PCB_MidWidth);
    G4ThreeVector WaferShiftZ = G4ThreeVector(0,0,CentralZOffset);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        WaferShiftY+WaferShiftZ
        +HoleShift21, // Shift along Z to put wafer in the 1st hole 
        logicWafer1,"Strasse_Outer_Wafer1",m_OuterDetector,
        false,0);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        WaferShiftY+WaferShiftZ
        +HoleShift22, // Shift along Z to put wafer in the 1st hole 
        logicWafer2,"Strasse_Outer_Wafer2",m_OuterDetector,
        false,0);

    // Sub volume Active Wafer
    G4Box*  ActiveWaferShape = new G4Box("OuterActiveWaferShape",
        0.5*m_Active_OuterWafer_Width,
        0.5*Outer_Wafer_Thickness,
        0.5*m_Active_OuterWafer_Length);

    G4LogicalVolume* logicActiveWafer1 =
      new G4LogicalVolume(ActiveWaferShape,m_MaterialSilicon,"logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(SiliconVisAtt);
    logicActiveWafer1->SetSensitiveDetector(m_OuterScorer1);

    G4LogicalVolume* logicActiveWafer2 =
      new G4LogicalVolume(ActiveWaferShape,m_MaterialSilicon,"logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(SiliconVisAtt);
    logicActiveWafer2->SetSensitiveDetector(m_OuterScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,0.5*(Outer_Wafer_PADExternal-Outer_Wafer_PADInternal)),
        logicActiveWafer1,"Strasse_Outer_ActiveWafer1",logicWafer1,
        false,1);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),
        G4ThreeVector(0,0,-0.5*(Outer_Wafer_PADExternal-Outer_Wafer_PADInternal)),
        logicActiveWafer2,"Strasse_Outer_ActiveWafer2",logicWafer2,
        false,2);

  }
  return m_OuterDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildChamber(){
  if(!m_Chamber){
    // Needed Element
    G4Material* Material = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4RotationMatrix* Rot = new G4RotationMatrix();

    // Main Cylinder
    G4Tubs* Cylinder = new G4Tubs("StrasseCylinderVolume",
        Chamber_Radius-Chamber_Thickness,
        Chamber_Radius,Chamber_Cylinder_Length*0.5,
        0,360*deg);
    // for substraction
    G4Tubs* DummyCyl = new G4Tubs("StrasseDummyCylVolume",
        0,
        Chamber_Sphere_Radius*1.1,Chamber_Cylinder_Length*0.5,
        0,360*deg);


  //  G4LogicalVolume* ChamberCyl = new G4LogicalVolume(Cyl,Material,"logic_Strasse_Chamber",0,0,0);

    // Entrance Flange
    G4Tubs* Flange = new G4Tubs("StrasseFlangeVolume",
        Chamber_Flange_Inner_Radius,
        Chamber_Radius,1*cm,
        0,360*deg);

   // G4LogicalVolume* ChamberFlange = new G4LogicalVolume(Flange,Material,"logic_Strasse_Flange",0,0,0);

    // Spherial Portion
    G4Sphere* Sphere= new G4Sphere("StrasseSphereVolume",
        Chamber_Sphere_Radius-Chamber_Thickness,
        Chamber_Sphere_Radius,
        0,360*deg,
        0,180*deg);
    
    // Exit tube portion
    G4Tubs* Tube = new G4Tubs("StrasseTubeVolume",
        Chamber_ExitTube_Radius-Chamber_Thickness,
        Chamber_ExitTube_Radius,Chamber_ExitTube_Length*0.5,
        0,360*deg);
    G4Tubs* DummyTube = new G4Tubs("StrasseDummyTubeVolume",
        0,
        Chamber_ExitTube_Radius*0.99,Chamber_ExitTube_Length*0.5,
        0,360*deg);
    
    //Partial Sphere
    
    G4SubtractionSolid* Sphere1= new G4SubtractionSolid("Sphere1",Sphere,DummyCyl,
      Rot,G4ThreeVector(0,0,-Chamber_Sphere_Shift));
    G4SubtractionSolid* Sphere2= new G4SubtractionSolid("Sphere2",Sphere1,DummyTube,
      Rot,G4ThreeVector(0,0,Chamber_Sphere_Radius+Chamber_ExitTube_Length*0.5-2*cm));
    
    // Building the whole chamber
    G4UnionSolid* Chamber1= new G4UnionSolid("Chamber1",Cylinder,Flange,
      Rot,G4ThreeVector(0,0,-Chamber_Cylinder_Length*0.5));

    G4UnionSolid* Chamber2= new G4UnionSolid("Chamber2",Chamber1,Sphere2,
      Rot,G4ThreeVector(0,0,Chamber_Sphere_Shift));

    G4UnionSolid* Chamber3= new G4UnionSolid("Chamber3",Chamber2,Tube,
      Rot,G4ThreeVector(0,0,Chamber_Sphere_Shift+Chamber_Sphere_Radius+Chamber_ExitTube_Length*0.5-2*cm));

    m_Chamber = new G4LogicalVolume(Chamber3,Material,"logic_Strasse_Chamber",0,0,0);


    m_Chamber->SetVisAttributes(ChamberVisAtt);
  }


  return m_Chamber;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildChamberFromCAD(string path){
    if(!m_Chamber){
        auto mesh = CADMesh::TessellatedMesh::FromSTL((char*) path.c_str());
        mesh->SetScale(mm);
        //mesh->SetOffset(offset);
        mesh->SetReverse(false);

        auto cad_solid = mesh->GetSolid();
        m_Chamber = new G4LogicalVolume(cad_solid,
            m_MaterialAl,
            "Strasse_Chamber", 0, 0, 0);

       m_Chamber->SetVisAttributes(ChamberVisAtt);
    }
    return m_Chamber;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildBlades(string path){
    if(!m_Blades){
        auto mesh = CADMesh::TessellatedMesh::FromSTL((char*) path.c_str());
        mesh->SetScale(mm);
        //mesh->SetOffset(offset);
        mesh->SetReverse(false);

        auto cad_solid = mesh->GetSolid();
        m_Blades = new G4LogicalVolume(cad_solid,
            m_MaterialCu,
            "Strasse_Blades", 0, 0, 0);

       m_Blades->SetVisAttributes(BladeVisAtt);
    }
    return m_Blades;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildStars(string path){
    if(!m_Stars){
        auto mesh = CADMesh::TessellatedMesh::FromSTL((char*) path.c_str());
        mesh->SetScale(mm);
        //mesh->SetOffset(offset);
        mesh->SetReverse(false);

        auto cad_solid = mesh->GetSolid();
        m_Stars = new G4LogicalVolume(cad_solid,
            m_MaterialAl,
            "Strasse_Stars", 0, 0, 0);

       m_Stars->SetVisAttributes(StarsVisAtt);
    }
    return m_Stars;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Strasse::BuildBase(string path){
    if(!m_Base){
        auto mesh = CADMesh::TessellatedMesh::FromSTL((char*) path.c_str());
        mesh->SetScale(mm);
        //mesh->SetOffset(offset);
        mesh->SetReverse(false);

        auto cad_solid = mesh->GetSolid();
        m_Base = new G4LogicalVolume(cad_solid,
            m_MaterialAl,
            "Strasse_Base", 0, 0, 0);

       m_Base->SetVisAttributes(StarsVisAtt);
    }
    return m_Base;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Strasse::ReadConfiguration(NPL::InputParser parser){
  // Info block
  vector<NPL::InputBlock*> blocks_info = parser.GetAllBlocksWithTokenAndValue("Strasse","Info");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_info.size() << " info block founds " << endl; 

  if(blocks_info.size()>1){
    cout << "ERROR: can only accepte one info block, " << blocks_info.size() << " info block founds." << endl; 
    exit(1); 
  }

  vector<string> info = {
    "Inner_Wafer_Length",         
    "Inner_Wafer_Width",          
    "Inner_Wafer_Thickness",     
    "Inner_Wafer_AlThickness",    
    "Inner_Wafer_PADExternal",    
    "Inner_Wafer_PADInternal",  
    "Inner_Wafer_GuardRing",    
    "Inner_PCB_PortWidth",      
    "Inner_PCB_StarboardWidth", 
    "Inner_PCB_BevelAngle",     
    "Inner_PCB_UpstreamWidth",  
    "Inner_PCB_DownstreamWidth",
    "Inner_PCB_MidWidth",       
    "Inner_PCB_Thickness",      
    "Inner_PCB_Ledge",      
    "Inner_PCB_Step",      
    "Inner_Wafer_TransverseStrips",
    "Inner_Wafer_LongitudinalStrips",
    "Outer_Wafer_Length",       
    "Outer_Wafer_Width",        
    "Outer_Wafer_Thickness",    
    "Outer_Wafer_AlThickness",  
    "Outer_Wafer_PADExternal",  
    "Outer_Wafer_PADInternal",  
    "Outer_Wafer_GuardRing",    
    "Outer_PCB_PortWidth",      
    "Outer_PCB_StarboardWidth", 
    "Outer_PCB_BevelAngle",     
    "Outer_PCB_UpstreamWidth",  
    "Outer_PCB_DownstreamWidth",
    "Outer_PCB_MidWidth",       
    "Outer_PCB_Thickness",      
    "Outer_PCB_Ledge",      
    "Outer_PCB_Step",      
    "Outer_Wafer_TransverseStrips",
    "Outer_Wafer_LongitudinalStrips",
    "Chamber_Thickness",
    "Chamber_Cylinder_Length",
    "Chamber_Radius",
    "Chamber_ExitTube_Radius",
    "Chamber_ExitTube_Length",
    "Chamber_Flange_Inner_Radius",
    "Chamber_Sphere_Radius",
    "Chamber_Sphere_Shift"
  };

  if(blocks_info[0]->HasTokenList(info)){
    cout << endl << "////  Strasse info block" <<  endl;
    Inner_Wafer_Length = blocks_info[0]->GetDouble("Inner_Wafer_Length","mm");
    Inner_Wafer_Width = blocks_info[0]->GetDouble("Inner_Wafer_Width","mm");          
    Inner_Wafer_Thickness = blocks_info[0]->GetDouble("Inner_Wafer_Thickness","micrometer");      
    Inner_Wafer_AlThickness = blocks_info[0]->GetDouble("Inner_Wafer_AlThickness","micrometer");     
    Inner_Wafer_PADExternal = blocks_info[0]->GetDouble("Inner_Wafer_PADExternal","mm");     
    Inner_Wafer_PADInternal = blocks_info[0]->GetDouble("Inner_Wafer_PADInternal","mm");   
    Inner_Wafer_GuardRing = blocks_info[0]->GetDouble("Inner_Wafer_GuardRing","mm");     
    Inner_Wafer_TransverseStrips = blocks_info[0]->GetInt("Inner_Wafer_TransverseStrips");        
    Inner_Wafer_LongitudinalStrips = blocks_info[0]->GetInt("Inner_Wafer_LongitudinalStrips");       
    Inner_PCB_PortWidth = blocks_info[0]->GetDouble("Inner_PCB_PortWidth","mm");       
    Inner_PCB_StarboardWidth = blocks_info[0]->GetDouble("Inner_PCB_StarboardWidth","mm");  
    Inner_PCB_BevelAngle = blocks_info[0]->GetDouble("Inner_PCB_BevelAngle","mm");      
    Inner_PCB_UpstreamWidth = blocks_info[0]->GetDouble("Inner_PCB_UpstreamWidth","mm");   
    Inner_PCB_DownstreamWidth = blocks_info[0]->GetDouble("Inner_PCB_DownstreamWidth","mm"); 
    Inner_PCB_MidWidth = blocks_info[0]->GetDouble("Inner_PCB_MidWidth","mm");        
    Inner_PCB_Thickness = blocks_info[0]->GetDouble("Inner_PCB_Thickness","mm");       
    Inner_PCB_Ledge = blocks_info[0]->GetDouble("Inner_PCB_Ledge","mm");       
    Inner_PCB_Step = blocks_info[0]->GetDouble("Inner_PCB_Step","mm");       
    Outer_Wafer_Length = blocks_info[0]->GetDouble("Outer_Wafer_Length","mm");        
    Outer_Wafer_Width = blocks_info[0]->GetDouble("Outer_Wafer_Width","mm");         
    Outer_Wafer_Thickness = blocks_info[0]->GetDouble("Outer_Wafer_Thickness","mm");     
    Outer_Wafer_AlThickness = blocks_info[0]->GetDouble("Outer_Wafer_AlThickness","micrometer");   
    Outer_Wafer_PADExternal = blocks_info[0]->GetDouble("Outer_Wafer_PADExternal","mm");   
    Outer_Wafer_PADInternal = blocks_info[0]->GetDouble("Outer_Wafer_PADInternal","mm");   
    Outer_Wafer_GuardRing = blocks_info[0]->GetDouble("Outer_Wafer_GuardRing","mm");     
    Outer_Wafer_TransverseStrips = blocks_info[0]->GetInt("Outer_Wafer_TransverseStrips");        
    Outer_Wafer_LongitudinalStrips = blocks_info[0]->GetInt("Outer_Wafer_LongitudinalStrips");       
    Outer_PCB_PortWidth = blocks_info[0]->GetDouble("Outer_PCB_PortWidth","mm");       
    Outer_PCB_StarboardWidth = blocks_info[0]->GetDouble("Outer_PCB_StarboardWidth","mm");  
    Outer_PCB_BevelAngle = blocks_info[0]->GetDouble("Outer_PCB_BevelAngle","deg");      
    Outer_PCB_UpstreamWidth = blocks_info[0]->GetDouble("Outer_PCB_UpstreamWidth","mm");   
    Outer_PCB_DownstreamWidth = blocks_info[0]->GetDouble("Outer_PCB_DownstreamWidth","mm"); 
    Outer_PCB_MidWidth = blocks_info[0]->GetDouble("Outer_PCB_MidWidth","mm");        
    Outer_PCB_Thickness = blocks_info[0]->GetDouble("Outer_PCB_Thickness","mm");       
    Outer_PCB_Ledge = blocks_info[0]->GetDouble("Outer_PCB_Ledge","mm");       
    Outer_PCB_Step = blocks_info[0]->GetDouble("Outer_PCB_Step","mm");       
    Chamber_Thickness= blocks_info[0]->GetDouble("Chamber_Thickness","mm"); 
    Chamber_Cylinder_Length= blocks_info[0]->GetDouble("Chamber_Cylinder_Length","mm");        
    Chamber_Radius= blocks_info[0]->GetDouble("Chamber_Radius","mm");       
    Chamber_ExitTube_Radius=blocks_info[0]->GetDouble("Chamber_ExitTube_Radius","mm");
    Chamber_ExitTube_Length=blocks_info[0]->GetDouble("Chamber_ExitTube_Length","mm");
    Chamber_Flange_Inner_Radius=blocks_info[0]->GetDouble("Chamber_Flange_Inner_Radius","mm");
    Chamber_Sphere_Radius=blocks_info[0]->GetDouble("Chamber_Sphere_Radius","mm");
    Chamber_Sphere_Shift=blocks_info[0]->GetDouble("Chamber_Sphere_Shift","mm");
  }

  else{
    cout << "ERROR: check your input file formatting " << endl;
    exit(1);
  }


  // Inner Barrel
  vector<NPL::InputBlock*> blocks_inner = parser.GetAllBlocksWithTokenAndValue("Strasse","Inner");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_inner.size() << " inner detectors found " << endl; 

  vector<string> coord = {"Radius","Z","Phi","Shift","Ref"};

  for(unsigned int i = 0 ; i < blocks_inner.size() ; i++){
    if(blocks_inner[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse inner detector" << i+1 <<  endl;

      double R = blocks_inner[i]->GetDouble("Radius","mm");
      double Z= blocks_inner[i]->GetDouble("Z","mm");
      double Phi = blocks_inner[i]->GetDouble("Phi","deg");
      double Shift = blocks_inner[i]->GetDouble("Shift","mm");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_inner[i]->GetTVector3("Ref","mm"));
      AddInnerDetector(R,Z,Phi,Shift,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting on " << i+1 << " inner block " <<endl;
      exit(1);
    }
  }

  // Outer barrel
  vector<NPL::InputBlock*> blocks_outer = parser.GetAllBlocksWithTokenAndValue("Strasse","Outer");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_outer.size() << " outer detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_outer.size() ; i++){
    if(blocks_outer[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse outer detector" << i+1 <<  endl;

      double R = blocks_outer[i]->GetDouble("Radius","mm");
      double Z= blocks_outer[i]->GetDouble("Z","mm");
      double Phi = blocks_outer[i]->GetDouble("Phi","deg");
      double Shift = blocks_outer[i]->GetDouble("Shift","mm");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_inner[i]->GetTVector3("Ref","mm"));
      AddOuterDetector(R,Z,Phi,Shift,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " outer block " <<endl;
      exit(1);
    }
  }

  // Chamber
  vector<std::string> token = {"Z"};
  vector<NPL::InputBlock*> blocks_chamber = parser.GetAllBlocksWithTokenAndValue("Strasse","Chamber");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_chamber.size() << " chamber detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_chamber.size() ; i++){
    if(blocks_chamber[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse chamber detector" << i+1 <<  endl;

      double Z= blocks_chamber[i]->GetDouble("Z","mm");
      AddChamber(Z);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " chamber block " <<endl;
      exit(1);
    }
  }


  // Inactive material inside chamber imported form CAD drawings
  vector<NPL::InputBlock*> blocks_material = parser.GetAllBlocksWithTokenAndValue("Strasse","InactiveMaterial");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_material.size() << " inactive material found " << endl; 

  for(unsigned int i = 0 ; i < blocks_material.size() ; i++){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse Inactive material from CAD " << i+1 <<  endl;

      if(blocks_material[i]->HasToken("Chamber")){
          ChamberPath= blocks_material[i]->GetString("Chamber");
          found_chamber = true;
      }
      if(blocks_material[i]->HasToken("Stars")){
          StarsPath= blocks_material[i]->GetString("Stars");
          found_stars = true;
      }
      if(blocks_material[i]->HasToken("Blades")){
          BladesPath= blocks_material[i]->GetString("Blades");
          found_blades = true;
      }
      if(blocks_material[i]->HasToken("Base")){
          BasePath= blocks_material[i]->GetString("Base");
          found_base = true;
      }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Strasse::ConstructDetector(G4LogicalVolume* world){

  // Inner Barrel
  for (unsigned short i = 0 ; i < m_Inner_R.size() ; i++) {
    G4ThreeVector Det_pos = G4ThreeVector(m_Inner_Shift[i],m_Inner_R[i]+0.5*Inner_PCB_Thickness+0.001,m_Inner_Z[i]) ; //0.001 offset just to avoid overlap with frame
    Det_pos.rotate(-m_Inner_Phi[i],G4ThreeVector(0,0,1));
    G4RotationMatrix* Rot =  new G4RotationMatrix(0*deg,0*deg,m_Inner_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos+m_Inner_Ref[i]),
        BuildInnerDetector(),
        "Strasse",world,false,i+1);
  }

  // Outer Barrel 
  for (unsigned short i = 0 ; i < m_Outer_R.size() ; i++) {
    G4ThreeVector Det_pos = G4ThreeVector(m_Outer_Shift[i],m_Outer_R[i]+0.5*Inner_PCB_Thickness+0.001,m_Outer_Z[i]) ;//0.001 offset just to avoid overlap with frame
    Det_pos.rotate(-m_Outer_Phi[i],G4ThreeVector(0,0,1));
    G4RotationMatrix* Rot =  new G4RotationMatrix(0*deg,0*deg,m_Outer_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos+m_Outer_Ref[i]),
        BuildOuterDetector(),
        "Strasse",world,false,i+1);
  }

  // Chamber 
  /*
  for (unsigned short i = 0 ; i < m_Chamber_Z.size() ; i++) {
    G4ThreeVector Det_pos = G4ThreeVector(0,0,-m_Chamber_Z[i]) ;
    G4RotationMatrix* Rot =  new G4RotationMatrix();

    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
        BuildChamber(),
        "Strasse",world,false,i+1);
  }
  */


    //G4ThreeVector Det_pos = G4ThreeVector(0,0,+11.5) ;
    G4ThreeVector Det_pos = G4ThreeVector(0,0,0) ;
    G4RotationMatrix* Rot =  new G4RotationMatrix();
    Rot->rotateY(270.*deg);
    Rot->rotateX(0.*deg);

   if(found_chamber){ 
    new G4PVPlacement(Rot, Det_pos, BuildChamberFromCAD(ChamberPath),
        "Strasse_Chamber",world, false, 0);
   }

   if(found_blades){ 
    new G4PVPlacement(Rot, Det_pos, BuildBlades(BladesPath),
        "Strasse_Blades",world, false, 0);
   }
   if(found_stars){ 
    G4ThreeVector Det_pos2 = G4ThreeVector(0,0,0) ;
    new G4PVPlacement(Rot, Det_pos2, BuildStars(StarsPath),
        "Strasse_Stars",world, false, 0);
   }
   if(found_base){ 
    G4ThreeVector Det_pos3 = G4ThreeVector(0,0,0) ;
    new G4PVPlacement(Rot, Det_pos3, BuildBase(BasePath),
        "Strasse_Base",world, false, 0);
   }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Strasse::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Strasse")){
    pTree->Branch("Strasse", "TStrasseData", &m_Event) ;
  }
  pTree->SetBranchAddress("Strasse", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Strasse::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Inner barrel scorer
  DSSDScorers::PS_Rectangle* InnerScorer1= (DSSDScorers::PS_Rectangle*) m_InnerScorer1->GetPrimitive(0);

  unsigned int size = InnerScorer1->GetWidthMult(); 
  for(unsigned int i = 0 ; i < size; i++){
    double Energy = RandGauss::shoot(InnerScorer1->GetEnergyWidth(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = InnerScorer1->GetDetectorWidth(i);
      int StripTransverse = InnerScorer1->GetStripWidth(i);
      m_Event->SetInnerTE(DetNbr, StripTransverse, Energy);
    }
  }

  size = InnerScorer1->GetLengthMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    double Energy = RandGauss::shoot(InnerScorer1->GetEnergyLength(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = InnerScorer1->GetDetectorLength(i);
      int StripLongitudinal= InnerScorer1->GetStripLength(i);
      m_Event->SetInnerLE(DetNbr, StripLongitudinal, Energy);
    }
  }
  InnerScorer1->clear();

  // second silicon
  DSSDScorers::PS_Rectangle* InnerScorer2= (DSSDScorers::PS_Rectangle*) m_InnerScorer2->GetPrimitive(0);

  size = InnerScorer2->GetWidthMult(); 
  for(unsigned int i = 0 ; i < size; i++){
    double Energy = RandGauss::shoot(InnerScorer2->GetEnergyWidth(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = InnerScorer2->GetDetectorWidth(i);
      int StripTransverse = InnerScorer2->GetStripWidth(i)+Inner_Wafer_TransverseStrips;
      m_Event->SetInnerTE(DetNbr, StripTransverse, Energy);
    }
  }
  size = InnerScorer2->GetLengthMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    double Energy = RandGauss::shoot(InnerScorer2->GetEnergyLength(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = InnerScorer2->GetDetectorLength(i);
      int StripLongitudinal= InnerScorer2->GetStripLength(i);
      m_Event->SetInnerLE(DetNbr, StripLongitudinal, Energy);
    }
  }
  InnerScorer2->clear();



  ///////////
  // Outer barrel scorer
  DSSDScorers::PS_Rectangle* OuterScorer1= (DSSDScorers::PS_Rectangle*) m_OuterScorer1->GetPrimitive(0);

  size = OuterScorer1->GetWidthMult(); 
  for(unsigned int i = 0 ; i < size; i++){
    double Energy = RandGauss::shoot(OuterScorer1->GetEnergyWidth(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = OuterScorer1->GetDetectorWidth(i);
      int StripTransverse = OuterScorer1->GetStripWidth(i);
      m_Event->SetOuterTE(DetNbr, StripTransverse, Energy);
    }
  }
  size = OuterScorer1->GetLengthMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    double Energy = RandGauss::shoot(OuterScorer1->GetEnergyLength(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = OuterScorer1->GetDetectorLength(i);
      int StripLongitudinal= OuterScorer1->GetStripLength(i);
      m_Event->SetOuterLE(DetNbr, StripLongitudinal, Energy);
    }
  }
  OuterScorer1->clear();

  // Second silicon
  DSSDScorers::PS_Rectangle* OuterScorer2= (DSSDScorers::PS_Rectangle*) m_OuterScorer2->GetPrimitive(0);

  size = OuterScorer2->GetWidthMult(); 
  for(unsigned int i = 0 ; i < size; i++){
    double Energy = RandGauss::shoot(OuterScorer2->GetEnergyWidth(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = OuterScorer2->GetDetectorWidth(i);
      int StripTransverse = OuterScorer2->GetStripWidth(i)+Outer_Wafer_TransverseStrips;
      m_Event->SetOuterTE(DetNbr, StripTransverse, Energy);
    }
  }
  size = OuterScorer2->GetLengthMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    double Energy = RandGauss::shoot(OuterScorer2->GetEnergyLength(i), ResoEnergy);   
    if(Energy>EnergyThreshold){
      int DetNbr  = OuterScorer2->GetDetectorLength(i);
      int StripLongitudinal= OuterScorer2->GetStripLength(i);
      m_Event->SetOuterLE(DetNbr, StripLongitudinal, Energy);
    }
  }
  OuterScorer2->clear();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Strasse::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_InnerScorer1 = CheckScorer("InnerScorer1",already_exist) ;
  m_OuterScorer1 = CheckScorer("OuterScorer1",already_exist) ;
  m_InnerScorer2 = CheckScorer("InnerScorer2",already_exist) ;
  m_OuterScorer2 = CheckScorer("OuterScorer2",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  m_Active_InnerWafer_Width= Inner_Wafer_Width-2.*Inner_Wafer_GuardRing;
  m_Active_InnerWafer_Length= 
    Inner_Wafer_Length-Inner_Wafer_PADExternal-Inner_Wafer_PADInternal-2*Inner_Wafer_GuardRing;


  G4VPrimitiveScorer* InnerScorer1 = new DSSDScorers::PS_Rectangle("InnerScorer1",2,
      m_Active_InnerWafer_Width,
      m_Active_InnerWafer_Length,
      Inner_Wafer_LongitudinalStrips,
      Inner_Wafer_TransverseStrips,0,"xz");

  G4VPrimitiveScorer* InnerScorer2 = new DSSDScorers::PS_Rectangle("InnerScorer2",2,
      m_Active_InnerWafer_Width,
      m_Active_InnerWafer_Length,
      Inner_Wafer_LongitudinalStrips,
      Inner_Wafer_TransverseStrips,0,"xz");


  m_Active_OuterWafer_Width=Outer_Wafer_Width-2.*Outer_Wafer_GuardRing;
  m_Active_OuterWafer_Length=
    Outer_Wafer_Length-Outer_Wafer_PADExternal-Outer_Wafer_PADInternal-2*Outer_Wafer_GuardRing;


  G4VPrimitiveScorer* OuterScorer1 = new DSSDScorers::PS_Rectangle("OuterScorer1",2,
      m_Active_OuterWafer_Width,
      m_Active_OuterWafer_Length,
      Outer_Wafer_LongitudinalStrips,
      Outer_Wafer_TransverseStrips,0,"xz");

  G4VPrimitiveScorer* OuterScorer2 = new DSSDScorers::PS_Rectangle("OuterScorer2",2,
      m_Active_OuterWafer_Width,
      m_Active_OuterWafer_Length,
      Outer_Wafer_LongitudinalStrips,
      Outer_Wafer_TransverseStrips,0,"xz");



  G4VPrimitiveScorer* InteractionInner1 = new InteractionScorers::PS_Interactions("InteractionInner1",ms_InterCoord,0);
  G4VPrimitiveScorer* InteractionOuter1 = new InteractionScorers::PS_Interactions("InteractionOuter1",ms_InterCoord,0);
  G4VPrimitiveScorer* InteractionInner2 = new InteractionScorers::PS_Interactions("InteractionInner2",ms_InterCoord,0);
  G4VPrimitiveScorer* InteractionOuter2 = new InteractionScorers::PS_Interactions("InteractionOuter2",ms_InterCoord,0);


  // Register it to the multifunctionnal detector
  m_InnerScorer1->RegisterPrimitive(InnerScorer1);
  m_InnerScorer1->RegisterPrimitive(InteractionInner1);
  m_OuterScorer1->RegisterPrimitive(OuterScorer1);
  m_OuterScorer1->RegisterPrimitive(InteractionOuter1);
  m_InnerScorer2->RegisterPrimitive(InnerScorer2);
  m_InnerScorer2->RegisterPrimitive(InteractionInner2);
  m_OuterScorer2->RegisterPrimitive(OuterScorer2);
  m_OuterScorer2->RegisterPrimitive(InteractionOuter2);


  G4SDManager::GetSDMpointer()->AddNewDetector(m_InnerScorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterScorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_InnerScorer2);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterScorer2);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Strasse::Construct(){
  return  (NPS::VDetector*) new Strasse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Strasse::InitializeMaterial(){
  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialAl = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  m_MaterialCu = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Strasse{
    public:
      proxy_nps_Strasse(){
        NPS::DetectorFactory::getInstance()->AddToken("Strasse","Strasse");
        NPS::DetectorFactory::getInstance()->AddDetector("Strasse",Strasse::Construct);
      }
  };

  proxy_nps_Strasse p_nps_Strasse;
}
