/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny    contact : flavigny@lpccaen.in2p3.fr       *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Strasse Treated  data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TStrassePhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "Math/Transform3D.h"
#include "Math/RotationZYX.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
using namespace NPUNITS;
//   ROOT
#include "TChain.h"

ClassImp(TStrassePhysics)
  ///////////////////////////////////////////////////////////////////////////
  TStrassePhysics::TStrassePhysics(){
    EventMultiplicity = 0;
    m_EventData = new TStrasseData;
    m_PreTreatedData = new TStrasseData;
    m_EventPhysics = this;
    m_Spectra = NULL;
    m_E_RAW_Threshold = 0; // adc channels
    m_E_Threshold = 0;     // MeV
    m_NumberOfInnerDetectors = 0;
    m_NumberOfOuterDetectors = 0;
    m_MaximumStripMultiplicityAllowed = 10;
    m_StripEnergyMatching = 1.00;

    ////////////////////
    // Inner Detector //
    ////////////////////
    // Wafer parameter
    Inner_Wafer_Length=100*mm;
    Inner_Wafer_Width=50*mm;
    Inner_Wafer_Thickness=300*micrometer;
    Inner_Wafer_AlThickness=0.4*micrometer;
    Inner_Wafer_PADExternal=1*cm;
    Inner_Wafer_PADInternal=1*mm;
    Inner_Wafer_GuardRing=0.5*mm;

    // PCB parameter
    Inner_PCB_PortWidth=1*cm;
    Inner_PCB_StarboardWidth=2*mm;
    Inner_PCB_BevelAngle= 60*deg;
    Inner_PCB_UpstreamWidth=1*cm;
    Inner_PCB_DownstreamWidth=2*mm;
    Inner_PCB_MidWidth=1*mm;
    Inner_PCB_Thickness=3*mm;
    Inner_PCB_Ledge = 1*mm ;
    Inner_PCB_Step = 2*mm ;
    Inner_Wafer_TransverseStrips= 128;
    Inner_Wafer_LongitudinalStrips= 128;

    ////////////////////
    // Outer Detector //
    ////////////////////
    // Wafer parameter
    Outer_Wafer_Length=150*mm;
    Outer_Wafer_Width=75*mm;
    Outer_Wafer_Thickness=300*micrometer;
    Outer_Wafer_AlThickness=0.4*micrometer;
    Outer_Wafer_PADExternal=1*cm;
    Outer_Wafer_PADInternal=1*mm;
    Outer_Wafer_GuardRing=0.5*mm;

    // PCB parameter
    Outer_PCB_PortWidth=1*cm;
    Outer_PCB_StarboardWidth=2*mm;
    Outer_PCB_BevelAngle= 60*deg;
    Outer_PCB_UpstreamWidth=1*cm;
    Outer_PCB_DownstreamWidth=2*mm;
    Outer_PCB_MidWidth=1*mm;
    Outer_PCB_Thickness=3*mm;
    Outer_PCB_Ledge = 1*mm ;
    Outer_PCB_Step = 2*mm ;
    Outer_Wafer_TransverseStrips= 128;
    Outer_Wafer_LongitudinalStrips= 128;


  }

///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::AddInnerDetector(double R, double Z, double Phi, double Shift, TVector3 Ref){
  m_NumberOfInnerDetectors++;
  double ActiveWidth  = Inner_Wafer_Width-2.*Inner_Wafer_GuardRing;
  double ActiveLength = Inner_Wafer_Length-Inner_Wafer_PADExternal-Inner_Wafer_PADInternal-2*Inner_Wafer_GuardRing;
  double LongitudinalPitch = ActiveWidth/Inner_Wafer_LongitudinalStrips;
  double TransversePitch = ActiveLength/Inner_Wafer_TransverseStrips;
//cout << ActiveWidth << " " << ActiveLength << " " << LongitudinalPitch << " " << TransversePitch << endl;

  // Vector C position of detector face center
  double Recess = (Inner_PCB_Thickness-Inner_PCB_Step-Inner_Wafer_Thickness);
  TVector3 C(Shift,R+Recess,Z);// center of the whole detector, including PCB
  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips) 
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  // Adding position for downstream silicon:
  // Moving to corner of the silicon
  TVector3 P_1_1 = C
      +U*0.5*(Inner_PCB_UpstreamWidth-Inner_PCB_DownstreamWidth) // In between wafer
      -U*0.5*Inner_PCB_MidWidth // Internal wafer edge
      -U*Inner_Wafer_Length // External wafer edge
      +U*(Inner_Wafer_GuardRing+Inner_Wafer_PADExternal) // External active wafer edge
      +U*0.5*TransversePitch // middle of strip
      -V*0.5*(Inner_PCB_StarboardWidth-Inner_PCB_PortWidth)
      -V*0.5*Inner_Wafer_Width
      +V*Inner_Wafer_GuardRing
      +V*0.5*LongitudinalPitch; // middle of strip
  
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<Inner_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<Inner_Wafer_LongitudinalStrips; j++){
      P = P_1_1 + Ref + i*U*TransversePitch + j*V*LongitudinalPitch;
      lineX.push_back(P.X());
      lineY.push_back(P.Y());
      lineZ.push_back(P.Z());
      }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_InnerStripPositionX.push_back(OneDetectorStripPositionX);
  m_InnerStripPositionY.push_back(OneDetectorStripPositionY);
  m_InnerStripPositionZ.push_back(OneDetectorStripPositionZ);

  // Adding position for upstream silicon:
  // Moving to corner of the silicon
  P_1_1 = C
      +U*0.5*(Inner_PCB_UpstreamWidth-Inner_PCB_DownstreamWidth) // In between wafer
      +U*0.5*Inner_PCB_MidWidth // Internal wafer edge
      +U*(Inner_Wafer_GuardRing+Inner_Wafer_PADInternal) // Internal active wafer edge
      +U*0.5*TransversePitch// middle of strip
      -V*0.5*(Inner_PCB_StarboardWidth-Inner_PCB_PortWidth)
      -V*0.5*Inner_Wafer_Width
      +V*Inner_Wafer_GuardRing
      +V*0.5*LongitudinalPitch; // middle of strip

  for(int i=0; i<Inner_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for(int j=0; j<Inner_Wafer_LongitudinalStrips; j++){
      P = P_1_1 + Ref + i*U*TransversePitch + j*V*LongitudinalPitch;
      lineX.push_back(P.X());
      lineY.push_back(P.Y());
      lineZ.push_back(P.Z());

    }

    m_InnerStripPositionX[m_NumberOfInnerDetectors-1].push_back(lineX);
    m_InnerStripPositionY[m_NumberOfInnerDetectors-1].push_back(lineY);
    m_InnerStripPositionZ[m_NumberOfInnerDetectors-1].push_back(lineZ);

  }
}

///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::AddOuterDetector(double R, double Z, double Phi, double Shift, TVector3 Ref){
  m_NumberOfOuterDetectors++;
  double ActiveWidth  = Outer_Wafer_Width-2.*Outer_Wafer_GuardRing;
  double ActiveLength = Outer_Wafer_Length-Outer_Wafer_PADExternal-Outer_Wafer_PADInternal-2*Outer_Wafer_GuardRing;
  double LongitudinalPitch = ActiveWidth/Outer_Wafer_LongitudinalStrips;
  double TransversePitch = ActiveLength/Outer_Wafer_TransverseStrips;
//cout << ActiveWidth << " " << ActiveLength << " " << LongitudinalPitch << " " << TransversePitch << endl;

  // Vector C position of detector face center
  double Recess = (Inner_PCB_Thickness-Inner_PCB_Step-Inner_Wafer_Thickness);
  TVector3 C(Shift,R+Recess,Z);// center of the whole detector, including PCB
  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips) 
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  // Adding position for downstream silicon:
  // Moving to corner of the silicon
  TVector3 P_1_1 = C
      +U*0.5*(Outer_PCB_UpstreamWidth-Outer_PCB_DownstreamWidth) // In between wafer
      -U*0.5*Outer_PCB_MidWidth // Internal wafer edge
      -U*Outer_Wafer_Length // External wafer edge
      +U*(Outer_Wafer_GuardRing+Outer_Wafer_PADExternal) // External active wafer edge
      +U*0.5*TransversePitch // middle of strip
      -V*0.5*(Outer_PCB_StarboardWidth-Outer_PCB_PortWidth)
      -V*0.5*Outer_Wafer_Width
      +V*Outer_Wafer_GuardRing
      +V*0.5*LongitudinalPitch; // middle of strip
  
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<Outer_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<Outer_Wafer_LongitudinalStrips; j++){
      P = P_1_1 + Ref + i*U*TransversePitch + j*V*LongitudinalPitch;
      lineX.push_back(P.X());
      lineY.push_back(P.Y());
      lineZ.push_back(P.Z());
      }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_OuterStripPositionX.push_back(OneDetectorStripPositionX);
  m_OuterStripPositionY.push_back(OneDetectorStripPositionY);
  m_OuterStripPositionZ.push_back(OneDetectorStripPositionZ);

  // Adding position for upstream silicon:
  // Moving to corner of the silicon
  P_1_1 = C
      +U*0.5*(Outer_PCB_UpstreamWidth-Outer_PCB_DownstreamWidth) // In between wafer
      +U*0.5*Outer_PCB_MidWidth // Internal wafer edge
      +U*(Outer_Wafer_GuardRing+Outer_Wafer_PADInternal) // Internal active wafer edge
      +U*0.5*TransversePitch// middle of strip
      -V*0.5*(Outer_PCB_StarboardWidth-Outer_PCB_PortWidth)
      -V*0.5*Outer_Wafer_Width
      +V*Outer_Wafer_GuardRing
      +V*0.5*LongitudinalPitch; // middle of strip

  for(int i=0; i<Outer_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for(int j=0; j<Outer_Wafer_LongitudinalStrips; j++){
      P = P_1_1 + Ref + i*U*TransversePitch + j*V*LongitudinalPitch;
      lineX.push_back(P.X());
      lineY.push_back(P.Y());
      lineZ.push_back(P.Z());

    }

    m_OuterStripPositionX[m_NumberOfOuterDetectors-1].push_back(lineX);
    m_OuterStripPositionY[m_NumberOfOuterDetectors-1].push_back(lineY);
    m_OuterStripPositionZ[m_NumberOfOuterDetectors-1].push_back(lineZ);

  }
}


///////////////////////////////////////////////////////////////////////////
TVector3 TStrassePhysics::GetInnerPositionOfInteraction(const int i){
  TVector3 Position = TVector3(
      GetInnerStripPositionX(DetectorNumber[i], InnerStripT[i], InnerStripL[i]),
      GetInnerStripPositionY(DetectorNumber[i], InnerStripT[i], InnerStripL[i]),
      GetInnerStripPositionZ(DetectorNumber[i], InnerStripT[i], InnerStripL[i]));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TStrassePhysics::GetOuterPositionOfInteraction(const int i){
  TVector3 Position = TVector3(GetOuterStripPositionX(DetectorNumber[i], OuterStripT[i], OuterStripL[i]),
      GetOuterStripPositionY(DetectorNumber[i], OuterStripT[i], OuterStripL[i]),
      GetOuterStripPositionZ(DetectorNumber[i], OuterStripT[i], OuterStripL[i]));

  return Position;
}


///////////////////////////////////////////////////////////////////////////
TVector3 TStrassePhysics::GetDetectorNormal(const int i){
  //  TVector3 U = TVector3(GetStripPositionX(DetectorNumber[i],128,1),
  //      GetStripPositionY(DetectorNumber[i],128,1),
  //      GetStripPositionZ(DetectorNumber[i],128,1))
  //
  //    -TVector3(GetStripPositionX(DetectorNumber[i],1,1),
  //        GetStripPositionY(DetectorNumber[i],1,1),
  //        GetStripPositionZ(DetectorNumber[i],1,1));
  //
  //  TVector3 V = TVector3(GetStripPositionX(DetectorNumber[i],128,128),
  //      GetStripPositionY(DetectorNumber[i],128,128),
  //      GetStripPositionZ(DetectorNumber[i],128,128))
  //
  //    -TVector3(GetStripPositionX(DetectorNumber[i],128,1),
  //        GetStripPositionY(DetectorNumber[i],128,1),
  //        GetStripPositionZ(DetectorNumber[i],128,1));
  //
  //  TVector3 Normal = U.Cross(V);
  //
  //  return (Normal.Unit());
  return TVector3(0,0,0);
}
///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();
 
    if(1 /*CheckEvent() == 1*/){
      vector<TVector2> inner = MatchInner();
      vector<TVector2> outer = MatchOuter();

      for(unsigned int i=0; i<inner.size(); i++){
        int N = m_PreTreatedData->GetInner_TE_DetectorNbr(inner[i].X());
        int innerT = m_PreTreatedData->GetInner_TE_StripNbr(inner[i].X());
        int innerL = m_PreTreatedData->GetInner_LE_StripNbr(inner[i].Y());
  
        double TE = m_PreTreatedData->GetInner_TE_Energy(inner[i].X());
               // look for outer  
        double outerE = 0;
        int outerT=0;
        int outerL=0;
        for(unsigned int j=0; j<outer.size(); j++){
          if(m_PreTreatedData->GetInner_TE_DetectorNbr(outer[j].X())==N){
            outerE =  m_PreTreatedData->GetOuter_TE_Energy(outer[j].X()); 
            outerT = m_PreTreatedData->GetOuter_TE_StripNbr(outer[j].X());
            outerL = m_PreTreatedData->GetOuter_LE_StripNbr(outer[j].Y());
            }
        }

        if(outerE){
          EventMultiplicity++;
          DetectorNumber.push_back(N);
          InnerStripT.push_back(innerT);
          InnerStripL.push_back(innerL);
          DE.push_back(TE);
          InnerPosX.push_back(GetInnerPositionOfInteraction(EventMultiplicity-1).x());
          InnerPosY.push_back(GetInnerPositionOfInteraction(EventMultiplicity-1).y());
          InnerPosZ.push_back(GetInnerPositionOfInteraction(EventMultiplicity-1).z());

          OuterStripT.push_back(outerT);
          OuterStripL.push_back(outerL);
          E.push_back(outerE);
          OuterPosX.push_back(GetOuterPositionOfInteraction(EventMultiplicity-1).x());
          OuterPosY.push_back(GetOuterPositionOfInteraction(EventMultiplicity-1).y());
          OuterPosZ.push_back(GetOuterPositionOfInteraction(EventMultiplicity-1).z());
        }

      }
    }
   
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TStrassePhysics::MatchInner(){
  vector<TVector2> ArrayOfGoodCouple;

  static unsigned int m_TEMult, m_LEMult;
  m_TEMult = m_PreTreatedData->GetInnerMultTEnergy();
  m_LEMult = m_PreTreatedData->GetInnerMultLEnergy();

  if(m_TEMult>m_MaximumStripMultiplicityAllowed || m_LEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_TEMult; i++){
    for(unsigned int j=0; j<m_LEMult; j++){

      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetInner_TE_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetInner_LE_DetectorNbr(j);

      // if same detector check energy
      if(XDetNbr == YDetNbr){
        // Declaration of variable for clarity
        double TE = m_PreTreatedData->GetInner_TE_Energy(i);
        double LE = m_PreTreatedData->GetInner_LE_Energy(j);

        // look if energy matches
        if(abs(TE-LE)/2.<m_StripEnergyMatching){
          ArrayOfGoodCouple.push_back(TVector2(i,j));
        }
      }
    }
  }

  return ArrayOfGoodCouple;
}
///////////////////////////////////////////////////////////////////////////
vector<TVector2> TStrassePhysics::MatchOuter(){
  vector<TVector2> ArrayOfGoodCouple;

  static unsigned int m_TEMult, m_LEMult;
  m_TEMult = m_PreTreatedData->GetOuterMultTEnergy();
  m_LEMult = m_PreTreatedData->GetOuterMultLEnergy();

  if(m_TEMult>m_MaximumStripMultiplicityAllowed || m_LEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_TEMult; i++){
    for(unsigned int j=0; j<m_LEMult; j++){

      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetOuter_TE_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetOuter_LE_DetectorNbr(j);

      // if same detector check energy
      if(XDetNbr == YDetNbr){
        // Declaration of variable for clarity
        double TE = m_PreTreatedData->GetOuter_TE_Energy(i);
        double LE = m_PreTreatedData->GetOuter_LE_Energy(j);

        // look if energy matches
        if(abs(TE-LE)/2.<m_StripEnergyMatching){
          ArrayOfGoodCouple.push_back(TVector2(i,j));
        }
      }
    }
  }

  return ArrayOfGoodCouple;
}


///////////////////////////////////////////////////////////////////////////
int TStrassePhysics::CheckEvent(){
  // Check the size of the different elements
  if(m_PreTreatedData->GetInnerMultTEnergy() == m_PreTreatedData->GetInnerMultLEnergy() )
    return 1;

  else
    return -1;
}


///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  //////
  // First Stage Energy
  unsigned int sizeFront = m_EventData->GetInnerMultTEnergy();
  for (UShort_t i = 0; i < sizeFront ; ++i) {
    if (m_EventData->GetInner_TE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetInner_TE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetInner_TE_DetectorNbr(i)),m_EventData->GetInner_TE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetInnerTE(m_EventData->GetInner_TE_DetectorNbr(i), m_EventData->GetInner_TE_StripNbr(i), Energy);
      }
    }
  }
  unsigned int sizeBack = m_EventData->GetInnerMultTEnergy();
  for (UShort_t i = 0; i < sizeBack ; ++i) {
    if (m_EventData->GetInner_LE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetInner_LE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetInner_LE_DetectorNbr(i)),m_EventData->GetInner_LE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetInnerLE(m_EventData->GetInner_LE_DetectorNbr(i), m_EventData->GetInner_LE_StripNbr(i), Energy);
      }
    }
  }

  //////
  // Second Stage Energy
  sizeFront = m_EventData->GetOuterMultTEnergy();
  for (UShort_t i = 0; i < sizeFront ; ++i) {
    if (m_EventData->GetOuter_TE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetOuter_TE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetOuter_TE_DetectorNbr(i)),m_EventData->GetOuter_TE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetOuterTE(m_EventData->GetOuter_TE_DetectorNbr(i), m_EventData->GetOuter_TE_StripNbr(i), Energy);
      }
    }
  }
  sizeBack = m_EventData->GetOuterMultTEnergy();
  for (UShort_t i = 0; i < sizeBack ; ++i) {
    if (m_EventData->GetOuter_LE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetOuter_LE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("Strasse/ENERGY"+NPL::itoa(m_EventData->GetOuter_LE_DetectorNbr(i)),m_EventData->GetOuter_LE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetOuterLE(m_EventData->GetOuter_LE_DetectorNbr(i), m_EventData->GetOuter_LE_StripNbr(i), Energy);
      }
    }
  }

}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigStrasse.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigStrasse.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigStrasse.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigStrasse.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigStrasse";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::Clear() {
  EventMultiplicity = 0;
  // DSSD
  DetectorNumber.clear();
  E.clear();
  InnerStripT.clear();
  InnerStripL.clear();
  OuterStripT.clear();
  OuterStripL.clear();
  DE.clear();

  // Position Information
  InnerPosX.clear();
  InnerPosY.clear();
  InnerPosZ.clear();
  OuterPosX.clear();
  OuterPosY.clear();
  OuterPosZ.clear();


}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::ReadConfiguration(NPL::InputParser parser) {
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
  }

  else{
    cout << "ERROR: check your input file formatting " << endl;
    exit(1);
  }


  // Inner Barrel
  vector<NPL::InputBlock*> blocks_inner = parser.GetAllBlocksWithTokenAndValue("Strasse","Inner");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_inner.size() << " inner detectors found " << endl; 

  vector<string> coord = {"Radius","Z","Phi","Shift"};

  for(unsigned int i = 0 ; i < blocks_inner.size() ; i++){
    if(blocks_inner[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Strasse inner detector" << i+1 <<  endl;

      double R = blocks_inner[i]->GetDouble("Radius","mm");
      double Z= blocks_inner[i]->GetDouble("Z","mm");
      double Phi = blocks_inner[i]->GetDouble("Phi","deg");
      double Shift = blocks_inner[i]->GetDouble("Shift","mm");
      TVector3 Ref = blocks_inner[i]->GetTVector3("Ref","mm");
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
      TVector3 Ref = blocks_outer[i]->GetTVector3("Ref","mm");
      AddOuterDetector(R,Z,Phi,Shift,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " outer block " <<endl;
      exit(1);
    }
  }

}

///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitSpectra() {
  //  m_Spectra = new TStrasseSpectra(m_NumberOfInnerDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TStrassePhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  /*  for (int i = 0; i < m_NumberOfInnerDetectors; ++i) {
      Cal->AddParameter("Strasse", "INNER"+ NPL::itoa(i+1)+"_ENERGY","Strasse_INNER"+ NPL::itoa(i+1)+"_ENERGY");
      }
      for (int i = 0; i < m_NumberOfInnerDetectors; ++i) {
      Cal->AddParameter("Strasse", "OUTER"+ NPL::itoa(i+1)+"_ENERGY","Strasse_OUTER"+ NPL::itoa(i+1)+"_ENERGY");
      }
      */
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Strasse",  true );
  inputChain->SetBranchAddress("Strasse", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Strasse", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TStrassePhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Strasse", "TStrassePhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TStrassePhysics::Construct() {
  return (NPL::VDetector*) new TStrassePhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_Strasse{
    public:
      proxy_Strasse(){
        NPL::DetectorFactory::getInstance()->AddToken("Strasse","Strasse");
        NPL::DetectorFactory::getInstance()->AddDetector("Strasse",TStrassePhysics::Construct);
      }
  };

  proxy_Strasse p_Strasse;
}

