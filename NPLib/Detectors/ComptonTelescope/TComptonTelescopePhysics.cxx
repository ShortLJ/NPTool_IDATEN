/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ComptonTelescope treated data                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TComptonTelescopePhysics.h"
using namespace ComptonTelescope_LOCAL;

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <numeric>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "NPOptionManager.h"
#include "NPDetectorFactory.h"

//   ROOT
#include "TChain.h"


ClassImp(TComptonTelescopePhysics)


  ///////////////////////////////////////////////////////////////////////////
TComptonTelescopePhysics::TComptonTelescopePhysics()
  : m_EventData(new TComptonTelescopeData),
  m_PreTreatedData(new TComptonTelescopeData),
  m_EventPhysics(this),
  m_Spectra(0),
  m_nCounterEvt(50),
  m_nCounterHit(50),
  m_MaximumStripMultiplicityAllowed(32),
  m_MultOneOnly(false),
  m_StripEnergyMatchingSigma(0.006),      // MeV
  m_StripEnergyMatchingNumberOfSigma(3),  
  m_StripFront_E_RAW_Threshold(0),
  m_StripFront_E_Threshold(0),  // MeV
  m_StripBack_E_RAW_Threshold(0),
  m_StripBack_E_Threshold(0),  // MeV
  m_Take_E_Front(true), // p-side
  m_NumberOfDetectors(0),
  m_NumberOfStrips(32),
  m_NPixels(64)
{
  EventMultiplicity   = 0;
  for (Int_t i = 0; i < m_nCounterEvt; i++) m_CounterEvt[i] = 0;
  for (Int_t i = 0; i < m_nCounterHit; i++) m_CounterHit[i] = 0;
}

///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::BuildPhysicalEvent()
{
  BuildSimplePhysicalEvent();
}


///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::BuildSimplePhysicalEvent()
{
  //cout << "\nBegin event treatment" << endl;
  //m_CounterEvt[0] = 1; // total nb of events

  // select active channels and apply threhsolds
  PreTreat();

  //// DSSSD analysis ////

  // Check event type
  Int_t evtType = CheckEvent();
  //cout << "event type = " << evtType << endl;
 
  // check possible interstrip
/*  if (evtType == 2) {
    if (m_PreTreatedData->GetCTTrackerFrontEMult() == 2 && m_PreTreatedData->GetCTTrackerBackEMult() == 1 || m_PreTreatedData->GetCTTrackerFrontEMult() == 1 && m_PreTreatedData->GetCTTrackerBackEMult() == 2) {
      //cout << "event type 2 possible interstrip" << endl;
      m_CounterEvt[21] = 1; // possible interstrip mult 2-1
      if (m_PreTreatedData->GetCTTrackerFrontEMult() == 2) {
        if (m_PreTreatedData->GetCTTrackerFrontEStripNbr(0) == m_PreTreatedData->GetCTTrackerFrontEStripNbr(1)+1 || m_PreTreatedData->GetCTTrackerFrontEStripNbr(0) == m_PreTreatedData->GetCTTrackerFrontEStripNbr(1)-1) {
          m_CounterEvt[22] = 1; // possible interstrip close SF         
          //cout << "SFi = " << m_PreTreatedData->GetCTTrackerFrontEStripNbr(0) << " SFj = " << m_PreTreatedData->GetCTTrackerFrontEStripNbr(1) << endl;
          //cout << "EB " << m_PreTreatedData->GetCTTrackerBackEEnergy(0) << " EFi " << m_PreTreatedData->GetCTTrackerFrontEEnergy(0) << " EFj " << m_PreTreatedData->GetCTTrackerFrontEEnergy(1) << endl;
          if (abs((m_PreTreatedData->GetCTTrackerBackEEnergy(0)-(m_PreTreatedData->GetCTTrackerFrontEEnergy(0)+m_PreTreatedData->GetCTTrackerFrontEEnergy(1)))/2.) < m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma) {
            m_CounterEvt[24] = 1; // interstrip close SF
            //cout << "interstrip" << endl;
          }
        }
      }
      else if (m_PreTreatedData->GetCTTrackerBackEMult() == 2) {
        if (m_PreTreatedData->GetCTTrackerBackEStripNbr(0) == m_PreTreatedData->GetCTTrackerBackEStripNbr(1)+1 || m_PreTreatedData->GetCTTrackerBackEStripNbr(0) == m_PreTreatedData->GetCTTrackerBackEStripNbr(1)-1) {
          m_CounterEvt[23] = 1; // possible interstrip close SB
          //cout << "SBi = " << m_PreTreatedData->GetCTTrackerBackEStripNbr(0) << " SBj = " << m_PreTreatedData->GetCTTrackerBackEStripNbr(1) << endl;
          if (abs((m_PreTreatedData->GetCTTrackerFrontEEnergy(0)-(m_PreTreatedData->GetCTTrackerBackEEnergy(0)+m_PreTreatedData->GetCTTrackerBackEEnergy(1)))/2.) < m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma) {
            m_CounterEvt[25] = 1; // interstrip close SB
          }
        }
      }
    }
  }
*/

  // Remove: do general case
  //if (CheckEvent() == 1) {   // case where multiplicity front = multiplicity back
  
  vector<TVector2> couple = Match_Front_Back();
  EventMultiplicity = couple.size();
 //cout << "event multiplicity = " << couple.size() << endl;

  // keep only mult 1 couples
  //if (couple.size() ==  1) { // pb if done here, so done in Match_Front_Back()
  for (UShort_t i = 0; i < couple.size(); ++i) { // loop on selected events
    //m_CounterEvt[16] = 1; // nb of physics events
    //m_CounterHit[10] += 1; // nb of physics hits

    // Energy
    Int_t Tower = m_PreTreatedData->GetCTTrackerFrontETowerNbr(couple[i].X());
    Int_t    N       = m_PreTreatedData->GetCTTrackerFrontEDetectorNbr(couple[i].X());
    Int_t    Front   = m_PreTreatedData->GetCTTrackerFrontEStripNbr(couple[i].X());
    Int_t    Back    = m_PreTreatedData->GetCTTrackerBackEStripNbr(couple[i].Y());
    Double_t Front_E = m_PreTreatedData->GetCTTrackerFrontEEnergy(couple[i].X());
    Double_t Back_E  = m_PreTreatedData->GetCTTrackerBackEEnergy(couple[i].Y());

    //cout << "mult " << couple.size() << " event type " << evtType << endl;
    // Event type
/*    if (evtType == 1) {
      m_CounterEvt[17] = 1; // nb of physics events with mult F = mult B
      //cout << "event type 1: couple size " << couple.size() << " couple " <<i << " SF" << Front << " SB" << Back << " EF " << Front_E << " EB " << Back_E << endl;
    }
    if (evtType == 2) {
      if (couple.size() == 1) {
        m_CounterEvt[18] = 1; // nb of physics events with mult F = mult B +- 1
        //cout << "event type 2: couple size " << couple.size() << " couple " << i << " SF" << Front << " SB" << Back << " EF " << Front_E << " EB " << Back_E << endl;
      }
    }
    if (evtType == -1) m_CounterEvt[19] = 1; // nb of physics events with mult F != mult B or mult B +- 1
*/

    // Time
    // Front
    //cout << "time front multiplicity = " << m_PreTreatedData->GetCTTrackerFrontTMult() << endl;
    Double_t Front_T = -1000.;
    for (UShort_t t = 0; t < m_PreTreatedData->GetCTTrackerFrontTMult(); t++) {
      if (m_PreTreatedData->GetCTTrackerFrontETowerNbr(couple[i].X()) == m_PreTreatedData->GetCTTrackerFrontTTowerNbr(t) && m_PreTreatedData->GetCTTrackerFrontEDetectorNbr(couple[i].X()) == m_PreTreatedData->GetCTTrackerFrontTDetectorNbr(t)) {
        //m_CounterEvt[42] = 1; // nb of physics events with front time
        //m_CounterHit[42] += 1; // nb of physics hits with front time
        Front_T = m_PreTreatedData->GetCTTrackerFrontTTime(t);
        //cout << "Time F: Tower" << m_PreTreatedData->GetCTTrackerFrontTTowerNbr(t) << " D" << m_PreTreatedData->GetCTTrackerFrontTDetectorNbr(t) << " T = " << Front_T << endl;
      }
    }
    // Back
    //cout << "time back  multiplicity = " << m_PreTreatedData->GetCTTrackerBackTMult() << endl;
    Double_t Back_T = -1000;
    for (UShort_t t = 0; t < m_PreTreatedData->GetCTTrackerBackTMult(); t++) {
      if (m_PreTreatedData->GetCTTrackerBackETowerNbr(couple[i].Y()) == m_PreTreatedData->GetCTTrackerBackTTowerNbr(t) && m_PreTreatedData->GetCTTrackerBackEDetectorNbr(couple[i].Y()) == m_PreTreatedData->GetCTTrackerBackTDetectorNbr(t)) {
        //m_CounterEvt[43] = 1; // nb of physics events with back time
        //m_CounterHit[43] += 1; // nb of physics hits with back time
        Back_T = m_PreTreatedData->GetCTTrackerBackTTime(t);
        //cout << "Time B: Tower" << m_PreTreatedData->GetCTTrackerBackTTowerNbr(t) << " D" << m_PreTreatedData->GetCTTrackerBackTDetectorNbr(t) << " T = " << Back_T << endl;
      }
    }
    // check time
/*    bool Same_T = false;
    if (Front_T == Back_T) {
      //cout << "same T" << endl;
      Same_T = true;
      m_CounterEvt[44] = 1; // nb of physics events with same FB time
      m_CounterHit[44] += 1; // nb of physics hits with same FB time
    }
*/
    //cout << "couple " << i << " CT" << Tower << " D" << N << " SF" << Front << " SB" << Back << " EF " << Front_E << " EB " << Back_E << " TF " << Front_T << " TB " << Back_T << endl;

    // Fill TComptonTelescopePhysics members
    EventType.push_back(evtType);
    TowerNumber.push_back(Tower);
    DetectorNumber.push_back(N);
    Strip_Front.push_back(Front);
    Strip_Back.push_back(Back);
    Front_Energy.push_back(Front_E);
    Back_Energy.push_back(Back_E);
    Half_Energy.push_back((Front_E+Back_E)/2);
    Front_Time.push_back(Front_T);
    Back_Time.push_back(Back_T);
    //Same_FBTime.push_back(Same_T);

/*    if (m_Take_E_Front)
      Strip_E.push_back(Front_E);
    else
      Strip_E.push_back(Back_E);
*/
    //}
  }
//  } // end check event


  //// Calorimeter analysis ////
  int nCalorTriggered = m_PreTreatedData -> GetCTCalorimeterTMult();

/*  double charge = 0;
  unsigned int maxIndex = 0;
  int max = 0;*/
  double sumADC = 0;
  double E_calib = 0;
  unsigned int cursor = 0;
  UShort_t detectorNumber = 0;
  for (int j = 0; j < nCalorTriggered; j++) {
    cursor = j*m_NPixels;
    // Calculate an approximate position of interaction
/*    maxIndex = 0;
    max = m_PreTreatedData->GetCTCalorimeterEEnergy(maxIndex);
    for (unsigned int i = 1; i < m_NPixels; i++) {
      if (max < m_PreTreatedData->GetCTCalorimeterEEnergy(cursor+i)) {
        maxIndex = i;
      }
    }//  int maxIndex = max_element(mat.begin(), mat.end()) - mat.begin(); cout << maxIndex << "x, y: " << 6*(maxIndex/8)-21 << ", " << 6*(maxIndex%8)-21 << endl;
    CalorPosX.push_back(6*(maxIndex/8)-21);
    CalorPosY.push_back(6*(maxIndex%8)-21);*/
  
    // Export pretreated data for NN analysis
    detectorNumber = m_PreTreatedData->GetCTCalorimeterEDetectorNbr(j*m_NPixels);
    vector<int> data;
    data.clear();
    data.resize(m_NPixels, 0);
    vector<int> data_align; // for online analysis
    data_align.clear();
    data_align.resize(m_NPixels, 0);
    for (int i = 0; i < m_NPixels; i++) {
      if (m_PreTreatedData->GetCTCalorimeterEChannelNbr(i+cursor) == i) {
        data[i] = m_PreTreatedData->GetCTCalorimeterEEnergy(i+cursor);
        data_align[i] = fCalorimeter_coeffAlign(m_PreTreatedData, i+cursor);
        //cout << "data pretreat " << data[i] << " ; align " << data_align[i] << endl;
      } else {
        cout << "Inconsistency found in channel ordering while filling CalorData: field #" << cursor+i << endl;
      }
    }
    CalorData.insert(pair<int, vector<int>>(detectorNumber, data));
/*    CalorID.push_back(detectorNumber);
    for (int i = cursor; i < cursor+m_NPixels; i++) {
      CalorData.push_back(m_PreTreatedData->GetCTCalorimeterEEnergy(i));
    }*/
  
    // Calculate a corrected energy for the calorimeter
    /*for (UShort_t i = 0; i < m_PreTreatedData->GetCTCalorimeterEMult(); ++i) {
      charge += fCalorimeter_Q(m_PreTreatedData, i);//Apply calibration other than pedestal and sum anodes
    }*/
/*    for (UShort_t i = cursor; i < cursor+m_NPixels; ++i) {
      charge += fCalorimeter_Q(m_EventData, i);//Apply full calibration and sum anodes
    }
    Calor_E.push_back(fCalorimeter_E(charge, detectorNumber));*/

    Calor_E.push_back(accumulate(data.begin(), data.end(), 0));// Uncorrected uncalibrated energy
  
    Calor_T.push_back(m_PreTreatedData->GetCTCalorimeterTTime(j));

    // for online analysis: calibrated energy
    sumADC = accumulate(data_align.begin(), data_align.end(), 0);
    //cout << "sum ADC " << accumulate(data.begin(), data.end(), 0) << " align " << sumADC << endl;

    E_calib = 0;
    E_calib = fCalorimeter_calibE(sumADC, detectorNumber);
    Calor_E_calib.push_back(E_calib);
    //cout << "det " << detectorNumber << " E " << E_calib << endl;
    

  }//End of loop on triggered calorimeter detectors

  /// Delta T analysis ///
  if (EventMultiplicity+nCalorTriggered > 1) {
    for (int i = 0; i < EventMultiplicity; i++) {
      for (int j = 0; j < nCalorTriggered; j++) {
        deltaT.push_back(Front_Time[i]-Calor_T[j]);
      }
    }
  }
  resetCount = m_EventData -> GetResetCount();

  //   if (DetectorNumber.size() == 1) return;
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::PreTreat()
{
  // Clear pre treated object
  ClearPreTreatedData();

  // Front, energy
  //cout << m_EventData->GetCTTrackerFrontEMult() << endl;
  //m_EventData->Dump();
  for (UShort_t i = 0; i < m_EventData->GetCTTrackerFrontEMult(); ++i) {
    //m_CounterEvt[1] = 1; // nb of events with at least one EF raw recorded
    //m_CounterHit[0] += 1; // nb of hits with EF raw
    //cout << "Raw: DetF = " << m_EventData->GetCTTrackerFrontEDetectorNbr(i) << " ; stripF = " << m_EventData->GetCTTrackerFrontEStripNbr(i) << " ; EF raw = " << m_EventData->GetCTTrackerFrontEEnergy(i) << endl;

    if (m_EventData->GetCTTrackerFrontEEnergy(i) > m_StripFront_E_RAW_Threshold &&
        IsValidChannel("Front", m_EventData->GetCTTrackerFrontEDetectorNbr(i), m_EventData->GetCTTrackerFrontEStripNbr(i))) {     
      //m_CounterEvt[2] = 1; // nb of events with at least one EF raw > threshold
      //m_CounterHit[1] += 1; // nb of hits with EF raw > threshold
      Double_t E = fStrip_Front_E(m_EventData, i)*1e3;//Calibration happens here, E in keV
      
      if (E > m_StripFront_E_Threshold) {
        //cout << "CalF: Det = " << m_EventData->GetCTTrackerFrontEDetectorNbr(i) << " ; strip = " << m_EventData->GetCTTrackerFrontEStripNbr(i) << " ; E cal = " << E << endl;
        //m_CounterEvt[3] = 1; // nb of events with at least one EF cal > threshold
        //m_CounterHit[2] += 1; // nb of hits with EF cal > threshold
        m_PreTreatedData->SetFrontE(
            m_EventData->GetCTTrackerFrontETowerNbr(i),
            m_EventData->GetCTTrackerFrontEDetectorNbr(i),
            m_EventData->GetCTTrackerFrontEStripNbr(i),
            E);
      }
    }
  }

  // Back, energy
  for (UShort_t i = 0; i < m_EventData->GetCTTrackerBackEMult(); ++i) {
    //m_CounterEvt[4] = 1; // nb of events with at least one EB raw recorded
    //m_CounterHit[3] += 1; // nb of hits with EB raw
    //cout << "Raw: DetB = " << m_EventData->GetCTTrackerBackEDetectorNbr(i) << " ; stripB = " << m_EventData->GetCTTrackerBackEStripNbr(i) << " ; EB raw = " << m_EventData->GetCTTrackerBackEEnergy(i) << endl;

    if (m_EventData->GetCTTrackerBackEEnergy(i) > m_StripBack_E_RAW_Threshold && 
        IsValidChannel("Back", m_EventData->GetCTTrackerBackEDetectorNbr(i), m_EventData->GetCTTrackerBackEStripNbr(i))) {
      //m_CounterEvt[5] = 1; // nb of events with at least one EB raw > threshold
      //m_CounterHit[4] += 1; // nb of hits with EB raw > threshold
      Double_t E = fStrip_Back_E(m_EventData, i)*1e3;//Calibration happens here, E in keV

      if (E > m_StripBack_E_Threshold) {
        //cout << "CalB: Det = " << m_EventData->GetCTTrackerBackEDetectorNbr(i) << " ; stripB = " << m_EventData->GetCTTrackerBackEStripNbr(i) << " ; EB cal = " << E << endl;
        //m_CounterEvt[6] = 1; // nb of events with at least one EB cal > threshold
        //m_CounterHit[5] += 1; // nb of hits with EB cal > threshold
        m_PreTreatedData->SetBackE(
            m_EventData->GetCTTrackerBackETowerNbr(i),
            m_EventData->GetCTTrackerBackEDetectorNbr(i),
            m_EventData->GetCTTrackerBackEStripNbr(i),
            E);
      }
    }
  }

  // DSSSD time information and calorimeter still have to be done...
  // Front, time
  for (UShort_t i = 0; i < m_EventData->GetCTTrackerFrontTMult(); ++i) {
    //m_CounterEvt[40] = 1; // nb of events with FT
    //m_CounterHit[40] += 1; // nb of hits with FT
    m_PreTreatedData->SetFrontT(
        m_EventData->GetCTTrackerFrontTTowerNbr(i),
        m_EventData->GetCTTrackerFrontTDetectorNbr(i),
        m_EventData->GetCTTrackerFrontTStripNbr(i),
        m_EventData->GetCTTrackerFrontTTime(i));
    //cout << "Pretreat time front : T" << m_EventData->GetCTTrackerFrontTTowerNbr(i) << " D" << m_EventData->GetCTTrackerFrontTDetectorNbr(i) << " Strip " << m_EventData->GetCTTrackerFrontTStripNbr(i) << " time " << m_EventData->GetCTTrackerFrontTTime(i) << endl;
  }

  // Back, time
  for (UShort_t i = 0; i < m_EventData->GetCTTrackerBackTMult(); ++i) {
    //m_CounterEvt[41] = 1; // nb of events with BT
    //m_CounterHit[41] += 1; // nb of hits with BT
    m_PreTreatedData->SetBackT(
        m_EventData->GetCTTrackerBackTTowerNbr(i),
        m_EventData->GetCTTrackerBackTDetectorNbr(i),
        m_EventData->GetCTTrackerBackTStripNbr(i),
        m_EventData->GetCTTrackerBackTTime(i));
    //cout << "Pretreat time back : T" << m_EventData->GetCTTrackerBackTTowerNbr(i) << " D" << m_EventData->GetCTTrackerBackTDetectorNbr(i) << " Strip " << m_EventData->GetCTTrackerBackTStripNbr(i) << " time " << m_EventData->GetCTTrackerBackTTime(i) << endl;
  }


  // Calorimeter
  int data[m_NPixels];
  //cout << "mult calo " << m_EventData->GetCTCalorimeterTMult() << endl;
  for (UShort_t i = 0; i < m_EventData->GetCTCalorimeterTMult(); ++i) {
    for (int j = 0; j < m_NPixels; j++) {
      //cout << "data " << m_EventData->GetCTCalorimeterEEnergy(j) << endl;
      data[j] = fCalorimeter_ped(m_EventData, i*m_NPixels+j);
      //cout << "data - ped " << data[j] << endl;
    }
    m_PreTreatedData -> SetCTCalorimeter(m_EventData->GetCTCalorimeterTTowerNbr(i), m_EventData->GetCTCalorimeterEDetectorNbr(i*m_NPixels), m_EventData->GetCTCalorimeterTChannelNbr(i), m_EventData->GetCTCalorimeterTTime(i), data, m_NPixels);
  }

  m_PreTreatedData -> SetResetCount(m_EventData -> GetResetCount());
}



///////////////////////////////////////////////////////////////////////////
int TComptonTelescopePhysics::CheckEvent()
{

  // same multiplicity on front and back side 
  if (m_PreTreatedData->GetCTTrackerBackEMult() == m_PreTreatedData->GetCTTrackerFrontEMult()) {
    //cout << "mult event type 1 = " << m_PreTreatedData->GetCTTrackerBackEMult() << endl;
    //m_CounterEvt[8] = 1; // nb of pretreated events with mult F = mult B
    return 1 ; // Regular Event
  }

  // possibly interstrip
  else if (m_PreTreatedData->GetCTTrackerFrontEMult() == m_PreTreatedData->GetCTTrackerBackEMult()+1 || 
      m_PreTreatedData->GetCTTrackerFrontEMult() == m_PreTreatedData->GetCTTrackerBackEMult()-1) {
    //m_CounterEvt[9] = 1; // nb of pretreated events with mult F = mult B +- 1
    return 2;
  }

  else {
    //m_CounterEvt[10] = 1; // nb of pretreated events with mult F != mult B or mult B +- 1
    return -1 ; // Rejected Event
  }
}



///////////////////////////////////////////////////////////////////////////
vector<TVector2> TComptonTelescopePhysics::Match_Front_Back()
{
  vector<TVector2> ArrayOfGoodCouple;

  // Select allowed multiplicity events. If multiplicity is too high, then return "empty" vector
/*  if (m_PreTreatedData->GetCTTrackerFrontEMult() > m_MaximumStripMultiplicityAllowed || 
      m_PreTreatedData->GetCTTrackerBackEMult() > m_MaximumStripMultiplicityAllowed)
    return ArrayOfGoodCouple;
*/

  // Loop on front multiplicity
  for (UShort_t i = 0; i < m_PreTreatedData->GetCTTrackerFrontEMult(); i++) {
    // Loop on back multiplicity
    for (UShort_t j = 0; j < m_PreTreatedData->GetCTTrackerBackEMult(); j++) {
      // if same tower and same detector check energy
      if ((m_PreTreatedData->GetCTTrackerFrontETowerNbr(i) == m_PreTreatedData->GetCTTrackerBackETowerNbr(j)) && 
          (m_PreTreatedData->GetCTTrackerFrontEDetectorNbr(i) == m_PreTreatedData->GetCTTrackerBackEDetectorNbr(j))) {
        //m_CounterEvt[12] = 1; // nb of pretreated events with EF and EB in same tower and detector
        //m_CounterHit[7] += 1; // nb of hits with EF and EB in same tower and detector
        
        // equal energy
        if (abs((m_PreTreatedData->GetCTTrackerFrontEEnergy(i) - m_PreTreatedData->GetCTTrackerBackEEnergy(j))/2.) < m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma) {
          //m_CounterEvt[13] = 1; // nb of pretreated events with same D and E
          //m_CounterHit[8] += 1; // nb of hits with same D and E
          ArrayOfGoodCouple.push_back(TVector2(i,j));
        } // end test energy

        // add interstrip
        // only for mult 2 - front case
/*        if (m_PreTreatedData->GetCTTrackerFrontEMult() == 2 && m_PreTreatedData->GetCTTrackerBackEMult() == 1) {
          // if close strip
          if (m_PreTreatedData->GetCTTrackerFrontEStripNbr(0) == m_PreTreatedData->GetCTTrackerFrontEStripNbr(1)+1 || m_PreTreatedData->GetCTTrackerFrontEStripNbr(0) == m_PreTreatedData->GetCTTrackerFrontEStripNbr(1)-1) {
            // if same energy
            if (abs((m_PreTreatedData->GetCTTrackerBackEEnergy(0)-(m_PreTreatedData->GetCTTrackerFrontEEnergy(0)+m_PreTreatedData->GetCTTrackerFrontEEnergy(1)))/2.) < m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma) {
              m_CounterEvt[26] = 1; // interstrip close SF
              ArrayOfGoodCouple.push_back(TVector2(i,j));
            }
          }
        }*/
      } // end test same tower and detector
    } // end loop back multiplicity
  } // end loop front multiplicity
 
  //cout << "ArrayOfGoodCouple initial size = " << ArrayOfGoodCouple.size() << endl;

  // prevent treating event with ambiguous matching beetween X and Y
  // not done
  if (ArrayOfGoodCouple.size() > m_PreTreatedData->GetCTTrackerFrontEMult()) {
    //m_CounterEvt[14] = 1; // nb of pretreated events with ambiguous matching
    //ArrayOfGoodCouple.clear();
  }

  // keep only mult = 1
  if (m_MultOneOnly) {
    if (ArrayOfGoodCouple.size() > 1) {
      //m_CounterEvt[15] = 1; // nb of events with couple size > 1
      ArrayOfGoodCouple.clear();
    }
  }

  return ArrayOfGoodCouple;
}



////////////////////////////////////////////////////////////////////////////
bool TComptonTelescopePhysics::IsValidChannel(const string DetectorType, const int detector, const int channel)
{
  if (DetectorType == "Front")
    return *(m_FrontChannelStatus[detector-1].begin()+channel);

  else if (DetectorType == "Back")
    return *(m_BackChannelStatus[detector-1].begin()+channel);

  else return false;
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::ReadAnalysisConfig()
{
  bool ReadingStatus = false;

  cout << "\t/////////// Reading ConfigComptonTelescope.dat file ///////////" << endl;

  // path to file
  string FileName = "./configs/ConfigComptonTelescope.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << "\tNo ConfigComptonTelescope.dat found: default parameters loaded for Analysis " << FileName << endl;
    return;
  }
  cout << "\tLoading user parameters from ConfigComptonTelescope.dat " << endl;

  // storing config file in the ROOT output file
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigComptonTelescope.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");

  // read analysis config file
  string LineBuffer, DataBuffer, whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    if (LineBuffer.compare(0, 22, "ConfigComptonTelescope") == 0) ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {
      whatToDo = "";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo == "MAX_STRIP_MULTIPLICITY") {
        AnalysisConfigFile >> DataBuffer;
        m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str());
        cout << "\t" << whatToDo << "\t" << m_MaximumStripMultiplicityAllowed << endl;
      }
      
      else if (whatToDo == "ONLY_GOOD_COUPLE_WITH_MULTIPLICITY_ONE") {
        AnalysisConfigFile >> DataBuffer;
        if (DataBuffer == "ON") m_MultOneOnly = true;
        cout << "\t" << whatToDo << "\t" << DataBuffer << endl;
      }

      else if (whatToDo == "FRONT_BACK_ENERGY_MATCHING_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingSigma = atof(DataBuffer.c_str());
        cout << "\t" << whatToDo << "\t" << m_StripEnergyMatchingSigma << endl;
      }

      else if (whatToDo == "FRONT_BACK_ENERGY_MATCHING_NUMBER_OF_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingNumberOfSigma = atoi(DataBuffer.c_str());
        cout << "\t" << whatToDo << "\t" << m_StripEnergyMatchingNumberOfSigma << endl;
      }

      else if (whatToDo == "DISABLE_ALL") {
        AnalysisConfigFile >> DataBuffer;
        cout << "\t" << whatToDo << "\t" << DataBuffer << endl;
        Int_t Detector = atoi(DataBuffer.substr(18,1).c_str());
        vector< bool > ChannelStatus;
        ChannelStatus.resize(m_NumberOfStrips, false);
        m_FrontChannelStatus[Detector-1] = ChannelStatus;
        m_BackChannelStatus[Detector-1]  = ChannelStatus;
      }

      else if (whatToDo == "DISABLE_CHANNEL") {
        AnalysisConfigFile >> DataBuffer;
        cout << "\t" << whatToDo << "\t" << DataBuffer << endl;
        //Int_t Detector = atoi(DataBuffer.substr(2,1).c_str());
        Int_t Detector = atoi(DataBuffer.substr(18,1).c_str());
        Int_t channel = -1;
        //if (DataBuffer.compare(3,5,"FRONT") == 0) {
        if (DataBuffer.compare(26,5,"FRONT") == 0) {
          //channel = atoi(DataBuffer.substr(7).c_str());
          channel = atoi(DataBuffer.substr(31).c_str());
          *(m_FrontChannelStatus[Detector-1].begin()+channel) = false;
        }
        //else if (DataBuffer.compare(3,4,"BACK") == 0) {
        else if (DataBuffer.compare(26,4,"BACK") == 0) {
          //channel = atoi(DataBuffer.substr(7).c_str());
          channel = atoi(DataBuffer.substr(30).c_str());
          *(m_BackChannelStatus[Detector-1].begin()+channel) = false;
        }
      }

      else if (whatToDo=="TAKE_E_FRONT") {
        m_Take_E_Front = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo=="TAKE_E_BACK") {
        m_Take_E_Front = false;
        cout << whatToDo << endl;
      }

      else if (whatToDo=="STRIP_FRONT_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripFront_E_RAW_Threshold = atoi(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripFront_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="STRIP_BACK_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripBack_E_RAW_Threshold = atoi(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripBack_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="STRIP_FRONT_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripFront_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripFront_E_Threshold << endl;
      }

      else if (whatToDo=="STRIP_BACK_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripBack_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripBack_E_Threshold << endl;
      }


      else {
        ReadingStatus = false;
      }
    }
  }
  cout << "\t/////////////////////////////////////////////////" << endl;
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::Clear()
{
  EventMultiplicity = 0;

  //   Provide a Classification of Event
  EventType.clear();

  // Detector

  //   DSSD
  TowerNumber.clear();
  DetectorNumber.clear();
  Strip_Front.clear();
  Strip_Back.clear();
  //Strip_E.clear();
  Strip_T.clear();
  Front_Energy.clear();
  Back_Energy.clear();
  Half_Energy.clear();
  Front_Time.clear();
  Back_Time.clear();
  //Same_FBTime.clear();

  // counters
  for (Int_t i = 0; i < m_nCounterEvt; i++) m_CounterEvt[i] = 0;
  for (Int_t i = 0; i < m_nCounterHit; i++) m_CounterHit[i] = 0;

  // Calorimeter
  Calor_E.clear();
  Calor_E_calib.clear();
  Calor_T.clear();
  CalorPosX.clear();
  CalorPosY.clear();
  CalorData.clear();

  // all
  deltaT.clear();
}


///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::ReadConfiguration(NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ComptonTelescope");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> token = {"SIZE_DSSSD", "THICKNESS_DSSSD", "NUMBER_STRIPS", "X0_Y0", "X31_Y0", "X0_Y31", "X31_Y31", "NUMBER_DSSSD","DISTANCE_INTER_DSSSD", "DISTANCE_TRACKER_CALORIMETER", "THICKNESS_CALORIMETER", "NPIXELS_CALORIMETER", "TRACKER","CALORIMETER","VIS"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){      
      double     size = blocks[i]->GetDouble("SIZE_DSSSD","mm");
      double thickness = blocks[i]->GetDouble("THICKNESS_DSSSD","mm");
      int    nbr_strip = blocks[i]->GetInt("NUMBER_STRIPS");
      TVector3 A = blocks[i]->GetTVector3("X0_Y0", "mm");
      TVector3 B = blocks[i]->GetTVector3("X31_Y0", "mm");
      TVector3 C = blocks[i]->GetTVector3("X0_Y31", "mm");
      TVector3 D = blocks[i]->GetTVector3("X31_Y31", "mm");
      int     nbr_det = blocks[i]->GetInt("NUMBER_DSSSD");
      double inter = blocks[i]->GetDouble("DISTANCE_INTER_DSSSD","mm");
      double distance_cal = blocks[i]->GetDouble("DISTANCE_TRACKER_CALORIMETER","mm");
      double thickness_cal = blocks[i]->GetDouble("THICKNESS_CALORIMETER","mm");
      int    npixels_cal = blocks[i]->GetInt("NPIXELS_CALORIMETER");
      int    tracker = blocks[i]->GetInt("TRACKER");
      int    calorimeter = blocks[i]->GetInt("CALORIMETER");
      int    vis= blocks[i]->GetInt("VIS");
      AddComptonTelescope(A,B,C,D,size,nbr_strip);
    }

    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  InitializeStandardParameter();
  ReadAnalysisConfig();

/*  // DSSSD blocks
  vector<NPL::InputBlock*> blocksDSSSD = parser.GetAllBlocksWithToken("DSSSD");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocksDSSSD.size() << " DSSSD detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocksDSSSD.size() ; i++){

    vector<string> cartDSSSD = {"X0_Y0", "X31_Y0", "X0_Y31", "X31_Y31"};

    if(blocksDSSSD[i]->HasTokenList(cartDSSSD)){
      cout << endl << "//// DSSSD " << i+1 << endl;
      TVector3 A = blocksDSSSD[i]->GetTVector3("X0_Y0", "mm");
      TVector3 B = blocksDSSSD[i]->GetTVector3("X31_Y0", "mm");
      TVector3 C = blocksDSSSD[i]->GetTVector3("X0_Y31", "mm");
      TVector3 D = blocksDSSSD[i]->GetTVector3("X31_Y31", "mm");

      int nbr_detDSSSD = blocksDSSSD[i]->GetInt("NUMBER_DSSSD");
      double sizeDSSSD = blocksDSSSD[i]->GetDouble("SIZE_DSSSD","mm");
      double thicknessDSSSD = blocksDSSSD[i]->GetDouble("THICKNESS_DSSSD","mm");
      int nbr_stripDSSSD = blocksDSSSD[i]->GetInt("NUMBER_STRIPS");

      AddDetectorDSSSD(A,B,C,D,sizeDSSSD,nbr_stripDSSSD);
    }
    else {
      cout << "ERROR: Missing token for DSSSD blocks, check your input file" << endl;
      exit(1);
    }
  } // end DSSSD blocks
  
  // Calorimeter blocks
  vector<NPL::InputBlock*> blocksCalo = parser.GetAllBlocksWithToken("CALORIMETER");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
  cout <<  "//// " << blocksCalo.size() << " calorimeter detectors found " << endl;

  for (unsigned int i  = 0 ; i < blocksCalo.size() ; i++){
    int nbr_detCalorimeter = blocksCalo[i]->GetInt("NUMBER_CALORIMETER");
    double thickness_cal = blocksCalo[i]->GetDouble("THICKNESS_CALORIMETER","mm");
    int npixels_cal = blocksCalo[i]->GetInt("NPIXELS_CALORIMETER");
  } // end calorimeter blocks

  // ComptonTelescope blocks
  vector<NPL::InputBlock*> blocksCompton = parser.GetAllBlocksWithToken("ComptonTelescope");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout <<  "//// " << blocksCompton.size() << " ComptonTelescope found " << endl;

  for (unsigned int i  = 0 ; i < blocksCompton.size() ; i++){
    int nb_tracker = blocksCompton[i]->GetInt("TRACKER");
    double inter = blocksCompton[i]->GetDouble("DISTANCE_INTER_DSSSD","mm");
    double distance_cal = blocksCompton[i]->GetDouble("DISTANCE_TRACKER_CALORIMETER","mm");
    int vis= blocksCompton[i]->GetInt("VIS");
  }
*/

}

///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::AddParameterToCalibrationManager()
{
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    for (int j = 0; j < m_NumberOfStrips; ++j) {
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j)+"_E", "COMPTONTELESCOPE_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j)+"_E");
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j)+"_T", "COMPTONTELESCOPE_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j)+"_T");
    }
    for (int j = 0; j < m_NumberOfStrips; ++j) {
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j)+"_E",  "COMPTONTELESCOPE_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j)+"_E");
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j)+"_T",  "COMPTONTELESCOPE_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j)+"_T");
    }
    for (int j = 0; j < m_NPixels; ++j) {
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+1)+"_CHANNEL"+ NPL::itoa(j)+"_E",  "COMPTONTELESCOPE_D"+ NPL::itoa(i+1)+"_CHANNEL"+ NPL::itoa(j)+"_E");
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+3)+"_CHANNEL"+ NPL::itoa(j)+"_PED",  "COMPTONTELESCOPE_D"+ NPL::itoa(i+3)+"_CHANNEL"+ NPL::itoa(j)+"_PED");
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+4)+"_CHANNEL"+ NPL::itoa(j)+"_PED",  "COMPTONTELESCOPE_D"+ NPL::itoa(i+4)+"_CHANNEL"+ NPL::itoa(j)+"_PED");
      Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+4)+"_CHANNEL"+ NPL::itoa(j)+"_COEFFALIGN",  "COMPTONTELESCOPE_D"+ NPL::itoa(i+4)+"_CHANNEL"+ NPL::itoa(j)+"_COEFFALIGN");
    }
    Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+1)+"_Q_E", "COMPTONTELESCOPE_D"+ NPL::itoa(i+1)+"_Q_E");
    Cal->AddParameter("COMPTONTELESCOPE", "D"+ NPL::itoa(i+4)+"_E_CALO", "COMPTONTELESCOPE_D"+ NPL::itoa(i+4)+"_E_CALO");
  }

  return;  
}


///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::InitializeRootInputRaw()
{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("ComptonTelescope", true);
  inputChain->SetBranchStatus("fCT_*", true);
  inputChain->SetBranchAddress("ComptonTelescope", &m_EventData);
}


///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::InitializeRootInputPhysics()
{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("EventMultiplicity", true);
  inputChain->SetBranchStatus("EventType",         true);
  inputChain->SetBranchStatus("TowerNumber",    true);
  inputChain->SetBranchStatus("DetectorNumber",    true);
  inputChain->SetBranchStatus("Strip_Front",       true);
  inputChain->SetBranchStatus("Strip_Back",        true);
  //inputChain->SetBranchStatus("Strip_E",           true);
  inputChain->SetBranchStatus("Strip_T",           true);
  inputChain->SetBranchStatus("Front_Energy",      true);
  inputChain->SetBranchStatus("Back_Energy",      true);
  inputChain->SetBranchStatus("Half_Energy",      true);
  inputChain->SetBranchStatus("Front_Time",      true);
  inputChain->SetBranchStatus("Back_Time",       true);
  //inputChain->SetBranchStatus("Same_FBTime",     true);
  inputChain->SetBranchStatus("Calor_E",        true);
  inputChain->SetBranchStatus("Calor_E_calib",  true);
  inputChain->SetBranchStatus("Calor_T",        true);
  inputChain->SetBranchStatus("CalorPosX",      true);
  inputChain->SetBranchStatus("CalorPosY",      true);
  inputChain->SetBranchStatus("CalorData",      true);
}


///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::InitializeRootOutput()
{
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("ComptonTelescope", "TComptonTelescopePhysics", &m_EventPhysics);
}


/////   Specific to ComptonTelescopeArray   ////

//void TComptonTelescopePhysics::AddDetectorDSSSD(TVector3 C_X0_Y0, TVector3 C_X31_Y0, TVector3 C_X0_Y31, TVector3 C_X31_Y31, double size_dsssd, int nb_strip)
void TComptonTelescopePhysics::AddComptonTelescope(TVector3 C_X0_Y0, TVector3 C_X31_Y0, TVector3 C_X0_Y31, TVector3 C_X31_Y31, double size_dsssd, int nb_strip)
{
  m_NumberOfDetectors++;

  // remove warning using C_X31_Y31
  C_X31_Y31.Unit();

  // Vector U on Module Face (parallele to Y/Back Strip)
  // NB: Y strips are allong X axis
  TVector3 U = C_X31_Y0 - C_X0_Y0;
  U = U.Unit();

  // Vector V on Module Face (parallele to X Strip)
  TVector3 V = C_X0_Y31 - C_X0_Y0;
  V = V.Unit();

  // Position Vector of Strip Center
  TVector3 StripCenter = TVector3(0,0,0);

  // Position Vector of X=0 Y=0 strip
  TVector3 Strip_0_0;

  // Buffer object to fill Position Array
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector< vector< double > >   OneModuleStripPositionX;
  vector< vector< double > >   OneModuleStripPositionY;
  vector< vector< double > >   OneModuleStripPositionZ;

  // strip pitch
  double stripPitch = size_dsssd/(double)nb_strip;
  // Moving StripCenter to strip center of 0.0 corner
  Strip_0_0 = C_X0_Y0 + (U+V) * (stripPitch/2.);

  for (int i = 0; i < nb_strip; i++) {
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < nb_strip; j++) {
      StripCenter = Strip_0_0 + stripPitch*(i*U + j*V);
      lineX.push_back( StripCenter.X() );
      lineY.push_back( StripCenter.Y() );
      lineZ.push_back( StripCenter.Z() );
    }

    OneModuleStripPositionX.push_back(lineX);
    OneModuleStripPositionY.push_back(lineY);
    OneModuleStripPositionZ.push_back(lineZ);
  }

  m_StripPositionX.push_back(OneModuleStripPositionX);
  m_StripPositionY.push_back(OneModuleStripPositionY);
  m_StripPositionZ.push_back(OneModuleStripPositionZ);
}

void TComptonTelescopePhysics::AddComptonTelescope(double Z)
{
  m_NumberOfDetectors++;
  // empty at the moment
  // needed if solid angle analysis are needed
}

TVector3 TComptonTelescopePhysics::GetDetectorNormal( const int i) const{
  /*  TVector3 U =    TVector3 ( GetStripPositionX( DetectorNumber[i] , 24 , 1 ) ,
      GetStripPositionY( DetectorNumber[i] , 24 , 1 ) ,
      GetStripPositionZ( DetectorNumber[i] , 24 , 1 ) )

      -TVector3 ( GetStripPositionX( DetectorNumber[i] , 1 , 1 ) ,
      GetStripPositionY( DetectorNumber[i] , 1 , 1 ) ,
      GetStripPositionZ( DetectorNumber[i] , 1 , 1 ) );

      TVector3 V =    TVector3 ( GetStripPositionX( DetectorNumber[i] , 24 , 48 ) ,
      GetStripPositionY( DetectorNumber[i] , 24 , 48 ) ,
      GetStripPositionZ( DetectorNumber[i] , 24 , 48 ) )

      -TVector3 ( GetStripPositionX( DetectorNumber[i] , 24 , 1 ) ,
      GetStripPositionY( DetectorNumber[i] , 24 , 1 ) ,
      GetStripPositionZ( DetectorNumber[i] , 24 , 1 ) );

      TVector3 Normal = U.Cross(V);

      return(Normal.Unit()) ;*/

  return (TVector3(0,0,i));

}

TVector3 TComptonTelescopePhysics::GetPositionOfInteractionDSSSD(const int i) const{
  TVector3    Position = TVector3 (  GetStripPositionX( DetectorNumber[i] , Strip_Front[i] , Strip_Back[i] )    ,
      GetStripPositionY( DetectorNumber[i] , Strip_Front[i] , Strip_Back[i] )    ,
      GetStripPositionZ( DetectorNumber[i] , Strip_Front[i] , Strip_Back[i] )    ) ;

  return(Position) ;

}

/*double TComptonTelescopePhysics::GetCalor_E() {
  return Calor_E;
}*/


///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::InitializeStandardParameter()
{
  // Enable all channels
  vector<bool> ChannelStatus;
  m_FrontChannelStatus.clear();
  m_BackChannelStatus.clear();

  ChannelStatus.resize(m_NumberOfStrips, true);
  for(Int_t i = 0; i < m_NumberOfDetectors; ++i) {
    m_FrontChannelStatus[i] = ChannelStatus;
    m_BackChannelStatus[i]  = ChannelStatus;
  }

  m_MaximumStripMultiplicityAllowed = m_NumberOfDetectors;
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::InitSpectra()
{
  m_Spectra = new TComptonTelescopeSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::FillSpectra()
{
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::CheckSpectra()
{
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::ClearSpectra()
{
  // m_Spectra -> Clear();
  // To be done
}



///////////////////////////////////////////////////////////////////////////
void TComptonTelescopePhysics::WriteSpectra()
{
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TComptonTelescopePhysics::GetSpectra() 
{
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}



///////////////////////////////////////////////////////////////////////////
namespace ComptonTelescope_LOCAL
{
  // DSSD
  // Front
  Double_t fStrip_Front_E(const TComptonTelescopeData* m_EventData, const int i)
  {
    return CalibrationManager::getInstance()->ApplyCalibration("COMPTONTELESCOPE/D" + NPL::itoa(m_EventData->GetCTTrackerFrontEDetectorNbr(i)) + "_STRIP_FRONT" + NPL::itoa(m_EventData->GetCTTrackerFrontEStripNbr(i)) + "_E", m_EventData->GetCTTrackerFrontEEnergy(i));
  }

  // Back
  Double_t fStrip_Back_E(const TComptonTelescopeData* m_EventData, const int i)
  {
    return CalibrationManager::getInstance()->ApplyCalibration("COMPTONTELESCOPE/D" + NPL::itoa(m_EventData->GetCTTrackerBackEDetectorNbr(i)) + "_STRIP_BACK" + NPL::itoa(m_EventData->GetCTTrackerBackEStripNbr(i)) + "_E", m_EventData->GetCTTrackerBackEEnergy(i));
  }

  //Calorimeter
  Double_t fCalorimeter_ped(const TComptonTelescopeData* m_EventData, const int i)
  {
    return CalibrationManager::getInstance()->ApplyCalibration("COMPTONTELESCOPE/D" + NPL::itoa(m_EventData->GetCTCalorimeterEDetectorNbr(i)) + "_CHANNEL" + NPL::itoa(m_EventData->GetCTCalorimeterEChannelNbr(i)) + "_PED", m_EventData->GetCTCalorimeterEEnergy(i));
  }

  Double_t fCalorimeter_coeffAlign(const TComptonTelescopeData* m_EventData, const int i) // coeff to correct gain dispersion between pixels ("calibration de Jean")
  {
    return CalibrationManager::getInstance()->ApplyCalibration("COMPTONTELESCOPE/D" + NPL::itoa(m_EventData->GetCTCalorimeterEDetectorNbr(i)) + "_CHANNEL" + NPL::itoa(m_EventData->GetCTCalorimeterEChannelNbr(i)) + "_COEFFALIGN", m_EventData->GetCTCalorimeterEEnergy(i));
  }

  Double_t fCalorimeter_calibE(double sumADC, int detectorNumber) // ADC - Energy calibration
  {
    return CalibrationManager::getInstance()->ApplyCalibration("COMPTONTELESCOPE/D" + NPL::itoa(detectorNumber) + "_E_CALO", sumADC);
  }

 Double_t fCalorimeter_Q(const TComptonTelescopeData* m_EventData, const int i) // Charge-ADC calibration
  {
    return CalibrationManager::getInstance()->ApplyCalibration("COMPTONTELESCOPE/D" + NPL::itoa(m_EventData->GetCTCalorimeterEDetectorNbr(i)) + "_CHANNEL" + NPL::itoa(m_EventData->GetCTCalorimeterEChannelNbr(i)) + "_E", m_EventData->GetCTCalorimeterEEnergy(i));
  }

  Double_t fCalorimeter_E(double charge, int detectorNumber) // Total charge-energy relation
  {
    return CalibrationManager::getInstance()->ApplyCalibration("COMPTONTELESCOPE/D" + NPL::itoa(detectorNumber) + "_Q_E", charge);
  }
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TComptonTelescopePhysics::Construct(){
  return (NPL::VDetector*) new TComptonTelescopePhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_comptontelescope{
  public:
    proxy_comptontelescope(){
      NPL::DetectorFactory::getInstance()->AddToken("ComptonTelescope","ComptonTelescope");
      NPL::DetectorFactory::getInstance()->AddDetector("ComptonTelescope",TComptonTelescopePhysics::Construct);
    }
};

proxy_comptontelescope p;
}

