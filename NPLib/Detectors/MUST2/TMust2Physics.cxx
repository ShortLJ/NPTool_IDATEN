/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : febuary 2009                                             *
 * Last update    : July 2021
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold must2 treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TMust2Physics.h"
using namespace MUST2_LOCAL;

//   STL
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdlib.h>

//   NPL
#include "NPDetectorFactory.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
using namespace NPUNITS;
//   ROOT
#include "TChain.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TMust2Physics)

    ///////////////////////////////////////////////////////////////////////////
    TMust2Physics::TMust2Physics() {
  EventMultiplicity                  = 0;
  m_EventData                        = new TMust2Data;
  m_PreTreatedData                   = new TMust2Data;
  m_EventPhysics                     = this;
  m_Spectra                          = NULL;
  m_NumberOfTelescope                = 0;
  m_MaximumStripMultiplicityAllowed  = 10;
  m_StripEnergyMatchingSigma         = 0.020;
  m_StripEnergyMatchingNumberOfSigma = 3;
  // Raw Threshold
  m_Si_X_E_RAW_Threshold = 8192;
  m_Si_Y_E_RAW_Threshold = 8192;
  m_SiLi_E_RAW_Threshold = 8192;
  m_CsI_E_RAW_Threshold  = 8192;
  // Calibrated Threshold
  m_Si_X_E_Threshold = 0;
  m_Si_Y_E_Threshold = 0;
  m_SiLi_E_Threshold = 0;
  m_CsI_E_Threshold  = 0;

  m_Ignore_not_matching_SiLi = false;
  m_Ignore_not_matching_CsI  = false;

  m_Take_E_Y = false;
  m_Take_T_Y = true;

  //////////////////
  // SiLi matching //
  //////////////////

  m_SiLi_Size = 32;
  m_SiLi_MatchingX.resize(16, 0);
  m_SiLi_MatchingY.resize(16, 0);

  m_SiLi_MatchingX[0] = 112;
  m_SiLi_MatchingY[0] = 112;

  m_SiLi_MatchingX[1] = 112;
  m_SiLi_MatchingY[1] = 80;

  m_SiLi_MatchingX[2] = 80;
  m_SiLi_MatchingY[2] = 112;

  m_SiLi_MatchingX[3] = 80;
  m_SiLi_MatchingY[3] = 80;
  //
  m_SiLi_MatchingX[4] = 48;
  m_SiLi_MatchingY[4] = 80;

  m_SiLi_MatchingX[5] = 48;
  m_SiLi_MatchingY[5] = 112;

  m_SiLi_MatchingX[6] = 16;
  m_SiLi_MatchingY[6] = 80;

  m_SiLi_MatchingX[7] = 16;
  m_SiLi_MatchingY[7] = 112;
  //
  m_SiLi_MatchingX[8] = 112;
  m_SiLi_MatchingY[8] = 16;

  m_SiLi_MatchingX[9] = 112;
  m_SiLi_MatchingY[9] = 48;

  m_SiLi_MatchingX[10] = 80;
  m_SiLi_MatchingY[10] = 16;

  m_SiLi_MatchingX[11] = 80;
  m_SiLi_MatchingY[11] = 48;
  //
  m_SiLi_MatchingX[12] = 48;
  m_SiLi_MatchingY[12] = 48;

  m_SiLi_MatchingX[13] = 48;
  m_SiLi_MatchingY[13] = 16;

  m_SiLi_MatchingX[14] = 16;
  m_SiLi_MatchingY[14] = 48;

  m_SiLi_MatchingX[15] = 16;
  m_SiLi_MatchingY[15] = 16;

  //////////////////
  // CsI matching //
  //////////////////

  m_CsI_Size = 32;
  m_CsI_MatchingX.resize(16, 0);
  m_CsI_MatchingY.resize(16, 0);

  m_CsI_MatchingX[0] = 112;
  m_CsI_MatchingY[0] = 112;

  m_CsI_MatchingX[1] = 112;
  m_CsI_MatchingY[1] = 80;

  m_CsI_MatchingX[2] = 112;
  m_CsI_MatchingY[2] = 48;

  m_CsI_MatchingX[3] = 112;
  m_CsI_MatchingY[3] = 16;
  //
  m_CsI_MatchingX[4] = 80;
  m_CsI_MatchingY[4] = 16;

  m_CsI_MatchingX[5] = 80;
  m_CsI_MatchingY[5] = 48;

  m_CsI_MatchingX[6] = 80;
  m_CsI_MatchingY[6] = 80;

  m_CsI_MatchingX[7] = 80;
  m_CsI_MatchingY[7] = 112;
  //
  m_CsI_MatchingX[8] = 48;
  m_CsI_MatchingY[8] = 16;

  m_CsI_MatchingX[9] = 48;
  m_CsI_MatchingY[9] = 48;

  m_CsI_MatchingX[10] = 48;
  m_CsI_MatchingY[10] = 80;

  m_CsI_MatchingX[11] = 48;
  m_CsI_MatchingY[11] = 112;
  //
  m_CsI_MatchingX[12] = 16;
  m_CsI_MatchingY[12] = 16;

  m_CsI_MatchingX[13] = 16;
  m_CsI_MatchingY[13] = 48;

  m_CsI_MatchingX[14] = 16;
  m_CsI_MatchingY[14] = 80;

  m_CsI_MatchingX[15] = 16;
  m_CsI_MatchingY[15] = 112;

  m_CsI_MatchingX[0] = 112;
  m_CsI_MatchingY[0] = 112;

  m_CsI_MatchingX[1] = 112;
  m_CsI_MatchingY[1] = 80;

  m_CsI_MatchingX[2] = 112;
  m_CsI_MatchingY[2] = 48;

  m_CsI_MatchingX[3] = 112;
  m_CsI_MatchingY[3] = 16;
  //
  m_CsI_MatchingX[4] = 80;
  m_CsI_MatchingY[4] = 16;

  m_CsI_MatchingX[5] = 80;
  m_CsI_MatchingY[5] = 48;

  m_CsI_MatchingX[6] = 80;
  m_CsI_MatchingY[6] = 80;

  m_CsI_MatchingX[7] = 80;
  m_CsI_MatchingY[7] = 112;
  //
  m_CsI_MatchingX[8] = 48;
  m_CsI_MatchingY[8] = 16;

  m_CsI_MatchingX[9] = 48;
  m_CsI_MatchingY[9] = 48;

  m_CsI_MatchingX[10] = 48;
  m_CsI_MatchingY[10] = 80;

  m_CsI_MatchingX[11] = 48;
  m_CsI_MatchingY[11] = 112;
  //
  m_CsI_MatchingX[12] = 16;
  m_CsI_MatchingY[12] = 16;

  m_CsI_MatchingX[13] = 16;
  m_CsI_MatchingY[13] = 48;

  m_CsI_MatchingX[14] = 16;
  m_CsI_MatchingY[14] = 80;

  // FIXME: Temporary fix... for zero degree
  // Particular strip matching required for zero degree telescope...
  // Might be a geometrical problem, only relevant if a telescope is placed at
  // zero degree
  // m_ZeroDegree_CsI_MatchingX[0].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[0].second = 128;
  // m_ZeroDegree_CsI_MatchingY[0].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[0].second = 128;

  // m_ZeroDegree_CsI_MatchingX[1].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[1].second = 128;
  // m_ZeroDegree_CsI_MatchingY[1].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[1].second = 102;

  // m_ZeroDegree_CsI_MatchingX[2].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[2].second = 128;
  // m_ZeroDegree_CsI_MatchingY[2].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[2].second = 64;

  // m_ZeroDegree_CsI_MatchingX[3].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[3].second = 128;
  // m_ZeroDegree_CsI_MatchingY[3].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[3].second = 26;

  // m_ZeroDegree_CsI_MatchingX[4].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[4].second = 102;
  // m_ZeroDegree_CsI_MatchingY[4].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[4].second = 26;

  // m_ZeroDegree_CsI_MatchingX[5].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[5].second = 102;
  // m_ZeroDegree_CsI_MatchingY[5].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[5].second = 64;

  // m_ZeroDegree_CsI_MatchingX[6].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[6].second = 102;
  // m_ZeroDegree_CsI_MatchingY[6].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[6].second = 102;

  // m_ZeroDegree_CsI_MatchingX[7].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[7].second = 102;
  // m_ZeroDegree_CsI_MatchingY[7].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[7].second = 128;

  // m_ZeroDegree_CsI_MatchingX[8].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[8].second = 64;
  // m_ZeroDegree_CsI_MatchingY[8].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[8].second = 26;

  // m_ZeroDegree_CsI_MatchingX[9].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[9].second = 64;
  // m_ZeroDegree_CsI_MatchingY[9].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[9].second = 64;

  // m_ZeroDegree_CsI_MatchingX[10].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[10].second = 64;
  // m_ZeroDegree_CsI_MatchingY[10].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[10].second = 102;

  // m_ZeroDegree_CsI_MatchingX[11].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[11].second = 64;
  // m_ZeroDegree_CsI_MatchingY[11].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[11].second = 128;

  // m_ZeroDegree_CsI_MatchingX[12].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[12].second = 26;
  // m_ZeroDegree_CsI_MatchingY[12].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[12].second = 26;

  // m_ZeroDegree_CsI_MatchingX[13].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[13].second = 26;
  // m_ZeroDegree_CsI_MatchingY[13].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[13].second = 64;

  // m_ZeroDegree_CsI_MatchingX[14].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[14].second = 26;
  // m_ZeroDegree_CsI_MatchingY[14].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[14].second = 102;

  // m_ZeroDegree_CsI_MatchingX[15].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[15].second = 26;
  // m_ZeroDegree_CsI_MatchingY[15].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[15].second = 128;
}

///////////////////////////////////////////////////////////////////////////
TMust2Physics::~TMust2Physics() {}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

///////////////////////////////////////////////////////////////////////////

void TMust2Physics::BuildPhysicalEvent() {

  PreTreat();

  // FIXME has to be set in the configuration file
  m_multimatch = false;

  bool check_SILI = false;
  bool check_CSI  = false;

  m_StripXEMult = m_PreTreatedData->GetMMStripXEMult();
  m_StripYEMult = m_PreTreatedData->GetMMStripYEMult();
  m_StripXTMult = m_PreTreatedData->GetMMStripXTMult();
  m_StripYTMult = m_PreTreatedData->GetMMStripYTMult();
  m_SiLiEMult   = m_PreTreatedData->GetMMSiLiEMult();
  m_SiLiTMult   = m_PreTreatedData->GetMMSiLiTMult();
  m_CsIEMult    = m_PreTreatedData->GetMMCsIEMult();
  m_CsITMult    = m_PreTreatedData->GetMMCsITMult();

  /////////////////////////////////////////////////
  // Returns the matching couples and set the type
  // of matching. In case of a multiple match in one
  // detector, m_match_type is set to 1, else it is
  // set to 2
  vector<TVector2> couple = Match_X_Y();
  /////////////////////////////////////////////////

  unsigned int couple_size = couple.size();

  for (unsigned int i = 0; i < couple_size; ++i) {

    // if (m_match_type[i] == 1 || m_multimatch) {
    if (m_match_type[i] == 1) {
      check_SILI = false;
      check_CSI  = false;

      int couple_X = couple[i].X() + 0.5;
      int couple_Y = couple[i].Y() + 0.5;
      int N        = m_PreTreatedData->GetMMStripXEDetectorNbr(couple_X);

      int X = m_PreTreatedData->GetMMStripXEStripNbr(couple_X);
      int Y = m_PreTreatedData->GetMMStripYEStripNbr(couple_Y);

      double Si_X_E = m_PreTreatedData->GetMMStripXEEnergy(couple_X);
      double Si_Y_E = m_PreTreatedData->GetMMStripYEEnergy(couple_Y);

      //  Search for associate Time
      double Si_X_T = -1000;
      for (unsigned int t = 0; t < m_StripXTMult; ++t) {
        if (m_PreTreatedData->GetMMStripXTStripNbr(couple_X)
            == m_PreTreatedData->GetMMStripXTStripNbr(t)) {
          Si_X_T = m_PreTreatedData->GetMMStripXTTime(t);
          break;
        }
      }

      double Si_Y_T = -1000;
      for (unsigned int t = 0; t < m_StripYTMult; ++t) {
        if (m_PreTreatedData->GetMMStripYTStripNbr(couple_Y)
            == m_PreTreatedData->GetMMStripYTStripNbr(t)) {
          Si_Y_T = m_PreTreatedData->GetMMStripYTTime(t);
          break;
        }
      }

      for (unsigned int j = 0; j < m_SiLiEMult; ++j) {
        if (m_PreTreatedData->GetMMSiLiEDetectorNbr(j) == N) {
          // pad vs strip number match
          if (Match_Si_SiLi(X, Y, m_PreTreatedData->GetMMSiLiEPadNbr(j))) {
            SiLi_N.push_back(m_PreTreatedData->GetMMSiLiEPadNbr(j));
            SiLi_E.push_back(m_PreTreatedData->GetMMSiLiEEnergy(j));
            SiLi_T.push_back(-1000);
            // Look for associate time
            for (unsigned int k = 0; k < m_SiLiTMult; ++k) {
              // Same Pad, same Detector
              if (m_PreTreatedData->GetMMSiLiEPadNbr(j)
                  == m_PreTreatedData->GetMMSiLiTPadNbr(k)) {
                SiLi_T[SiLi_T.size() - 1] = m_PreTreatedData->GetMMSiLiTTime(k);
                break;
              }
            }
            check_SILI = true;
          }
        }
      }

      for (unsigned int j = 0; j < m_CsIEMult; ++j) {
        if (m_PreTreatedData->GetMMCsIEDetectorNbr(j) == N) {
          if (Match_Si_CsI(X, Y, m_PreTreatedData->GetMMCsIECristalNbr(j),
                           m_PreTreatedData->GetMMCsIEDetectorNbr(j))) {
            CsI_N.push_back(m_PreTreatedData->GetMMCsIECristalNbr(j));
            CsI_E.push_back(m_PreTreatedData->GetMMCsIEEnergy(j));
            CsI_T.push_back(-1000);
            // Look for associate Time
            for (unsigned int k = 0; k < m_CsITMult; ++k) {
              // Same Cristal, Same Detector
              if (m_PreTreatedData->GetMMCsIECristalNbr(j)
                  == m_PreTreatedData->GetMMCsITCristalNbr(k)) {
                CsI_T[CsI_T.size() - 1] = m_PreTreatedData->GetMMCsITTime(j);
                break;
              }
            }
            check_CSI = true;
          }
        }
      }

      /////////////////////////////////////////////////

      TelescopeNumber.push_back(N);
      Si_X.push_back(X);
      Si_Y.push_back(Y);

      // Store the other value for checking purpose
      Si_EX.push_back(Si_X_E);
      Si_TX.push_back(Si_X_T);

      Si_EY.push_back(Si_Y_E);
      Si_TY.push_back(Si_Y_T);

      if (m_Take_E_Y)
        Si_E.push_back(Si_Y_E);
      else
        Si_E.push_back(Si_X_E);

      if (m_Take_T_Y)
        Si_T.push_back(Si_Y_T);
      else
        Si_T.push_back(Si_X_T);

      if (!check_SILI) {
        SiLi_N.push_back(0);
        SiLi_E.push_back(-1000);
        SiLi_T.push_back(-1000);
      }

      if (!check_CSI) {
        CsI_N.push_back(0);
        CsI_E.push_back(-1000);
        CsI_T.push_back(-1000);
      }

      EventType.push_back(m_match_type[i]);
    }
  } // loop on event multiplicity
  EventMultiplicity = TelescopeNumber.size();

  return;
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::PreTreat() {
  ClearPreTreatedData();
  m_StripXEMult = m_EventData->GetMMStripXEMult();
  m_StripYEMult = m_EventData->GetMMStripYEMult();
  m_StripXTMult = m_EventData->GetMMStripXTMult();
  m_StripYTMult = m_EventData->GetMMStripYTMult();
  m_SiLiEMult   = m_EventData->GetMMSiLiEMult();
  m_SiLiTMult   = m_EventData->GetMMSiLiTMult();
  m_CsIEMult    = m_EventData->GetMMCsIEMult();
  m_CsITMult    = m_EventData->GetMMCsITMult();
  map<int, int> hitX;
  map<int, int> hitY;

  //   X
  //   E
  for (unsigned int i = 0; i < m_StripXEMult; ++i) {
    if (m_EventData->GetMMStripXEEnergy(i) > m_Si_X_E_RAW_Threshold
        && IsValidChannel(0, m_EventData->GetMMStripXEDetectorNbr(i),
                          m_EventData->GetMMStripXEStripNbr(i))) {
      double EX = fSi_X_E(m_EventData, i);
      if (EX > m_Si_X_E_Threshold)
        m_PreTreatedData->SetStripXE(m_EventData->GetMMStripXEDetectorNbr(i),
                                     m_EventData->GetMMStripXEStripNbr(i), EX);
    }
  }

  //   T
  for (unsigned int i = 0; i < m_StripXTMult; ++i) {
    if (IsValidChannel(0, m_EventData->GetMMStripXTDetectorNbr(i),
                       m_EventData->GetMMStripXTStripNbr(i)))
      m_PreTreatedData->SetStripXT(m_EventData->GetMMStripXTDetectorNbr(i),
                                   m_EventData->GetMMStripXTStripNbr(i),
                                   fSi_X_T(m_EventData, i));
  }

  //   Y
  //   E
  for (unsigned int i = 0; i < m_StripYEMult; ++i) {
    if (m_EventData->GetMMStripYEEnergy(i) < m_Si_Y_E_RAW_Threshold
        && IsValidChannel(1, m_EventData->GetMMStripYEDetectorNbr(i),
                          m_EventData->GetMMStripYEStripNbr(i))) {
      double EY = fSi_Y_E(m_EventData, i);
      if (EY > m_Si_Y_E_Threshold)
        m_PreTreatedData->SetStripYE(m_EventData->GetMMStripYEDetectorNbr(i),
                                     m_EventData->GetMMStripYEStripNbr(i), EY);
    }
  }

  //   T
  for (unsigned int i = 0; i < m_StripYTMult; ++i) {
    if (IsValidChannel(1, m_EventData->GetMMStripYTDetectorNbr(i),
                       m_EventData->GetMMStripYTStripNbr(i)))
      m_PreTreatedData->SetStripYT(m_EventData->GetMMStripYTDetectorNbr(i),
                                   m_EventData->GetMMStripYTStripNbr(i),
                                   fSi_Y_T(m_EventData, i));
  }

  //   CsI
  //   E
  for (unsigned int i = 0; i < m_CsIEMult; ++i) {
    if (m_EventData->GetMMCsIEEnergy(i) > m_CsI_E_RAW_Threshold
        && IsValidChannel(3, m_EventData->GetMMCsIEDetectorNbr(i),
                          m_EventData->GetMMCsIECristalNbr(i))) {
      double ECsI = fCsI_E(m_EventData, i);
      if (ECsI > m_CsI_E_Threshold)
        m_PreTreatedData->SetCsIE(m_EventData->GetMMCsIEDetectorNbr(i),
                                  m_EventData->GetMMCsIECristalNbr(i), ECsI);
    }
  }

  //   T
  for (unsigned int i = 0; i < m_CsITMult; ++i) {
    if (IsValidChannel(3, m_EventData->GetMMCsITDetectorNbr(i),
                       m_EventData->GetMMCsITCristalNbr(i)))
      m_PreTreatedData->SetCsIT(m_EventData->GetMMCsITDetectorNbr(i),
                                m_EventData->GetMMCsITCristalNbr(i),
                                fCsI_T(m_EventData, i));
  }

  //   SiLi
  //   E
  for (unsigned int i = 0; i < m_SiLiEMult; ++i) {
    if (m_EventData->GetMMSiLiEEnergy(i) > m_SiLi_E_RAW_Threshold
        && IsValidChannel(2, m_EventData->GetMMSiLiEDetectorNbr(i),
                          m_EventData->GetMMSiLiEPadNbr(i))) {
      double ESiLi = fSiLi_E(m_EventData, i);
      if (ESiLi > m_SiLi_E_Threshold)
        m_PreTreatedData->SetSiLiE(m_EventData->GetMMSiLiEDetectorNbr(i),
                                   m_EventData->GetMMSiLiEPadNbr(i), ESiLi);
    }
  }

  //   T
  for (unsigned int i = 0; i < m_SiLiTMult; ++i) {
    if (IsValidChannel(2, m_EventData->GetMMSiLiTDetectorNbr(i),
                       m_EventData->GetMMSiLiTPadNbr(i)))
      m_PreTreatedData->SetSiLiT(m_EventData->GetMMSiLiTDetectorNbr(i),
                                 m_EventData->GetMMSiLiTPadNbr(i),
                                 fSiLi_T(m_EventData, i));
  }

  return;
}

///////////////////////////////////////////////////////////////////////////
bool TMust2Physics::ResolvePseudoEvent() { return false; }

///////////////////////////////////////////////////////////////////////////

vector<TVector2> TMust2Physics::Match_X_Y() {
  vector<TVector2> ArrayOfGoodCouple;
  ArrayOfGoodCouple.clear();

  m_StripXEMult = m_PreTreatedData->GetMMStripXEMult();
  m_StripYEMult = m_PreTreatedData->GetMMStripYEMult();

  double matchSigma  = m_StripEnergyMatchingSigma;
  double NmatchSigma = m_StripEnergyMatchingNumberOfSigma;

  // Prevent code from treating very high multiplicity Event
  // Those event are not physical anyway and that improve speed.
  if (m_StripXEMult > m_MaximumStripMultiplicityAllowed
      || m_StripYEMult > m_MaximumStripMultiplicityAllowed) {
    return ArrayOfGoodCouple;
  }

  // Get Detector multiplicity
  for (unsigned int i = 0; i < m_StripXEMult; i++) {
    int N = m_PreTreatedData->GetMMStripXEDetectorNbr(i);
    m_StripXMultDet[N] += 1;
  }

  for (unsigned int j = 0; j < m_StripYEMult; j++) {
    int N = m_PreTreatedData->GetMMStripYEDetectorNbr(j);
    m_StripYMultDet[N] += 1;
  }
  for (unsigned int i = 0; i < m_StripXEMult; i++) {
    for (unsigned int j = 0; j < m_StripYEMult; j++) {

      // Declaration of variable for clarity
      int StripXDetNbr = m_PreTreatedData->GetMMStripXEDetectorNbr(i);
      int StripYDetNbr = m_PreTreatedData->GetMMStripYEDetectorNbr(j);

      //   if same detector check energy
      if (StripXDetNbr == StripYDetNbr) {

        int DetNbr = StripXDetNbr;

        // Declaration of variable for clarity
        double StripXEnergy = m_PreTreatedData->GetMMStripXEEnergy(i);
        double StripXNbr    = m_PreTreatedData->GetMMStripXEStripNbr(i);

        double StripYEnergy = m_PreTreatedData->GetMMStripYEEnergy(j);
        double StripYNbr    = m_PreTreatedData->GetMMStripYEStripNbr(j);

        //   Look if energy match
        if (abs((StripXEnergy - StripYEnergy) / 2.)
            < NmatchSigma * matchSigma) {

          // Special Option, if the event is between two CsI
          // cristal, it is rejected.
          if (m_Ignore_not_matching_CsI) {
            bool check_validity = false;
            for (unsigned int hh = 0; hh < 16; ++hh) {
              if (Match_Si_CsI(StripXNbr, StripYNbr, hh + 1, DetNbr)) {
                check_validity = true;
              }
            }
            if (check_validity) {
              ArrayOfGoodCouple.push_back(TVector2(i, j));
            }
          }

          // Special Option, if the event is between two SiLi pad ,
          // it is rejected.
          else if (m_Ignore_not_matching_SiLi) {
            bool check_validity = false;
            for (unsigned int hh = 0; hh < 16; ++hh) {
              if (Match_Si_SiLi(StripXNbr, StripYNbr, hh + 1))
                check_validity = true;
            }
            if (check_validity)
              ArrayOfGoodCouple.push_back(TVector2(i, j));
          } else {
            // Regular case, keep the event
            ArrayOfGoodCouple.push_back(TVector2(i, j));
            m_NMatchX[i] += 1;
            m_NMatchY[j] += 1;

            m_NMatchDet[DetNbr] += 1;
          }
        }
      } // if same detector
    } // loop on StripY Mult
  } // loop on StripX Mult

  unsigned int couple_size = ArrayOfGoodCouple.size();
  for (unsigned int i = 0; i < couple_size; ++i) {
    int N = m_PreTreatedData->GetMMStripXEDetectorNbr(ArrayOfGoodCouple[i].X());
    int Xi = ArrayOfGoodCouple[i].X();
    int Yj = ArrayOfGoodCouple[i].Y();
    if (m_NMatchX[Xi] > 1 || m_NMatchY[Yj] > 1) {
      m_match_type.push_back(2);
    } else {
      m_match_type.push_back(CheckEvent(N));
    }
  }

  return ArrayOfGoodCouple;
}

////////////////////////////////////////////////////////////////////////////
int TMust2Physics::CheckEvent(int N) {
  if (m_NMatchDet[N] > m_StripXMultDet[N]
      || m_NMatchDet[N] > m_StripYMultDet[N]) {
    // Bad event
    return 2;
  } else {
    // Good event
    return 1;
  }
}

////////////////////////////////////////////////////////////////////////////
bool TMust2Physics::IsValidChannel(const int& DetectorType,
                                   const int& telescope, const int& channel) {
  if (DetectorType == 0)
    return *(m_XChannelStatus[telescope - 1].begin() + channel - 1);

  else if (DetectorType == 1)
    return *(m_YChannelStatus[telescope - 1].begin() + channel - 1);

  else if (DetectorType == 2)
    return *(m_SiLiChannelStatus[telescope - 1].begin() + channel - 1);

  if (DetectorType == 3)
    return *(m_CsIChannelStatus[telescope - 1].begin() + channel - 1);

  else
    return false;
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::ReadAnalysisConfig() {
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

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig
      = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigMust2.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
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

      else if (whatToDo == "MAX_STRIP_MULTIPLICITY") {
        AnalysisConfigFile >> DataBuffer;
        m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str());
        cout << "MAXIMUN STRIP MULTIPLICITY "
             << m_MaximumStripMultiplicityAllowed << endl;
      }

      else if (whatToDo == "STRIP_ENERGY_MATCHING_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingSigma = atof(DataBuffer.c_str());
        cout << "STRIP ENERGY MATCHING SIGMA " << m_StripEnergyMatchingSigma
             << endl;
      }

      else if (whatToDo == "STRIP_ENERGY_MATCHING_NUMBER_OF_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingNumberOfSigma = atoi(DataBuffer.c_str());
        cout << "STRIP ENERGY MATCHING NUMBER OF SIGMA "
             << m_StripEnergyMatchingNumberOfSigma << endl;
      }

      else if (whatToDo == "DISABLE_ALL") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer << endl;
        int          telescope = atoi(DataBuffer.substr(2, 1).c_str());
        vector<bool> ChannelStatus;
        ChannelStatus.resize(128, false);
        m_XChannelStatus[telescope - 1] = ChannelStatus;
        m_YChannelStatus[telescope - 1] = ChannelStatus;
        ChannelStatus.resize(16, false);
        m_SiLiChannelStatus[telescope - 1] = ChannelStatus;
        m_CsIChannelStatus[telescope - 1]  = ChannelStatus;
      }

      else if (whatToDo == "DISABLE_CHANNEL") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer << endl;
        int telescope = atoi(DataBuffer.substr(2, 1).c_str());
        int channel   = -1;
        if (DataBuffer.compare(3, 4, "STRX") == 0) {
          channel = atoi(DataBuffer.substr(7).c_str());
          *(m_XChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else if (DataBuffer.compare(3, 4, "STRY") == 0) {
          channel = atoi(DataBuffer.substr(7).c_str());
          *(m_YChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else if (DataBuffer.compare(3, 4, "SILI") == 0) {
          channel = atoi(DataBuffer.substr(7).c_str());
          *(m_SiLiChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else if (DataBuffer.compare(3, 3, "CSI") == 0) {
          channel = atoi(DataBuffer.substr(6).c_str());
          *(m_CsIChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else
          cout << "Warning: detector type for Must2 unknown!" << endl;

      }

      else if (whatToDo == "TAKE_E_Y") {
        m_Take_E_Y = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "TAKE_T_Y") {
        m_Take_T_Y = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "TAKE_E_X") {
        m_Take_E_Y = false;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "TAKE_T_X") {
        m_Take_T_Y = false;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "IGNORE_NOT_MATCHING_CSI") {
        m_Ignore_not_matching_CsI = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "CSI_SIZE") {
        AnalysisConfigFile >> DataBuffer;
        m_CsI_Size = atoi(DataBuffer.c_str());
        cout << whatToDo << " " << m_CsI_Size << endl;
      }

      else if (whatToDo == "IGNORE_NOT_MATCHING_SILI") {
        m_Ignore_not_matching_SiLi = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "SILI_SIZE") {
        AnalysisConfigFile >> DataBuffer;
        m_SiLi_Size = atoi(DataBuffer.c_str());
        cout << whatToDo << " " << m_SiLi_Size << endl;
      }

      else if (whatToDo == "SI_X_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_X_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_X_E_RAW_Threshold << endl;
      }

      else if (whatToDo == "SI_Y_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_Y_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_Y_E_RAW_Threshold << endl;
      }

      else if (whatToDo == "SILI_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_SiLi_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_SiLi_E_RAW_Threshold << endl;
      }

      else if (whatToDo == "CSI_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_CsI_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_CsI_E_Threshold << endl;
      }

      else if (whatToDo == "SI_X_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_X_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_X_E_Threshold << endl;
      }

      else if (whatToDo == "SI_Y_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_Y_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_Y_E_Threshold << endl;
      }

      else if (whatToDo == "SILI_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_SiLi_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_SiLi_E_Threshold << endl;
      }

      else if (whatToDo == "CSI_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_CsI_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_CsI_E_Threshold << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
bool TMust2Physics::Match_Si_SiLi(int X, int Y, int PadNbr) {

  // remove the central part and surrounding
  if (X < 8 || X > 120 || (Y < 68 && Y > 60))
    return false;

  if (abs(m_SiLi_MatchingX[PadNbr - 1] - X) < m_SiLi_Size / 2.
      && abs(m_SiLi_MatchingY[PadNbr - 1] - Y) < m_SiLi_Size / 2.)
    return true;

  else
    return false;
}

///////////////////////////////////////////////////////////////////////////
bool TMust2Physics::Match_Si_CsI(int X, int Y, int CristalNbr,
                                 int TelescopeNumber) {

  // FIXME: Temporary fix... for zero degree
  // Added to correct visible gaps at zero degree. This seem to be a geometrical
  // issue, this is not required in most cases (ZeroDegreeTelescopeNumber should
  // be added to TMust2Physics.h or to the input file if more experiment use a
  // telescope at zero degree
  // if (TelescopeNumber == ZeroDegreeTelescopeNumber) {
  //   if ((X >= m_ZeroDegree_CsI_MatchingX[CristalNbr - 1].first - 1
  //        && X <= m_ZeroDegree_CsI_MatchingX[CristalNbr - 1].second + 1)
  //       && (Y >= m_ZeroDegree_CsI_MatchingY[CristalNbr - 1].first - 1
  //           && Y <= m_ZeroDegree_CsI_MatchingY[CristalNbr - 1].second + 1)) {
  //     return true;
  //   } else {
  //     return false;
  //   }
  // } else {
  if (abs(m_CsI_MatchingX[CristalNbr - 1] - X) <= (double)m_CsI_Size / 2.
      && abs(m_CsI_MatchingY[CristalNbr - 1] - Y) <= (double)m_CsI_Size / 2.) {
    return true;
  } else {
    return false;
  }
  // }
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::Clear() {
  EventMultiplicity = 0;

  m_match_type.clear();

  m_StripXMultDet.clear();
  m_StripYMultDet.clear();
  m_NMatchDet.clear();
  m_NMatchX.clear();
  m_NMatchY.clear();

  TelescopeNumber.clear();
  EventType.clear();
  TotalEnergy.clear();

  // Si X
  Si_E.clear();
  Si_T.clear();
  Si_X.clear();
  Si_Y.clear();

  // Si(Li)
  SiLi_E.clear();
  SiLi_T.clear();
  SiLi_N.clear();

  // CsI
  CsI_E.clear();
  CsI_T.clear();
  CsI_N.clear();

  Si_EX.clear();
  Si_TX.clear();
  Si_EY.clear();
  Si_TY.clear();
  TelescopeNumber_X.clear();
  TelescopeNumber_Y.clear();
}
///////////////////////////////////////////////////////////////////////////

void TMust2Physics::ReadCalibrationRun() {
  m_StripXEMult = m_EventData->GetMMStripXEMult();
  m_StripYEMult = m_EventData->GetMMStripYEMult();
  m_StripXTMult = m_EventData->GetMMStripXTMult();
  m_StripYTMult = m_EventData->GetMMStripYTMult();
  m_SiLiEMult   = m_EventData->GetMMSiLiEMult();
  m_SiLiTMult   = m_EventData->GetMMSiLiTMult();
  m_CsIEMult    = m_EventData->GetMMCsIEMult();
  m_CsITMult    = m_EventData->GetMMCsITMult();

  //   X
  //   E
  for (unsigned int i = 0; i < m_StripXEMult; ++i) {
    TelescopeNumber_X.push_back(m_EventData->GetMMStripXEDetectorNbr(i));
    Si_EX.push_back(fSi_X_E(m_EventData, i));
    Si_X.push_back(m_EventData->GetMMStripXEStripNbr(i));
  }
  //   T
  for (unsigned int i = 0; i < m_StripXTMult; ++i) {
    TelescopeNumber_X.push_back(m_EventData->GetMMStripXTDetectorNbr(i));
    Si_TX.push_back(fSi_X_T(m_EventData, i));
    Si_X.push_back(m_EventData->GetMMStripXTStripNbr(i));
  }

  //   Y
  //   E
  for (unsigned int i = 0; i < m_StripYEMult; ++i) {
    TelescopeNumber_Y.push_back(m_EventData->GetMMStripYEDetectorNbr(i));
    Si_EY.push_back(fSi_Y_E(m_EventData, i));
    Si_Y.push_back(m_EventData->GetMMStripYEStripNbr(i));
  }

  //   T
  for (unsigned int i = 0; i < m_StripYTMult; ++i) {
    TelescopeNumber_Y.push_back(m_EventData->GetMMStripYTDetectorNbr(i));
    Si_TY.push_back(fSi_Y_T(m_EventData, i));
    Si_Y.push_back(m_EventData->GetMMStripYTStripNbr(i));
  }

  // CsI Energy
  for (unsigned int j = 0; j < m_CsIEMult; ++j) {
    CsI_N.push_back(m_EventData->GetMMCsIECristalNbr(j));
    CsI_E.push_back(fCsI_E(m_EventData, j));
  }

  // CsI Time
  for (unsigned int j = 0; j < m_CsITMult; ++j) {
    CsI_T.push_back(fCsI_T(m_EventData, j));
  }
}

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("M2Telescope");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " Telescope found " << endl;

  // Cartesian Case
  vector<string> cart
      = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128", "SI", "SILI", "CSI"};
  // Spherical Case
  vector<string> sphe = {"R", "THETA", "PHI", "BETA", "SI", "SILI", "CSI"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Must2 Telescope " << i + 1 << endl;
      TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
      TVector3 B = blocks[i]->GetTVector3("X128_Y1", "mm");
      TVector3 C = blocks[i]->GetTVector3("X1_Y128", "mm");
      TVector3 D = blocks[i]->GetTVector3("X128_Y128", "mm");
      AddTelescope(A, B, C, D);
      m_CsIPresent[i + 1]  = blocks[i]->GetInt("CSI");
      m_SiLiPresent[i + 1] = blocks[i]->GetInt("SILI");
    }

    else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Must2 Telescope " << i + 1 << endl;

      double         Theta = blocks[i]->GetDouble("THETA", "deg");
      double         Phi   = blocks[i]->GetDouble("PHI", "deg");
      double         R     = blocks[i]->GetDouble("R", "mm");
      vector<double> beta  = blocks[i]->GetVectorDouble("BETA", "deg");
      cout << Theta << " " << Phi << endl;
      AddTelescope(Theta, Phi, R, beta[0], beta[1], beta[2]);

      m_CsIPresent[i + 1]  = blocks[i]->GetInt("CSI");
      m_SiLiPresent[i + 1] = blocks[i]->GetInt("SILI");
    }

    else {
      cout << "ERROR: Missing token for M2Telescope blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }
  // m_MultEvt = new TMust2MultTelescope[m_NumberOfTelescope];

  InitializeStandardParameter();
  ReadAnalysisConfig();
}
//////////////////////////////////////////////////////////////////////////
void TMust2Physics::InitSpectra() {
  m_Spectra = new TMust2Spectra(m_NumberOfTelescope);
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::FillSpectra() {
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TMust2Physics::CheckSpectra() { m_Spectra->CheckSpectra(); }
///////////////////////////////////////////////////////////////////////////
void TMust2Physics::ClearSpectra() {
  // To be done
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::WriteSpectra() {
  if (m_Spectra)
    m_Spectra->WriteSpectra();
}

///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TMust2Physics::GetSpectra() {
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  // Good for simulation, close to typical values
  vector<double> standardX    = {-63, 63. / 8192.};
  vector<double> standardY    = {63, -63. / 8192.};
  vector<double> standardCsI  = {-250, 250. / 8192.};
  vector<double> standardSiLi = {-250, 250. / 8192.};
  vector<double> standardT    = {-1000, 1000. / 8192.};

  for (int i = 0; i < m_NumberOfTelescope; ++i) {

    for (int j = 0; j < 128; ++j) {
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_E",
          "MUST2_T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_E",
          standardX);
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_E",
          "MUST2_T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_E",
          standardY);
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_T",
          "MUST2_T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_T",
          standardT);
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_T",
          "MUST2_T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_T",
          standardT);
    }

    for (int j = 0; j < 16; ++j) {
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_E",
          "MUST2_T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_E",
          standardSiLi);
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_T",
          "MUST2_T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_T",
          standardT);
    }

    for (int j = 0; j < 16; ++j) {
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E",
          "MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E",
          standardCsI);
      Cal->AddParameter(
          "MUST2", "T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T",
          "MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T",
          standardT);
    }
  }

  return;
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("MUST2", true);
  inputChain->SetBranchStatus("fMM_*", true);
  inputChain->SetBranchAddress("MUST2", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("MUST2", true);
  inputChain->SetBranchStatus("EventMultiplicity", true);
  inputChain->SetBranchStatus("EventType", true);
  inputChain->SetBranchStatus("TelescopeNumber", true);
  inputChain->SetBranchStatus("Si_E", true);
  inputChain->SetBranchStatus("Si_T", true);
  inputChain->SetBranchStatus("Si_X", true);
  inputChain->SetBranchStatus("Si_Y", true);
  inputChain->SetBranchStatus("Si_EX", true);
  inputChain->SetBranchStatus("Si_TX", true);
  inputChain->SetBranchStatus("Si_EY", true);
  inputChain->SetBranchStatus("Si_TY", true);
  inputChain->SetBranchStatus("TelescopeNumber_X", true);
  inputChain->SetBranchStatus("TelescopeNumber_Y", true);
  inputChain->SetBranchStatus("SiLi_E", true);
  inputChain->SetBranchStatus("SiLi_T", true);
  inputChain->SetBranchStatus("SiLi_N", true);
  inputChain->SetBranchStatus("CsI_E", true);
  inputChain->SetBranchStatus("CsI_T", true);
  inputChain->SetBranchStatus("CsI_N", true);
  inputChain->SetBranchStatus("TotalEnergy", true);
  inputChain->SetBranchAddress("MUST2", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("MUST2", "TMust2Physics", &m_EventPhysics);
}

/////   Specific to MUST2Array   ////

void TMust2Physics::AddTelescope(TVector3 C_X1_Y1, TVector3 C_X128_Y1,
                                 TVector3 C_X1_Y128, TVector3 C_X128_Y128) {
  // To avoid warning
  C_X128_Y128 *= 1;

  m_NumberOfTelescope++;

  //   Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that
  //   Y strip are allong X axis)
  TVector3 U      = C_X128_Y1 - C_X1_Y1;
  double   Ushift = (U.Mag() - 98) / 2.;
  U               = U.Unit();
  //   Vector V on Telescope Face (parallele to X Strip)
  TVector3 V      = C_X1_Y128 - C_X1_Y1;
  double   Vshift = (V.Mag() - 98) / 2.;
  V               = V.Unit();

  //   Position Vector of Strip Center
  TVector3 StripCenter = TVector3(0, 0, 0);
  //   Position Vector of X=1 Y=1 Strip
  TVector3 Strip_1_1;

  //   Geometry Parameter
  double Face          = 98; // mm
  double NumberOfStrip = 128;
  double StripPitch    = Face / NumberOfStrip; // mm
  //   Buffer object to fill Position Array
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneTelescopeStripPositionX;
  vector<vector<double>> OneTelescopeStripPositionY;
  vector<vector<double>> OneTelescopeStripPositionZ;

  //   Moving StripCenter to 1.1 corner:
  Strip_1_1 = C_X1_Y1 + (U + V) * (StripPitch / 2.);
  Strip_1_1 += U * Ushift + V * Vshift;

  for (int i = 0; i < 128; ++i) {
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < 128; ++j) {
      StripCenter = Strip_1_1 + StripPitch * (i * U + j * V);
      lineX.push_back(StripCenter.X());
      lineY.push_back(StripCenter.Y());
      lineZ.push_back(StripCenter.Z());
    }

    OneTelescopeStripPositionX.push_back(lineX);
    OneTelescopeStripPositionY.push_back(lineY);
    OneTelescopeStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneTelescopeStripPositionX);
  m_StripPositionY.push_back(OneTelescopeStripPositionY);
  m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}

void TMust2Physics::InitializeStandardParameter() {
  //   Enable all channel
  vector<bool> ChannelStatus;
  m_XChannelStatus.clear();
  m_YChannelStatus.clear();
  m_SiLiChannelStatus.clear();
  m_CsIChannelStatus.clear();

  ChannelStatus.resize(128, true);
  for (int i = 0; i < m_NumberOfTelescope; ++i) {
    m_XChannelStatus[i] = ChannelStatus;
    m_YChannelStatus[i] = ChannelStatus;
  }

  ChannelStatus.resize(16, true);
  for (int i = 0; i < m_NumberOfTelescope; ++i) {
    m_SiLiChannelStatus[i] = ChannelStatus;
    m_CsIChannelStatus[i]  = ChannelStatus;
  }

  m_MaximumStripMultiplicityAllowed = m_NumberOfTelescope;

  return;
}

void TMust2Physics::AddTelescope(double theta, double phi, double distance,
                                 double beta_u, double beta_v, double beta_w) {

  m_NumberOfTelescope++;

  double Pi = 3.141592654;

  cout << theta << " " << phi << endl;
  // convert from degree to radian:
  // theta = theta * Pi / 180.;
  // phi   = phi * Pi / 180.;

  // Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that Y
  // strip are allong X axis)
  TVector3 U;
  // Vector V on Telescope Face (parallele to X Strip)
  TVector3 V;
  // Vector W normal to Telescope Face (pointing CsI)
  TVector3 W;
  // Vector position of Telescope Face center
  TVector3 C;

  C = TVector3(distance * sin(theta) * cos(phi),
               distance * sin(theta) * sin(phi), distance * cos(theta));

  //   std::cout << C.Theta() * 180. / Pi << " " << C.Phi() * 180. / Pi << " "
  //             << C.Mag() << std::endl;

  std::cout << C.X() << " " << C.Y() << std::endl;

  TVector3 P
      = TVector3(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
  // = TVector3(1, 1, 1);

  W = C.Unit();
  // U = W.Cross(P);
  V = W.Cross(P);
  // U = TVector3(0, -1, 0);
  // V = W.Cross(U);
  U = W.Cross(V);
  // V = TVector3(-1, 0, 0);

  U = U.Unit();
  // U.Rotate(90. * Pi / 180.,W);
  V = V.Unit();
  // V.Rotate(90. * Pi / 180.,W);

  std::cout << "U: " << U.X() << " " << U.Y() << " " << U.Z() << " "
            << "V: " << V.X() << " " << V.Y() << " " << V.Z() << " "
            << "W: " << W.X() << " " << W.Y() << " " << W.Z() << std::endl;

  // exit(1);

  U.Rotate(beta_u * Pi / 180., U);
  V.Rotate(beta_u * Pi / 180., U);

  U.Rotate(beta_v * Pi / 180., V);
  V.Rotate(beta_v * Pi / 180., V);

  U.Rotate(beta_w * Pi / 180., W);
  V.Rotate(beta_w * Pi / 180., W);

  double Face          = 98; // mm
  double NumberOfStrip = 128;
  double StripPitch    = Face / NumberOfStrip; // mm

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneTelescopeStripPositionX;
  vector<vector<double>> OneTelescopeStripPositionY;
  vector<vector<double>> OneTelescopeStripPositionZ;

  double X, Y, Z;

  // Moving C to the 1.1 corner:
  C.SetX(C.X() - (Face / 2 - StripPitch / 2) * (V.X() + U.X()));
  C.SetY(C.Y() - (Face / 2 - StripPitch / 2) * (V.Y() + U.Y()));
  C.SetZ(C.Z() - (Face / 2 - StripPitch / 2) * (V.Z() + U.Z()));

  for (int i = 0; i < 128; ++i) {

    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < 128; ++j) {
      X = C.X() + StripPitch * (U.X() * i + V.X() * j);
      Y = C.Y() + StripPitch * (U.Y() * i + V.Y() * j);
      Z = C.Z() + StripPitch * (U.Z() * i + V.Z() * j);

      lineX.push_back(X);
      lineY.push_back(Y);
      lineZ.push_back(Z);
    }

    OneTelescopeStripPositionX.push_back(lineX);
    OneTelescopeStripPositionY.push_back(lineY);
    OneTelescopeStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneTelescopeStripPositionX);
  m_StripPositionY.push_back(OneTelescopeStripPositionY);
  m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}

TVector3 TMust2Physics::GetPositionOfInteraction(const int i) const {
  TVector3 Position
      = TVector3(GetStripPositionX(TelescopeNumber[i], Si_X[i], Si_Y[i]),
                 GetStripPositionY(TelescopeNumber[i], Si_X[i], Si_Y[i]),
                 GetStripPositionZ(TelescopeNumber[i], Si_X[i], Si_Y[i]));

  return (Position);
}

TVector3 TMust2Physics::GetTelescopeNormal(const int i) const {
  TVector3 U = TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1),
                        GetStripPositionY(TelescopeNumber[i], 128, 1),
                        GetStripPositionZ(TelescopeNumber[i], 128, 1))

               - TVector3(GetStripPositionX(TelescopeNumber[i], 1, 1),
                          GetStripPositionY(TelescopeNumber[i], 1, 1),
                          GetStripPositionZ(TelescopeNumber[i], 1, 1));

  TVector3 V = TVector3(GetStripPositionX(TelescopeNumber[i], 128, 128),
                        GetStripPositionY(TelescopeNumber[i], 128, 128),
                        GetStripPositionZ(TelescopeNumber[i], 128, 128))

               - TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1),
                          GetStripPositionY(TelescopeNumber[i], 128, 1),
                          GetStripPositionZ(TelescopeNumber[i], 128, 1));

  TVector3 Normal = U.Cross(V);

  return (Normal.Unit());
}

///////////////////////////////////////////////////////////////////////////
namespace MUST2_LOCAL {
//   DSSD
//   X
double fSi_X_E(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMStripXEDetectorNbr(i));
  name += "_Si_X";
  name += NPL::itoa(m_EventData->GetMMStripXEStripNbr(i));
  name += "_E";
  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMStripXEEnergy(i));
}

double fSi_X_T(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMStripXTDetectorNbr(i));
  name += "_Si_X";
  name += NPL::itoa(m_EventData->GetMMStripXTStripNbr(i));
  name += "_T";
  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMStripXTTime(i));
}

//   Y
double fSi_Y_E(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMStripYEDetectorNbr(i));
  name += "_Si_Y";
  name += NPL::itoa(m_EventData->GetMMStripYEStripNbr(i));
  name += "_E";
  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMStripYEEnergy(i));
}

double fSi_Y_T(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMStripYTDetectorNbr(i));
  name += "_Si_Y";
  name += NPL::itoa(m_EventData->GetMMStripYTStripNbr(i));
  name += "_T";
  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMStripYTTime(i));
}

//   SiLi
double fSiLi_E(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMSiLiEDetectorNbr(i));
  name += "_SiLi";
  name += NPL::itoa(m_EventData->GetMMSiLiEPadNbr(i));
  name += "_E";

  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMSiLiEEnergy(i));
}

double fSiLi_T(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMSiLiTDetectorNbr(i));
  name += "_SiLi";
  name += NPL::itoa(m_EventData->GetMMSiLiTPadNbr(i));
  name += "_T";
  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMSiLiTTime(i));
}

//   CsI
double fCsI_E(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMCsIEDetectorNbr(i));
  name += "_CsI";
  name += NPL::itoa(m_EventData->GetMMCsIECristalNbr(i));
  name += "_E";
  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMCsIEEnergy(i));
}

double fCsI_T(const TMust2Data* m_EventData, const int& i) {
  static string name;
  name = "MUST2/T";
  name += NPL::itoa(m_EventData->GetMMCsITDetectorNbr(i));
  name += "_CsI";
  name += NPL::itoa(m_EventData->GetMMCsITCristalNbr(i));
  name += "_T";
  return CalibrationManager::getInstance()->ApplyCalibration(
      name, m_EventData->GetMMCsITTime(i));
}
} // namespace MUST2_LOCAL

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TMust2Physics::Construct() {
  return (NPL::VDetector*)new TMust2Physics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_must2 {
public:
  proxy_must2() {
    NPL::DetectorFactory::getInstance()->AddToken("M2Telescope", "MUST2");
    NPL::DetectorFactory::getInstance()->AddDetector("M2Telescope",
                                                     TMust2Physics::Construct);
  }
};

proxy_must2 p_must2;
}
