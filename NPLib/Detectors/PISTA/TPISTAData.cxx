/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TPISTAData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TPISTAData)


//////////////////////////////////////////////////////////////////////
TPISTAData::TPISTAData() {
}



//////////////////////////////////////////////////////////////////////
TPISTAData::~TPISTAData() {
}



//////////////////////////////////////////////////////////////////////
void TPISTAData::Clear() {
  // Energy X
  fFirstStage_XE_DetectorNbr.clear();
  fFirstStage_XE_StripNbr.clear();
  fFirstStage_XE_Energy.clear();
  // Energy Y
  fFirstStage_YE_DetectorNbr.clear();
  fFirstStage_YE_StripNbr.clear();
  fFirstStage_YE_Energy.clear();
  // Time X
  fFirstStage_XT_DetectorNbr.clear();
  fFirstStage_XT_StripNbr.clear();
  fFirstStage_XT_Time.clear();
  // Time Y
  fFirstStage_YT_DetectorNbr.clear();
  fFirstStage_YT_StripNbr.clear();
  fFirstStage_YT_Time.clear();

  // Energy X
  fSecondStage_XE_DetectorNbr.clear();
  fSecondStage_XE_StripNbr.clear();
  fSecondStage_XE_Energy.clear();
  // Energy Y
  fSecondStage_YE_DetectorNbr.clear();
  fSecondStage_YE_StripNbr.clear();
  fSecondStage_YE_Energy.clear();
  // Time X
  fSecondStage_XT_DetectorNbr.clear();
  fSecondStage_XT_StripNbr.clear();
  fSecondStage_XT_Time.clear();
  // Time Y
  fSecondStage_YT_DetectorNbr.clear();
  fSecondStage_YT_StripNbr.clear();
  fSecondStage_YT_Time.clear();

}



//////////////////////////////////////////////////////////////////////
void TPISTAData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TPISTAData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fFirstStage_XE_DetectorNbr.size();
  cout << "First Stage PISTA_XE_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "X-DetNbr: " << fFirstStage_XE_DetectorNbr[i]
         << " X-Energy: " << fFirstStage_XE_Energy[i];
  }
  
  // Time
  mysize = fFirstStage_XT_DetectorNbr.size();
  cout << "First Stage PISTA_XT_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "X-DetNbr: " << fFirstStage_XT_DetectorNbr[i]
         << " X-Time: " << fFirstStage_XT_Time[i];
  }
}
