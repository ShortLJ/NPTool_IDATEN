/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : March 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Scone Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSconeData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSconeData)


//////////////////////////////////////////////////////////////////////
TSconeData::TSconeData() {
}



//////////////////////////////////////////////////////////////////////
TSconeData::~TSconeData() {
}



//////////////////////////////////////////////////////////////////////
void TSconeData::Clear() {
  // Energy
  fScone_E_DetectorNbr.clear();
  fScone_E_PlasticNbr.clear();
  fScone_Energy.clear();
  // Time
  fScone_T_DetectorNbr.clear();
  fScone_T_PlasticNbr.clear();
  fScone_Time.clear();
  // Flah for simulation
  fScone_HasCaptured.clear();
  fScone_CaptureTime.clear();
  fScone_GammaEnergy.clear();
  fScone_ProtonEnergy.clear();
  fScone_ProtonTime.clear();
  fScone_FCProcess.clear();
}



//////////////////////////////////////////////////////////////////////
void TSconeData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSconeData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fScone_E_DetectorNbr.size();
  cout << "Scone_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fScone_E_DetectorNbr[i]
         << " Energy: " << fScone_Energy[i];
  }
  
  // Time
  mysize = fScone_T_DetectorNbr.size();
  cout << "Scone_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fScone_T_DetectorNbr[i]
         << " Time: " << fScone_Time[i];
  }
}
