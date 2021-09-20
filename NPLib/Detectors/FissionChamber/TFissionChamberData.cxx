/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : September 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FissionChamber Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TFissionChamberData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TFissionChamberData)


//////////////////////////////////////////////////////////////////////
TFissionChamberData::TFissionChamberData() {
}



//////////////////////////////////////////////////////////////////////
TFissionChamberData::~TFissionChamberData() {
}



//////////////////////////////////////////////////////////////////////
void TFissionChamberData::Clear() {
  // Energy
  fFissionChamber_E_DetectorNbr.clear();
  fFissionChamber_Energy.clear();
  // Time
  fFissionChamber_T_DetectorNbr.clear();
  fFissionChamber_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TFissionChamberData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TFissionChamberData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fFissionChamber_E_DetectorNbr.size();
  cout << "FissionChamber_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fFissionChamber_E_DetectorNbr[i]
         << " Energy: " << fFissionChamber_Energy[i];
  }
  
  // Time
  mysize = fFissionChamber_T_DetectorNbr.size();
  cout << "FissionChamber_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fFissionChamber_T_DetectorNbr[i]
         << " Time: " << fFissionChamber_Time[i];
  }
}
