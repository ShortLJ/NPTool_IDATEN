/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2019                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ChiNu Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TChiNuData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TChiNuData)


//////////////////////////////////////////////////////////////////////
TChiNuData::TChiNuData() {
}



//////////////////////////////////////////////////////////////////////
TChiNuData::~TChiNuData() {
}



//////////////////////////////////////////////////////////////////////
void TChiNuData::Clear() {
  // Energy
  fChiNu_E_DetectorNbr.clear();
  fChiNu_Energy.clear();
  // Time
  fChiNu_T_DetectorNbr.clear();
  fChiNu_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TChiNuData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TChiNuData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fChiNu_E_DetectorNbr.size();
  cout << "ChiNu_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fChiNu_E_DetectorNbr[i]
         << " Energy: " << fChiNu_Energy[i];
  }
  
  // Time
  mysize = fChiNu_T_DetectorNbr.size();
  cout << "ChiNu_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fChiNu_T_DetectorNbr[i]
         << " Time: " << fChiNu_Time[i];
  }
}
