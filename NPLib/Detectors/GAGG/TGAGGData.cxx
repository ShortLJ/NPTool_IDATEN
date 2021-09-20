/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Valerian Alcindor  contact address:                         *
 *                                                                           *
 * Creation Date  : October 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold GAGG Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TGAGGData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TGAGGData)


//////////////////////////////////////////////////////////////////////
TGAGGData::TGAGGData() {
}



//////////////////////////////////////////////////////////////////////
TGAGGData::~TGAGGData() {
}



//////////////////////////////////////////////////////////////////////
void TGAGGData::Clear() {
  // Energy
  fGAGG_E_DetectorNbr.clear();
  fGAGG_Energy.clear();
  // Time
  fGAGG_T_DetectorNbr.clear();
  fGAGG_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TGAGGData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TGAGGData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fGAGG_E_DetectorNbr.size();
  cout << "GAGG_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fGAGG_E_DetectorNbr[i]
         << " Energy: " << fGAGG_Energy[i];
  }
  
  // Time
  mysize = fGAGG_T_DetectorNbr.size();
  cout << "GAGG_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fGAGG_T_DetectorNbr[i]
         << " Time: " << fGAGG_Time[i];
  }
}
