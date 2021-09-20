/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Catana Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TCatanaData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TCatanaData)


//////////////////////////////////////////////////////////////////////
TCatanaData::TCatanaData() {
}



//////////////////////////////////////////////////////////////////////
TCatanaData::~TCatanaData() {
}



//////////////////////////////////////////////////////////////////////
void TCatanaData::Clear() {
  // Energy
  fCatana_E_DetectorNbr.clear();
  fCatana_Energy.clear();
  // Time
  fCatana_T_DetectorNbr.clear();
  fCatana_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TCatanaData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TCatanaData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fCatana_E_DetectorNbr.size();
  cout << "Catana_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fCatana_E_DetectorNbr[i]
         << " Energy: " << fCatana_Energy[i];
  }
  
  // Time
  mysize = fCatana_T_DetectorNbr.size();
  cout << "Catana_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fCatana_T_DetectorNbr[i]
         << " Time: " << fCatana_Time[i];
  }
}
