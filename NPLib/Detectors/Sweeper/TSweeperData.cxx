/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: B. Monteagudo  contact address: monteagu@frib.msu.edu                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Sweeper Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSweeperData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSweeperData)


//////////////////////////////////////////////////////////////////////
TSweeperData::TSweeperData() {
}



//////////////////////////////////////////////////////////////////////
TSweeperData::~TSweeperData() {
}



//////////////////////////////////////////////////////////////////////
void TSweeperData::Clear() {
  // Energy
  fSweeper_E_DetectorNbr.clear();
  fSweeper_Energy.clear();
  // Time
  fSweeper_T_DetectorNbr.clear();
  fSweeper_Time.clear();
  //Position
  fSweeper_DC_DetectorNbr.clear();
  fSweeper_X.clear();
  fSweeper_Y.clear();
  fSweeper_DriftTime.clear();
}



//////////////////////////////////////////////////////////////////////
void TSweeperData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSweeperData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fSweeper_E_DetectorNbr.size();
  cout << "Sweeper_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fSweeper_E_DetectorNbr[i]
         << " Energy: " << fSweeper_Energy[i];
  }
  
  // Time
  mysize = fSweeper_T_DetectorNbr.size();
  cout << "Sweeper_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fSweeper_T_DetectorNbr[i]
         << " Time: " << fSweeper_Time[i];
  }
  // Position
  mysize = fSweeper_DC_DetectorNbr.size();
  cout << "Sweeper_DC_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fSweeper_DC_DetectorNbr[i]
         << " Drift Time: " << fSweeper_DriftTime[i]
         << " X position: " << fSweeper_X[i]
	 << " Y position: " << fSweeper_Y[i];   
  }
  
}
