/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Warren Lynch  contact address: warren.lynch@york.ac.uk                        *
 *                                                                           *
 * Creation Date  : June 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TACTIC Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TTACTICData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TTACTICData)


//////////////////////////////////////////////////////////////////////
TTACTICData::TTACTICData() {
}



//////////////////////////////////////////////////////////////////////
TTACTICData::~TTACTICData() {
}



//////////////////////////////////////////////////////////////////////
void TTACTICData::Clear() {
  // Energy
  fTACTIC_E_DetectorNbr.clear();
  fTACTIC_Energy.clear();
  // Time
  fTACTIC_T_DetectorNbr.clear();
  fTACTIC_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TTACTICData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TTACTICData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fTACTIC_E_DetectorNbr.size();
  cout << "TACTIC_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fTACTIC_E_DetectorNbr[i]
         << " Energy: " << fTACTIC_Energy[i];
  }
  
  // Time
  mysize = fTACTIC_T_DetectorNbr.size();
  cout << "TACTIC_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fTACTIC_T_DetectorNbr[i]
         << " Time: " << fTACTIC_Time[i];
  }
}
