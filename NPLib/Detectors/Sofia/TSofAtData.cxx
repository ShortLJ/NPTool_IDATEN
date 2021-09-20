/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 * Creation Date  : May 2021                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SofAt Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSofAtData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofAtData)


//////////////////////////////////////////////////////////////////////
TSofAtData::TSofAtData() {
}



//////////////////////////////////////////////////////////////////////
TSofAtData::~TSofAtData() {
}



//////////////////////////////////////////////////////////////////////
void TSofAtData::Clear() {
  fAT_AnodeNbr.clear();
  fAT_Energy.clear();
  fAT_Time.clear();
  fAT_PileUp.clear();
  fAT_Overflow.clear();
}



//////////////////////////////////////////////////////////////////////
void TSofAtData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofAtData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fAT_AnodeNbr.size();
  cout << "AT_Mult: " << GetMultiplicity() << endl;
 
}
