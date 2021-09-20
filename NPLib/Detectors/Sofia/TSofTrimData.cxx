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
 *  This class hold SofTrim Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSofTrimData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofTrimData)


//////////////////////////////////////////////////////////////////////
TSofTrimData::TSofTrimData() {
}



//////////////////////////////////////////////////////////////////////
TSofTrimData::~TSofTrimData() {
}



//////////////////////////////////////////////////////////////////////
void TSofTrimData::Clear() {
  fTrim_SectionNbr.clear();
  fTrim_AnodeNbr.clear();
  fTrim_Energy.clear();
  fTrim_DriftTime.clear();
  fTrim_PileUp.clear();
  fTrim_Overflow.clear();
}



//////////////////////////////////////////////////////////////////////
void TSofTrimData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofTrimData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fTrim_AnodeNbr.size();
  cout << "Trim_Mult: " << GetMultiplicity() << endl;
 
}
