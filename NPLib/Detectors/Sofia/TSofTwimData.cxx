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
 *  This class hold SofTwim Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSofTwimData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofTwimData)


//////////////////////////////////////////////////////////////////////
TSofTwimData::TSofTwimData() {
}



//////////////////////////////////////////////////////////////////////
TSofTwimData::~TSofTwimData() {
}



//////////////////////////////////////////////////////////////////////
void TSofTwimData::Clear() {
  fTwim_SectionNbr.clear();
  fTwim_AnodeNbr.clear();
  fTwim_Energy.clear();
  fTwim_DriftTime.clear();
  fTwim_PileUp.clear();
  fTwim_Overflow.clear();
}



//////////////////////////////////////////////////////////////////////
void TSofTwimData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofTwimData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fTwim_AnodeNbr.size();
  cout << "Twim_Mult: " << GetMultiplicity() << endl;
 
}
