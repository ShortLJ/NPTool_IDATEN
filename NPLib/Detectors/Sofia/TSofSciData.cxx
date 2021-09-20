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
 *  This class hold SofSci Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSofSciData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofSciData)


//////////////////////////////////////////////////////////////////////
TSofSciData::TSofSciData() {
}



//////////////////////////////////////////////////////////////////////
TSofSciData::~TSofSciData() {
}



//////////////////////////////////////////////////////////////////////
void TSofSciData::Clear() {
  fSofSci_DetNbr.clear();
  fSofSci_Pmt.clear();
  fSofSci_CT.clear();
  fSofSci_FT.clear();
}



//////////////////////////////////////////////////////////////////////
void TSofSciData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofSciData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fSofSci_DetNbr.size();
  cout << "SofSci_Mult: " << GetMultiplicity() << endl;
 
}
