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
 *  This class hold SofMwpc Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSofMwpcData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofMwpcData)


//////////////////////////////////////////////////////////////////////
TSofMwpcData::TSofMwpcData() {
}



//////////////////////////////////////////////////////////////////////
TSofMwpcData::~TSofMwpcData() {
}



//////////////////////////////////////////////////////////////////////
void TSofMwpcData::Clear() {
  fMwpc_DetNbr.clear();
  fMwpc_Plane.clear();
  fMwpc_Pad.clear();
  fMwpc_Charge.clear();
}



//////////////////////////////////////////////////////////////////////
void TSofMwpcData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofMwpcData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fMwpc_Charge.size();
  cout << "MWPC_Mult: " << GetMultiplicity() << endl;
 
}
