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
 *  This class hold FissionFragment Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSofFissionFragment.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofFissionFragment)


//////////////////////////////////////////////////////////////////////
TSofFissionFragment::TSofFissionFragment() {
}



//////////////////////////////////////////////////////////////////////
TSofFissionFragment::~TSofFissionFragment() {
}



//////////////////////////////////////////////////////////////////////
void TSofFissionFragment::Clear() {
  fFF_Z.clear();
  fFF_Qmax.clear();
  fFF_AoQ.clear();
  fFF_A.clear();
  fFF_Beta.clear();
  fFF_TOF.clear();
  fFF_Gamma.clear();
  fFF_Brho.clear();
  
  fFF_Zsum = -1;
}



//////////////////////////////////////////////////////////////////////
void TSofFissionFragment::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofFissionFragment::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  for(int i=0; i<fFF_Z.size(); i++){
    cout << "fFF_Z: " << fFF_Z[i] << endl;
    cout << "fFF_AoQ: " << fFF_AoQ[i] << endl;
    cout << "fFF_A: " << fFF_A[i] << endl;
    cout << "fFF_Beta: " << fFF_Beta[i] << endl;
    cout << "fFF_Gamma: " << fFF_Gamma[i] << endl;
    cout << "fFF_Brho: " << fFF_Brho[i] << endl;
  } 
}
