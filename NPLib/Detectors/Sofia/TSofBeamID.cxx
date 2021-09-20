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
 *  This class hold SofBeamID Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TSofBeamID.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSofBeamID)


//////////////////////////////////////////////////////////////////////
TSofBeamID::TSofBeamID() {
}



//////////////////////////////////////////////////////////////////////
TSofBeamID::~TSofBeamID() {
}



//////////////////////////////////////////////////////////////////////
void TSofBeamID::Clear() {
  fBeam_Z     = -1;
  fBeam_Qmax  = -1;
  fBeam_AoQ   = -1;
  fBeam_A     = -1;
  fBeam_Beta  = -1;
  fBeam_Gamma = -1;
  fBeam_Brho  = -1;
  fBeam_XS2   = -1000;
  fBeam_XCC   = -1000;
  fBeam_YCC   = -1000;

}



//////////////////////////////////////////////////////////////////////
void TSofBeamID::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSofBeamID::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  cout << "fBeam_Z: " << fBeam_Z << endl;
  cout << "fBeam_AoQ: " << fBeam_AoQ << endl;
  cout << "fBeam_A: " << fBeam_A << endl;
  cout << "fBeam_Beta: " << fBeam_Beta << endl;
  cout << "fBeam_Gamma: " << fBeam_Gamma << endl;
  cout << "fBeam_Brho: " << fBeam_Brho << endl;
  cout << "fBeam_XS2: " << fBeam_XS2 << endl;
  cout << "fBeam_XCC: " << fBeam_XCC << endl;
  cout << "fBeam_YCC: " << fBeam_YCC << endl;
 
}
