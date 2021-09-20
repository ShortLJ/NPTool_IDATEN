/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold eAGanil Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TeAGanilData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TeAGanilData)


//////////////////////////////////////////////////////////////////////
TeAGanilData::TeAGanilData() {
}



//////////////////////////////////////////////////////////////////////
TeAGanilData::~TeAGanilData() {
}



//////////////////////////////////////////////////////////////////////
void TeAGanilData::Clear() {
  // Energy
  feAGanil_E_DetectorNbr.clear();
  feAGanil_Energy.clear();
  // Time
  feAGanil_T_DetectorNbr.clear();
  feAGanil_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TeAGanilData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TeAGanilData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = feAGanil_E_DetectorNbr.size();
  cout << "eAGanil_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << feAGanil_E_DetectorNbr[i]
         << " Energy: " << feAGanil_Energy[i];
  }
  
  // Time
  mysize = feAGanil_T_DetectorNbr.size();
  cout << "eAGanil_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << feAGanil_T_DetectorNbr[i]
         << " Time: " << feAGanil_Time[i];
  }
}
