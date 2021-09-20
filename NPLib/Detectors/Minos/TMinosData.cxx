/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin                                        *
 * Maintainer : Adrien Matta                                                 *
 * contact address:  matta@lpccaen.in2p3.fr                                  *
 *                                                                           *
 *                                                                           *
 * Creation Date  : October 2018                                             *
 * Last update    : April 2021                                               *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Raw data                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TMinosData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TMinosData)

//////////////////////////////////////////////////////////////////////
TMinosData::TMinosData() {
}

//////////////////////////////////////////////////////////////////////
TMinosData::~TMinosData() {
}

//////////////////////////////////////////////////////////////////////
void TMinosData::Clear() {
  // Minos_Pads
  fMinos_PadNumber.clear();
  fMinos_HitCount.clear();
  fMinos_Charge.clear();
  fMinos_Time.clear();
}
