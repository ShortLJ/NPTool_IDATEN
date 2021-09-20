/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Nebula Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TNebulaData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TNebulaData)


//////////////////////////////////////////////////////////////////////
TNebulaData::TNebulaData() {
}



//////////////////////////////////////////////////////////////////////
TNebulaData::~TNebulaData() {
}



//////////////////////////////////////////////////////////////////////
void TNebulaData::Clear() {
    // UP // 
    // Charge 
    fNebula_Qu_ID.clear();
    fNebula_Qu_Charge.clear();
    
    // Time
    fNebula_Tu_ID.clear();
    fNebula_Tu_Time.clear();
    
    // DOWN // 
    // Charge 
    fNebula_Qd_ID.clear();
    fNebula_Qd_Charge.clear();
    
    // Time
    fNebula_Td_ID.clear();
    fNebula_Td_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TNebulaData::Dump() const {
/*  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TNebulaData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fNebula_E_DetectorNbr.size();
  cout << "Nebula_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fNebula_E_DetectorNbr[i]
         << " Energy: " << fNebula_Energy[i];
  }
  
  // Time
  mysize = fNebula_T_DetectorNbr.size();
  cout << "Nebula_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fNebula_T_DetectorNbr[i]
         << " Time: " << fNebula_Time[i];
  }
  */
}
