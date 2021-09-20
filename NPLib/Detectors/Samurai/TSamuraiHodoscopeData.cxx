/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : April 2021                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiHodoscope Raw data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSamuraiHodoscopeData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TSamuraiHodoscopeData)


//////////////////////////////////////////////////////////////////////
TSamuraiHodoscopeData::TSamuraiHodoscopeData() {
}



//////////////////////////////////////////////////////////////////////
TSamuraiHodoscopeData::~TSamuraiHodoscopeData() {
}



//////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopeData::Clear() {
    // UP // 
    // Charge 
    fSamuraiHodoscope_Qu_ID.clear();
    fSamuraiHodoscope_Qu_Charge.clear();
    
    // Time
    fSamuraiHodoscope_Tu_ID.clear();
    fSamuraiHodoscope_Tu_Time.clear();
    
    // DOWN // 
    // Charge 
    fSamuraiHodoscope_Qd_ID.clear();
    fSamuraiHodoscope_Qd_Charge.clear();
    
    // Time
    fSamuraiHodoscope_Td_ID.clear();
    fSamuraiHodoscope_Td_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopeData::Dump() const {
/*  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSamuraiHodoscopeData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fSamuraiHodoscope_E_DetectorNbr.size();
  cout << "SamuraiHodoscope_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fSamuraiHodoscope_E_DetectorNbr[i]
         << " Energy: " << fSamuraiHodoscope_Energy[i];
  }
  
  // Time
  mysize = fSamuraiHodoscope_T_DetectorNbr.size();
  cout << "SamuraiHodoscope_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fSamuraiHodoscope_T_DetectorNbr[i]
         << " Time: " << fSamuraiHodoscope_Time[i];
  }
  */
}
