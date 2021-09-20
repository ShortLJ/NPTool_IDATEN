/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: M. Labiche    contact address: marc.labiche@stfc.ac.uk   *
 *                                                                           *
 * Creation Date  : 04/12/2009                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class described the raw data of the Khala detector                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
#include "TKhalaData.h"
ClassImp(TKhalaData)

////////////////////////////////////////////////////////////////////////////////
TKhalaData::TKhalaData(){
}

TKhalaData::~TKhalaData(){
}

////////////////////////////////////////////////////////////////////////////////
void TKhalaData::Clear(){
   fKHALA_LaBr3_E_DetectorNbr.clear();
   fKHALA_LaBr3_E_Energy.clear();
   // Time
   fKHALA_LaBr3_T_DetectorNbr.clear();
   fKHALA_LaBr3_T_Time.clear();
}

////////////////////////////////////////////////////////////////////////////////
void TKhalaData::Dump() const{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event XXXXXXXXXXXXXXXXX" << endl;
   // E
   for (UShort_t i = 0; i < fKHALA_LaBr3_E_DetectorNbr.size(); i++)
      cout << "DetNbr: " << fKHALA_LaBr3_E_DetectorNbr[i] << " Energy: " << fKHALA_LaBr3_E_Energy[i] << endl;
   
   // T
   for (UShort_t i = 0; i < fKHALA_LaBr3_T_DetectorNbr.size(); i++)
      cout << "DetNbr: " << fKHALA_LaBr3_T_DetectorNbr[i] << " Time: " << fKHALA_LaBr3_T_Time[i] << endl;
}
