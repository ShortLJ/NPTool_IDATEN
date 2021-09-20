/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny    contact : flavigny@lpccaen.in2p3.fr       *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Strasse Raw data                                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TStrasseData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TStrasseData)


//////////////////////////////////////////////////////////////////////
TStrasseData::TStrasseData() {
}



//////////////////////////////////////////////////////////////////////
TStrasseData::~TStrasseData() {
}



//////////////////////////////////////////////////////////////////////
void TStrasseData::Clear() {
  // Energy X
  fInner_TE_DetectorNbr.clear();
  fInner_TE_StripNbr.clear();
  fInner_TE_Energy.clear();
  // Energy L
  fInner_LE_DetectorNbr.clear();
  fInner_LE_StripNbr.clear();
  fInner_LE_Energy.clear();

  // Energy X
  fOuter_TE_DetectorNbr.clear();
  fOuter_TE_StripNbr.clear();
  fOuter_TE_Energy.clear();
  // Energy L
  fOuter_LE_DetectorNbr.clear();
  fOuter_LE_StripNbr.clear();
  fOuter_LE_Energy.clear();

}



//////////////////////////////////////////////////////////////////////
void TStrasseData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TStrasseData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fInner_TE_DetectorNbr.size();
  cout << "Inner Strasse_TE_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "T-DetNbr: " << fInner_TE_DetectorNbr[i]
         << " T-Energy: " << fInner_TE_Energy[i];
  }
  
}
