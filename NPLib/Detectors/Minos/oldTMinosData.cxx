/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: tronchin@lpccaen.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Raw data                                    *
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

/* TGraph* TMinosData::GetChargeAsGraph(const unsigned int& i){ */
/*   TGraph* res = new TGraph (fMinos_Charge[i].size(),&fMinos_Charge[i],&fMinos_Time[i]); */
/*   /1* res->Sort(); *1/ */
/*   return res; */
/* } */

//////////////////////////////////////////////////////////////////////
void TMinosData::Clear() {

  // Minos_Pads
  fMinos_PadNumber.clear();
  fMinos_PadX.clear();
  fMinos_PadY.clear();
  /* fMinos_DriftTime.clear(); */
  /* fMinos_PadCharge.clear(); */ 
  fMinos_Charge.clear();
  fMinos_Time.clear();

  MINOSx_0.clear();
  MINOSy_0.clear();
  MINOSz_0.clear();
  MINOS_D_min.clear();
  MINOS_Radius.clear();
  MINOS_NumberTracks.clear();
  theta0.clear();
  phi0.clear();
  energy0.clear();
}

//////////////////////////////////////////////////////////////////////
/* void TMinosData::Dump() const { */
/*   // This method is very useful for debuging and worth the dev. */
/*   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TMinosData::Dump()] XXXXXXXXXXXXXXXXX" << endl; */

/*   // Energy */
/*   size_t mysize = fMinos_E_DetectorNbr.size(); */
/*   cout << "Minos_E_Mult: " << mysize << endl; */
 
/*   for (size_t i = 0 ; i < mysize ; i++){ */
/*     cout << "DetNbr: " << fMinos_E_DetectorNbr[i] */
/*          << " Energy: " << fMinos_Energy[i]; */
/*   } */
 
