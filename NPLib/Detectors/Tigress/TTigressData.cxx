/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold the Tigress  raw data (Made for TIG10 card)              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "TTigressData.h"

ClassImp(TTigressData)

/////////////////////////
TTigressData::TTigressData(){
}

/////////////////////////
TTigressData::~TTigressData(){
}

/////////////////////////
void TTigressData::Clear(){
  fTIG_Ge_CloverNbr.clear();
  fTIG_Ge_CrystalNbr.clear();
  fTIG_Ge_SegmentNbr.clear();
  fTIG_Ge_Energy.clear();
  fTIG_Ge_TimeCFD.clear();
  fTIG_Ge_TimeLED.clear();

  fTIG_Ge_CloverNbrAddback.clear();
  fTIG_Ge_EnergyAddback.clear();

  fTIG_BGO_CloverNbr.clear();
  fTIG_BGO_CrystalNbr.clear();
  fTIG_BGO_PmNbr.clear();
  fTIG_BGO_Energy.clear();
  fTIG_BGO_TimeCFD.clear();
  fTIG_BGO_TimeLED.clear();
}

/////////////////////////
void TTigressData::Dump() const{
  // Energy
 // cout << "Tigress_Mult = " << fTIG_CloverNbr.size() << endl;
  
  // Front
 // for (UShort_t i = 0; i < fTIG_CloverNbr.size(); i++){
 //   cout << "Clover: " << fTIG_CloverNbr[i]
 //        << " Crystal: " << fTIG_CrystalNbr[i]
 //        << " Energy: " << fTIG_Energy[i]
 //        << " Time: " << fTIG_Time[i] << endl;
 // }
}

void TTigressData::Addback(){
	fTIG_Ge_CloverNbrAddback.clear();
	fTIG_Ge_EnergyAddback.clear();

	if(fTIG_Ge_CloverNbr.size() == 0) return;

	for(unsigned int i=0; i<fTIG_Ge_CloverNbr.size(); i++){
		bool iflag=true;
		for(unsigned int j=0; j<fTIG_Ge_CloverNbrAddback.size(); j++){
			if(fTIG_Ge_CloverNbr[i]==fTIG_Ge_CloverNbrAddback[j]){
				fTIG_Ge_EnergyAddback[j]+=fTIG_Ge_Energy[i];
				iflag=false; break;
			}
		}
		if(iflag){
			fTIG_Ge_CloverNbrAddback.push_back(fTIG_Ge_CloverNbr[i]);
			fTIG_Ge_EnergyAddback.push_back(fTIG_Ge_Energy[i]);
		}
	}

}












