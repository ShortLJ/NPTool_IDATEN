/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TPISTASpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TPISTASpectra::TPISTASpectra() 
   : fNumberOfDetectors(0) {
  SetName("PISTA");
}



////////////////////////////////////////////////////////////////////////////////
TPISTASpectra::TPISTASpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TPISTASpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("PISTA");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TPISTASpectra::~TPISTASpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TPISTASpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "PISTA"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "PISTA/RAW");
    // Time 
    name = "PISTA"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "PISTA/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TPISTASpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "PISTA"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "PISTA/CAL");
    // Time
    name = "PISTA"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "PISTA/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TPISTASpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "PISTA_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "PISTA/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TPISTASpectra::FillRawSpectra(TPISTAData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetFirstStageMultXEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "PISTA"+NPL::itoa(RawData->GetFirstStage_XE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "PISTA/RAW";

    FillSpectra(family,name,RawData->GetFirstStage_XE_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetFirstStageMultXTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "PISTA"+NPL::itoa(RawData->GetFirstStage_XT_DetectorNbr(i))+"_TIME_RAW";
    family = "PISTA/RAW";

    FillSpectra(family,name,RawData->GetFirstStage_XT_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TPISTASpectra::FillPreTreatedSpectra(TPISTAData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetFirstStageMultXEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "PISTA"+NPL::itoa(PreTreatedData->GetFirstStage_XE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "PISTA/CAL";

    FillSpectra(family,name,PreTreatedData->GetFirstStage_XE_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetFirstStageMultXTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "PISTA"+NPL::itoa(PreTreatedData->GetFirstStage_XT_DetectorNbr(i))+"_TIME_CAL";
    family = "PISTA/CAL";

    FillSpectra(family,name,PreTreatedData->GetFirstStage_XT_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TPISTASpectra::FillPhysicsSpectra(TPISTAPhysics* Physics) {
  static string name;
  static string family;
  family= "PISTA/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->E.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "PISTA_ENERGY_TIME";
    FillSpectra(family,name,Physics->E[i],Physics->Time[i]);
  }
}

