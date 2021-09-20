/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny    contact : flavigny@lpccaen.in2p3.fr       *
 *                                                                           *
 * Creation Date  : May 2020                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Strasse Spectra                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TStrasseSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TStrasseSpectra::TStrasseSpectra() 
   : fNumberOfDetectors(0) {
  SetName("Strasse");
}



////////////////////////////////////////////////////////////////////////////////
TStrasseSpectra::TStrasseSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TStrasseSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("Strasse");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TStrasseSpectra::~TStrasseSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TStrasseSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Strasse"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Strasse/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TStrasseSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Strasse"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Strasse/CAL");
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TStrasseSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "Strasse_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "Strasse/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TStrasseSpectra::FillRawSpectra(TStrasseData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetInnerMultTEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Strasse"+NPL::itoa(RawData->GetInner_TE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "Strasse/RAW";

    FillSpectra(family,name,RawData->GetInner_TE_Energy(i));
  }

}



////////////////////////////////////////////////////////////////////////////////
void TStrasseSpectra::FillPreTreatedSpectra(TStrasseData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetInnerMultTEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Strasse"+NPL::itoa(PreTreatedData->GetInner_TE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "Strasse/CAL";

    FillSpectra(family,name,PreTreatedData->GetInner_TE_Energy(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TStrasseSpectra::FillPhysicsSpectra(TStrassePhysics* Physics) {
  static string name;
  static string family;
  family= "Strasse/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->E.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    //name = "Strasse_ENERGY_TIME";
    //FillSpectra(family,name,Physics->E[i],Physics->Time[i]);
  }
}

