/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                  B. Le Crom                        lecrom@ipno.in2p3.fr   *
 * Creation Date  : April 2014                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for ComptonTelescope      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + first version (not complete yet)                                     *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// NPL
#include "TComptonTelescopeSpectra.h"
#include "NPOptionManager.h"
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

// STL
#include <stdexcept>
#include <iostream>  
#include <cstdlib>
#include <string>
using namespace std;


////////////////////////////////////////////////////////////////////////////////
TComptonTelescopeSpectra::TComptonTelescopeSpectra(){
  SetName("ComptonTelescope");
  fNumberOfTelescope = 0;
  fNumberOfDetectors = 1;
  fNumberOfStripsFront=32;
  fNumberOfStripsBack=32;
  fStripEnergyMatchingSigma = 0.006;
  fStripEnergyMatchingNumberOfSigma = 2;
  fNumberOfCounters=50;
  fCalorimeterNPixels=64;
}



////////////////////////////////////////////////////////////////////////////////
TComptonTelescopeSpectra::TComptonTelescopeSpectra(unsigned int NumberOfTelescope)
{
   if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
      cout << "************************************************" << endl
         << "TComptonTelescopeSpectra : Initialising control spectra for " 
         << NumberOfTelescope << " Telescopes" << endl
         << "************************************************" << endl ;

   SetName("ComptonTelescope");
   fNumberOfTelescope = NumberOfTelescope;
   fNumberOfDetectors = 1;
   fNumberOfStripsFront=32;
   fNumberOfStripsBack=32;
   fStripEnergyMatchingSigma = 0.006;
   fStripEnergyMatchingNumberOfSigma = 2;
   fNumberOfCounters=50;
   fCalorimeterNPixels=64;

   InitRawSpectra();
   InitPreTreatedSpectra();
   InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TComptonTelescopeSpectra::~TComptonTelescopeSpectra()
{
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::InitRawSpectra()
{
  string name;
  int ntot = (fNumberOfStripsFront + fNumberOfStripsBack);

  for (unsigned int i = 0; i < fNumberOfTelescope; i++) { // loop on number of telescope

    // DSSSD
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {

/*     // FRONT_E_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_E_RAW";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 512, 0, 1024, "COMPTONTELESCOPE/RAW/ENERGY");

      // BACK_E_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_E_RAW";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 512, 0, 1024, "COMPTONTELESCOPE/RAW/ENERGY");
*/
      // E Front and Back vs Strip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_BACK_E_RAW";
      AddHisto2D(name, name, ntot, 0, ntot, 512, 0, 1024, "COMPTONTELESCOPE/RAW/ENERGY");
 
      // FRONT_T_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_T_RAW";
      AddHisto1D(name, name, 5000, 0, 5e9, "COMPTONTELESCOPE/RAW/TIME");

      // BACK_T_RAW
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_T_RAW";
      AddHisto1D(name, name, 5000, 0, 5e9, "COMPTONTELESCOPE/RAW/TIME");

      // FRONT_RAW_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront+1, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/RAW/MULT");

      // BACK_RAW_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsBack+1, 0, fNumberOfStripsBack+1, "COMPTONTELESCOPE/RAW/MULT");

      // Time multiplicity
      // Front
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_FRONT_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront+1, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/RAW/MULT");
      // Back
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_BACK_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsBack+1, 0, fNumberOfStripsBack+1, "COMPTONTELESCOPE/RAW/MULT");

   } // end loop on number of detectors

    // CALORIMETER
    name = "CT"+NPL::itoa(i+1)+"_CALOR_RAW_TRIGGER";
    AddHisto1D(name, name, fCalorimeterNPixels, 1, fCalorimeterNPixels+1, "COMPTONTELESCOPE/RAW/CALORTRIGGER");
/*    for (int j = 0; j < 4; j++) {
      name = "CT"+NPL::itoa(i*4+j+1)+"_CALOR_RAW_TRIGGER";
      AddHisto1D(name, name, fCalorimeterNPixels, 1, fCalorimeterNPixels+1, "COMPTONTELESCOPE/RAW/CALORTRIGGER");
    }*/

    // BIDIM SUM
    name = "CT"+NPL::itoa(i+1)+"_RAW_SUM_BIDIM";
    AddHisto2D(name, name, 512, 0, 1024, 500, 500000, 1000000, "COMPTONTELESCOPE/RAW/SUM_BIDIM");
  } // end loop on number of telescopes
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::InitPreTreatedSpectra()
{
  string name;
  int ntot = (fNumberOfStripsFront + fNumberOfStripsBack);

  for (unsigned int i = 0; i < fNumberOfTelescope; i++) { // loop on number of telescope

    // DSSSD
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {

/*      // FRONT_E_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_E_CAL";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/ENERGY");

      // BACK_E_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_E_CAL";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/ENERGY");
*/
      // E Front and Back vs Strip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_BACK_E_CAL";
      AddHisto2D(name, name, ntot, 0, ntot, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/ENERGY");
 
      // Front-Back Energy Correlation
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FB_COR_CAL";
      AddHisto2D(name, name, 1400,0,1400, 1400,0,1400, "COMPTONTELESCOPE/CAL/ENERGYCOR");
   
      // Front E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONTECAL_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/ENERGYSPEC");

      // Back E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACKECAL_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/ENERGYSPEC");

      // FRONT_T_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_T_CAL";
      AddHisto1D(name, name, 5000, 0, 5e9, "COMPTONTELESCOPE/CAL/TIME");

/*      // BACK_T_CAL
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_T_CAL";
      AddHisto1D(name, name, 5000, 0, 5e9, "COMPTONTELESCOPE/CAL/TIME");
*/
      // FRONT_CAL_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront+1, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/CAL/MULT");

      // BACK_CAL_MULT
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsBack+1, 0, fNumberOfStripsBack+1, "COMPTONTELESCOPE/CAL/MULT");

/*      // Event type 1 multiplicity
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_EVENTTYPE1_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsBack+1, 0, fNumberOfStripsBack+1, "COMPTONTELESCOPE/CAL/EVENTTYPEMULT");

      // Event type 2 multiplicity
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_EVENTTYPE2_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsBack+1, 0, fNumberOfStripsBack+1, "COMPTONTELESCOPE/CAL/EVENTTYPEMULT");

      // Possible interstrip
      // Front
      // side by side strips
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_E_CLOSE_STRIP";
      AddHisto2D(name, name, 1400, 0, 1400, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/INTERSTRIP");
      // + energy sum = E Back
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_E_INTERSTRIP";
      AddHisto2D(name, name, 1400, 0, 1400, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/INTERSTRIP");
      // Back
      // side by side strips
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_E_CLOSE_STRIP";
      AddHisto2D(name, name, 1400, 0, 1400, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/INTERSTRIP");
      // + energy sum = E Front
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_E_INTERSTRIP";
      AddHisto2D(name, name, 1400, 0, 1400, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/INTERSTRIP");
      // FB correlation for interstrip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FB_COR_INTERSTRIP";
      AddHisto2D(name, name, 1400, 0, 1400, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/INTERSTRIP");      
      // Half E spectrum for interstrip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_HALFE_INTERSTRIP";
      AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/CAL/INTERSTRIP");
*/
      // Time multiplicity
      // Front
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_FRONT_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront+1, 0, fNumberOfStripsFront+1, "COMPTONTELESCOPE/CAL/MULT");
      // Back
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_BACK_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsBack+1, 0, fNumberOfStripsBack+1, "COMPTONTELESCOPE/CAL/MULT");

    }  // end loop on number of detectors
  }// end loop on number of telescopes
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::InitPhysicsSpectra()
{
  string name;
  int ntot = (fNumberOfStripsFront+fNumberOfStripsBack);

  // counters
/*  name = "CT_DSSSD_COUNTERS_EVTS";
  AddHisto1D(name, name, fNumberOfCounters, 0, fNumberOfCounters, "COMPTONTELESCOPE/PHY/COUNTER");

  name = "CT_DSSSD_COUNTERS_HITS";
  AddHisto1D(name, name, fNumberOfCounters, 0, fNumberOfCounters, "COMPTONTELESCOPE/PHY/COUNTER");
*/

  // loop on number of telescopes
  for (unsigned int i = 0 ; i < fNumberOfTelescope ; i++) {

    //// DSSSD
    // loop on number of detectors
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {

      // E Front + Back vs Strip
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_BACK_E_PHY";
      AddHisto2D(name, name, ntot, 0, ntot, 1400, 0, 1400, "COMPTONTELESCOPE/PHY/ENERGY");

      // front back Energy Correlation
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FB_COR_PHY";
      AddHisto2D(name, name, 1400,0,1400, 1400,0,1400, "COMPTONTELESCOPE/PHY/ENERGYCOR");
 
      // Front E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONTE_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/PHY/ENERGYSPEC");

      // Back E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACKE_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/PHY/ENERGYSPEC");

      // Half E spectrum
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_HALFE_SPECTRUM";
      AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/PHY/ENERGYSPEC");

      // Multiplicity
      int MultMax = fNumberOfStripsFront*fNumberOfStripsBack;
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_MULT_PHYS";
      AddHisto1D(name, name, MultMax+1, 0, MultMax+1, "COMPTONTELESCOPE/PHY/MULT");

      // Time
      // Front
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_T_PHY";
      AddHisto1D(name, name, 5000, -2000, 5e9, "COMPTONTELESCOPE/PHY/TIME");
/*      // Back
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_T_PHY";
      AddHisto1D(name, name, 5000, -2000, 5e9, "COMPTONTELESCOPE/PHY/TIME");
*/
      // Hit pattern
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_HIT_PHY";
      int ntotF = fNumberOfStripsFront * fNumberOfDetectors;
      int ntotB = fNumberOfStripsBack  * fNumberOfDetectors;
      AddHisto2D(name, name, ntotF, 0, ntotF, ntotB, 0, ntotB, "COMPTONTELESCOPE/PHY/HITPATTERN");
 
      // X-Y Impact Matrix
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_IMPACT_MATRIX";
      AddHisto2D(name, name, 40, -40, 40, 40, -40, 40, "COMPTONTELESCOPE/PHY/IMPACTMATRIX");

    }

    // X-Y Energy Correlation
//    name = "CT"+NPL::itoa(i+1)+"_XY_COR";
//    AddHisto2D(name, name,500,0,50,500,0,50, "COMPTONTELESCOPE/PHY");

    // Position on calorimeter
    name = "CT"+NPL::itoa(i+1)+"_CALOR_POS";
    AddHisto2D(name, name, 8, -24, 24, 8, -24, 24, "COMPTONTELESCOPE/PHY/CALOR_POS");

    // Calorimeter energy spectrum
    name = "CT"+NPL::itoa(i+1)+"_CALOR_SPECTRUM";
    AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/PHY/CALOR_SPECTRUM");

    // Sum spectrum
    name = "CT"+NPL::itoa(i+1)+"_SUM_SPECTRUM";
    AddHisto1D(name, name, 1400, 0, 1400, "COMPTONTELESCOPE/PHY/CALOR_SUM");

    // Bidim sum
    name = "CT"+NPL::itoa(i+1)+"_SUM_BIDIM";
    AddHisto2D(name, name, 1400, 0, 1400, 1000, 0, 500000, "COMPTONTELESCOPE/PHY/SUM_BIDIM");

    name = "CT"+NPL::itoa(i+1)+"_SUM_BIDIM_CAL";
    AddHisto2D(name, name, 1400, 0, 1400, 1400, 0, 1400, "COMPTONTELESCOPE/PHY/SUM_BIDIM");

    // DELTA T
    name = "CT"+NPL::itoa(i+1)+"_DELTA_T";
    AddHisto1D(name, name, 2000, -1000, 1000, "COMPTONTELESCOPE/PHY/DELTA_T");
  }
}



////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::FillRawSpectra(TComptonTelescopeData* RawData)
{
  string name;
  string family;

  // FRONT_E 
/*  for (unsigned int i = 0; i < RawData->GetCTTrackerFrontEMult(); i++) {    
    name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_E_RAW";
    family = "COMPTONTELESCOPE/RAW/ENERGY";
//    if (RawData->GetCTTrackerFrontTTime(i) > 3e9) {
//    if (RawData->GetCTTrackerBackEMult() == 0) {
    FillSpectra(family,name,
          RawData->GetCTTrackerFrontEStripNbr(i),
          RawData->GetCTTrackerFrontEEnergy(i));
//    }
  }

  // BACK_E
  for (unsigned int i = 0; i < RawData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerBackEDetectorNbr(i))+"_BACK_E_RAW";
    family = "COMPTONTELESCOPE/RAW/ENERGY";
//    if (RawData->GetCTTrackerBackTTime(i) > 3e9) {
//    if (RawData->GetCTTrackerFrontEMult() == 0) {
    FillSpectra(family,name,
          RawData->GetCTTrackerBackEStripNbr(i),
          RawData->GetCTTrackerBackEEnergy(i));
//    }
  }
*/
  // E Front Back vs Strip
  family = "COMPTONTELESCOPE/RAW/ENERGY";
  for (unsigned int i = 0; i < RawData->GetCTTrackerFrontEMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_BACK_E_RAW";
    FillSpectra(family, name,
          RawData->GetCTTrackerFrontEStripNbr(i),
          RawData->GetCTTrackerFrontEEnergy(i));
  }
  for (unsigned int i = 0; i < RawData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerBackEDetectorNbr(i))+"_FRONT_BACK_E_RAW";
    FillSpectra(family, name,
        fNumberOfStripsFront+RawData->GetCTTrackerBackEStripNbr(i),
        RawData->GetCTTrackerBackEEnergy(i));
  }

  // FRONT_T
  for (unsigned int i = 0; i < RawData->GetCTTrackerFrontTMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontTTowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerFrontTDetectorNbr(i))+"_FRONT_T_RAW";
    family = "COMPTONTELESCOPE/RAW/TIME";
    FillSpectra(family,name,RawData->GetCTTrackerFrontTTime(i));
  }
  // BACK_T
  for (unsigned int i = 0; i < RawData->GetCTTrackerBackTMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTTrackerBackTTowerNbr(i))+"_DSSSD"+NPL::itoa(RawData->GetCTTrackerBackTDetectorNbr(i))+"_BACK_T_RAW";
    family = "COMPTONTELESCOPE/RAW/TIME";
    FillSpectra(family,name,RawData->GetCTTrackerBackTTime(i));
  }

  // FRONT MULT
  int myMULT[fNumberOfTelescope][fNumberOfDetectors];
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }
  for (unsigned int i = 0 ; i < RawData->GetCTTrackerFrontEMult(); i++) { 
     myMULT[RawData->GetCTTrackerFrontETowerNbr(i)-1][RawData->GetCTTrackerFrontEDetectorNbr(i)-1] += 1;  
     //myMULT[RawData->GetCTTrackerFrontEDetectorNbr(i)-1] += 1;  
  }

  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_RAW_MULT";
      family= "COMPTONTELESCOPE/RAW/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // BACK MULT
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }
  for (unsigned int i = 0 ; i < RawData->GetCTTrackerBackEMult(); i++) {
     myMULT[RawData->GetCTTrackerBackETowerNbr(i)-1][RawData->GetCTTrackerBackEDetectorNbr(i)-1] += 1;  
  }
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_RAW_MULT";
      family= "COMPTONTELESCOPE/RAW/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // Time Multiplicity 
  // Front
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ;
    }
  }
  for (unsigned int i = 0 ; i < RawData->GetCTTrackerFrontTMult(); i++) {
    myMULT[RawData->GetCTTrackerFrontTTowerNbr(i)-1][RawData->GetCTTrackerFrontTDetectorNbr(i)-1] += 1;
  }
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_FRONT_RAW_MULT";
      family= "COMPTONTELESCOPE/RAW/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }
  // Back
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ;
    }
  }
  for (unsigned int i = 0 ; i < RawData->GetCTTrackerBackTMult(); i++) {
    myMULT[RawData->GetCTTrackerBackTTowerNbr(i)-1][RawData->GetCTTrackerBackTDetectorNbr(i)-1] += 1;
  }
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_BACK_RAW_MULT";
      family= "COMPTONTELESCOPE/RAW/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // CALORIMETER TRIGGERS
  for (unsigned int i = 0; i < RawData->GetCTCalorimeterTMult(); i++) {
    name = "CT"+NPL::itoa(RawData->GetCTCalorimeterTDetectorNbr(i))+"_CALOR_RAW_TRIGGER";
    family = "COMPTONTELESCOPE/RAW/CALORTRIGGER";
    FillSpectra(family, name, RawData->GetCTCalorimeterTChannelNbr(i)+1);
  }

  // SUM BIDIM
  if (RawData->GetCTTrackerFrontEMult() > 0 and RawData->GetCTCalorimeterTMult() > 0) {
    for (unsigned int i = 0 ; i < RawData->GetCTTrackerFrontEMult(); i++) {
      name = "CT"+NPL::itoa(RawData->GetCTTrackerFrontETowerNbr(i))+"_RAW_SUM_BIDIM";
      family = "COMPTONTELESCOPE/RAW/SUM_BIDIM";
      int sumE = 0;
      for (int j = 0; j<64; j++) {
        sumE += RawData->GetCTCalorimeterEEnergy(j);
      }
      FillSpectra(family, name,
          RawData->GetCTTrackerFrontEEnergy(i), sumE);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::FillPreTreatedSpectra(TComptonTelescopeData* PreTreatedData)
{
  string name;
  string family;
  
/*  // FRONT_E
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
//    if (PreTreatedData->GetCTTrackerFrontTTime(i) > 3e9) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_E_CAL";
    family = "COMPTONTELESCOPE/CAL/ENERGY";
    FillSpectra(family,name
      ,PreTreatedData->GetCTTrackerFrontEStripNbr(i), 
          PreTreatedData->GetCTTrackerFrontEEnergy(i));
//    }
  }

  // BACK_E
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerBackEMult(); i++) {
//    if (PreTreatedData->GetCTTrackerBackTTime(i) > 3e9) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackEDetectorNbr(i))+"_BACK_E_CAL";
    family = "COMPTONTELESCOPE/CAL/ENERGY";
    FillSpectra(family,name
      ,PreTreatedData->GetCTTrackerBackEStripNbr(i), 
          PreTreatedData->GetCTTrackerBackEEnergy(i));
//    }
  }
*/
  // E Front Back vs Strip
  family = "COMPTONTELESCOPE/CAL/ENERGY";
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONT_BACK_E_CAL";
    FillSpectra(family, name,
        PreTreatedData->GetCTTrackerFrontEStripNbr(i),
        PreTreatedData->GetCTTrackerFrontEEnergy(i));
  }
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackEDetectorNbr(i))+"_FRONT_BACK_E_CAL";
    FillSpectra(family, name,
        fNumberOfStripsFront+PreTreatedData->GetCTTrackerBackEStripNbr(i),
        PreTreatedData->GetCTTrackerBackEEnergy(i));
  }        

  // E Front Back correlation
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
    for (unsigned int j = 0; j < PreTreatedData->GetCTTrackerBackEMult(); j++) {
      if(PreTreatedData->GetCTTrackerFrontETowerNbr(i) == PreTreatedData->GetCTTrackerBackETowerNbr(j) &&
          PreTreatedData->GetCTTrackerFrontEDetectorNbr(i) == PreTreatedData->GetCTTrackerBackEDetectorNbr(j)) {
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(i))+"_FB_COR_CAL";
        family = "COMPTONTELESCOPE/CAL/ENERGYCOR";
        FillSpectra(family,name,PreTreatedData->GetCTTrackerFrontEEnergy(i), PreTreatedData->GetCTTrackerBackEEnergy(j));
      }
    }
  }

  // Energy spectrum
  // front
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(i))+"_FRONTECAL_SPECTRUM";
    family= "COMPTONTELESCOPE/CAL/ENERGYSPEC";
    FillSpectra(family, name, PreTreatedData->GetCTTrackerFrontEEnergy(i));
  }
  // back
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerBackEMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackEDetectorNbr(i))+"_BACKECAL_SPECTRUM";
    family= "COMPTONTELESCOPE/CAL/ENERGYSPEC";
    FillSpectra(family, name, PreTreatedData->GetCTTrackerBackEEnergy(i));
  }

  // FRONT_T
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerFrontTMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontTTowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontTDetectorNbr(i))+"_FRONT_T_CAL";
    family = "COMPTONTELESCOPE/CAL/TIME";
    FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontTTime(i));
  }
/*  // BACK_T
  for (unsigned int i = 0; i < PreTreatedData->GetCTTrackerBackTMult(); i++) {
    name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackTTowerNbr(i))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackTDetectorNbr(i))+"_BACK_T_CAL";
    family = "COMPTONTELESCOPE/CAL/TIME";
    FillSpectra(family,name, PreTreatedData->GetCTTrackerBackTTime(i));
  }
*/
  // FRONT MULT
  int myMULT[fNumberOfTelescope][fNumberOfDetectors];
  for( unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }
  for(unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerFrontEMult();i++){
    myMULT[PreTreatedData->GetCTTrackerFrontETowerNbr(i)-1][PreTreatedData->GetCTTrackerFrontEDetectorNbr(i)-1] += 1;  
  }
  for( unsigned int i = 0; i < fNumberOfTelescope; i++){
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_FRONT_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // BACK MULT
  for(unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }
  for(unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerBackEMult();i++){
    myMULT[PreTreatedData->GetCTTrackerBackETowerNbr(i)-1][PreTreatedData->GetCTTrackerBackEDetectorNbr(i)-1] += 1; 
  }
  for( unsigned int i = 0; i < fNumberOfTelescope; i++){
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_BACK_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

/*  // Event type
  // type 1 (mult front = mult back)
  for(unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }
  if (PreTreatedData->GetCTTrackerBackEMult() == PreTreatedData->GetCTTrackerFrontEMult()) {
    for(unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerBackEMult();i++){
      myMULT[PreTreatedData->GetCTTrackerBackETowerNbr(i)-1][PreTreatedData->GetCTTrackerBackEDetectorNbr(i)-1] += 1;
    }
  }
  for( unsigned int i = 0; i < fNumberOfTelescope; i++){
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_EVENTTYPE1_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }
  // type 2 (mult front = mult back +- 1)
  for(unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ; 
    }
  }
  if (PreTreatedData->GetCTTrackerBackEMult() == PreTreatedData->GetCTTrackerFrontEMult()+1 || PreTreatedData->GetCTTrackerBackEMult() == PreTreatedData->GetCTTrackerFrontEMult()-1) {
    if (PreTreatedData->GetCTTrackerBackEMult() > PreTreatedData->GetCTTrackerFrontEMult()) {
      for(unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerBackEMult();i++){
        myMULT[PreTreatedData->GetCTTrackerBackETowerNbr(i)-1][PreTreatedData->GetCTTrackerBackEDetectorNbr(i)-1] += 1;
      }
    }
    else if (PreTreatedData->GetCTTrackerFrontEMult() > PreTreatedData->GetCTTrackerBackEMult()) {
      for(unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerFrontEMult(); i++) {
        myMULT[PreTreatedData->GetCTTrackerFrontETowerNbr(i)-1][PreTreatedData->GetCTTrackerFrontEDetectorNbr(i)-1] += 1;
      }
    }
  }  
  for( unsigned int i = 0; i < fNumberOfTelescope; i++){
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_EVENTTYPE2_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/EVENTTYPEMULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

  // Possible interstrip
  // Front
  if (PreTreatedData->GetCTTrackerFrontEMult() == 2 && PreTreatedData->GetCTTrackerBackEMult() == 1) {
    if (PreTreatedData->GetCTTrackerFrontEStripNbr(0) == PreTreatedData->GetCTTrackerFrontEStripNbr(1)+1 || PreTreatedData->GetCTTrackerFrontEStripNbr(0) == PreTreatedData->GetCTTrackerFrontEStripNbr(1)-1) {
      name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(0))+"_FRONT_E_CLOSE_STRIP";
      family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
      FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontEEnergy(0),PreTreatedData->GetCTTrackerFrontEEnergy(1));

      // interstrip
      if (abs((PreTreatedData->GetCTTrackerBackEEnergy(0)-(PreTreatedData->GetCTTrackerFrontEEnergy(0)+PreTreatedData->GetCTTrackerFrontEEnergy(1)))/2.)< fStripEnergyMatchingNumberOfSigma*fStripEnergyMatchingSigma) {
        // EF(i+-1) vs EF(i)
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(0))+"_FRONT_E_INTERSTRIP";
        family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
        FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontEEnergy(0), PreTreatedData->GetCTTrackerFrontEEnergy(1));
      
        // EB vs EF(i) and EF(i+-1)
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(0))+"_FB_COR_INTERSTRIP";
        family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
        FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontEEnergy(0),PreTreatedData->GetCTTrackerBackEEnergy(0));
        FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontEEnergy(1),PreTreatedData->GetCTTrackerBackEEnergy(0));        

        // Spectrum
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(0))+"_HALFE_INTERSTRIP";
        family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
        FillSpectra(family,name, (PreTreatedData->GetCTTrackerBackEEnergy(0)+PreTreatedData->GetCTTrackerFrontEEnergy(0)+PreTreatedData->GetCTTrackerFrontEEnergy(1))/2.);
      }
    }
  }

  // Back
  if (PreTreatedData->GetCTTrackerBackEMult() == 2 && PreTreatedData->GetCTTrackerFrontEMult() == 1) {
    if (PreTreatedData->GetCTTrackerBackEStripNbr(0) == PreTreatedData->GetCTTrackerBackEStripNbr(1)+1 || PreTreatedData->GetCTTrackerBackEStripNbr(0) == PreTreatedData->GetCTTrackerBackEStripNbr(1)-1) {
      name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackEDetectorNbr(0))+"_BACK_E_CLOSE_STRIP";
      family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
      FillSpectra(family,name, PreTreatedData->GetCTTrackerBackEEnergy(0),PreTreatedData->GetCTTrackerBackEEnergy(1));
  
      // interstrip
      if (abs((PreTreatedData->GetCTTrackerFrontEEnergy(0)-(PreTreatedData->GetCTTrackerBackEEnergy(0)+PreTreatedData->GetCTTrackerBackEEnergy(1)))/2.)< fStripEnergyMatchingNumberOfSigma*fStripEnergyMatchingSigma) {
        // EB(i+-1) vs EB(i)
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerBackETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerBackEDetectorNbr(0))+"_BACK_E_INTERSTRIP";
        family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
        FillSpectra(family,name, PreTreatedData->GetCTTrackerBackEEnergy(0), PreTreatedData->GetCTTrackerBackEEnergy(1));
      
        // EF vs EB(i) and EB(i+-1)
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(0))+"_FB_COR_INTERSTRIP";
        family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
        FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontEEnergy(0),PreTreatedData->GetCTTrackerBackEEnergy(0));
        FillSpectra(family,name, PreTreatedData->GetCTTrackerFrontEEnergy(0),PreTreatedData->GetCTTrackerBackEEnergy(1));        
 
        // Spectrum
        name = "CT"+NPL::itoa(PreTreatedData->GetCTTrackerFrontETowerNbr(0))+"_DSSSD"+NPL::itoa(PreTreatedData->GetCTTrackerFrontEDetectorNbr(0))+"_HALFE_INTERSTRIP";
        family = "COMPTONTELESCOPE/CAL/INTERSTRIP";
        FillSpectra(family,name, (PreTreatedData->GetCTTrackerFrontEEnergy(0)+PreTreatedData->GetCTTrackerBackEEnergy(0)+PreTreatedData->GetCTTrackerBackEEnergy(1))/2.);
    }
    }
  }*/


  // Time Multiplicity 
  // Front
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ;
    }
  }
  for (unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerFrontTMult(); i++) {
    myMULT[PreTreatedData->GetCTTrackerFrontTTowerNbr(i)-1][PreTreatedData->GetCTTrackerFrontTDetectorNbr(i)-1] += 1;
  }
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_FRONT_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }
  // Back
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ;
    }
  }
  for (unsigned int i = 0 ; i < PreTreatedData->GetCTTrackerBackTMult(); i++) {
    myMULT[PreTreatedData->GetCTTrackerBackTTowerNbr(i)-1][PreTreatedData->GetCTTrackerBackTDetectorNbr(i)-1] += 1;
  }
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_T_BACK_CAL_MULT";
      family= "COMPTONTELESCOPE/CAL/MULT";
      FillSpectra(family,name,myMULT[i][j]);
    }
  }

 
}

////////////////////////////////////////////////////////////////////////////////
void TComptonTelescopeSpectra::FillPhysicsSpectra(TComptonTelescopePhysics* Physics){
  string name;
  string family;

  //// DSSSD

  // counters
/*  // for events
  name = "CT_DSSSD_COUNTERS_EVTS";
  family= "COMPTONTELESCOPE/PHY/COUNTER";
  for (unsigned int i = 0; i < fNumberOfCounters; i++) {
    FillSpectra(family, name, i, Physics->m_CounterEvt[i]);
  }
  // for hits 
  name = "CT_DSSSD_COUNTERS_HITS";
  family= "COMPTONTELESCOPE/PHY/COUNTER";
  for (unsigned int i = 0; i < fNumberOfCounters; i++) {
    FillSpectra(family, name, i, Physics->m_CounterHit[i]);
  }*/

  // Front + Back E per strip
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_FRONT_BACK_E_PHY";
    family= "COMPTONTELESCOPE/PHY/ENERGY";
    FillSpectra(family, name, Physics->GetFrontStrip(i), Physics->GetFrontEnergy(i));
    FillSpectra(family, name, fNumberOfStripsFront+Physics->GetBackStrip(i), Physics->GetBackEnergy(i));
  }

  // front back Energy correlation
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_FB_COR_PHY";
    family= "COMPTONTELESCOPE/PHY/ENERGYCOR";
    FillSpectra(family, name, Physics->GetFrontEnergy(i), Physics->GetBackEnergy(i)); 
  }

  // Energy spectrum
  // front
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_FRONTE_SPECTRUM";
    family= "COMPTONTELESCOPE/PHY/ENERGYSPEC";
    FillSpectra(family, name, Physics->GetFrontEnergy(i));
  }
  // back
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_BACKE_SPECTRUM";
    family= "COMPTONTELESCOPE/PHY/ENERGYSPEC";
    FillSpectra(family, name, Physics->GetBackEnergy(i));
  }
  // half
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_HALFE_SPECTRUM";
    family= "COMPTONTELESCOPE/PHY/ENERGYSPEC";
    FillSpectra(family, name, Physics->GetHalfEnergy(i));
  }

  // Multiplicity
  int myMULT[fNumberOfTelescope][fNumberOfDetectors];
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      myMULT[i][j] = 0 ;
    }
  }
  for (unsigned int i = 0 ; i < Physics->GetEventMultiplicity(); i++) {
    myMULT[Physics->GetTowerNumber(i)-1][Physics->GetDetectorNumber(i)-1] += 1;
  }
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    for (unsigned int j = 0; j < fNumberOfDetectors; j++) {
      name = "CT"+NPL::itoa(i+1)+"_DSSSD"+NPL::itoa(j+1)+"_MULT_PHYS";
      family= "COMPTONTELESCOPE/PHY/MULT";
      FillSpectra(family, name, myMULT[i][j]);
    }
  }

  // Time
  // Front
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_FRONT_T_PHY";
    family= "COMPTONTELESCOPE/PHY/TIME";
    FillSpectra(family,name,Physics->GetFrontTime(i));
  }
/*  // Back
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_BACK_T_PHY";
    family= "COMPTONTELESCOPE/PHY/TIME";
    FillSpectra(family,name,Physics->GetBackTime(i));
  }
*/
  // Hit pattern
  for (unsigned int i = 0; i < Physics->GetEventMultiplicity(); i++) {
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_HIT_PHY";
    family= "COMPTONTELESCOPE/PHY/HITPATTERN";
    FillSpectra(family,name,Physics->GetFrontStrip(i) + fNumberOfStripsFront*(Physics->GetDetectorNumber(i)-1), Physics->GetBackStrip(i) + fNumberOfStripsBack*(Physics->GetDetectorNumber(i)-1));    
  }

  // X-Y Impact Matrix
  for(unsigned int i = 0 ; i < Physics->GetEventMultiplicity(); i++){
    name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_DSSSD"+NPL::itoa(Physics->GetDetectorNumber(i))+"_IMPACT_MATRIX";
    family= "COMPTONTELESCOPE/PHY/IMPACTMATRIX";
    double x = Physics->GetPositionOfInteractionDSSSD(i).x();
    double y = Physics->GetPositionOfInteractionDSSSD(i).y();
    FillSpectra(family,name,x,y);
  }


  //// Calorimeter

  // Position on calorimeter
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    name = "CT"+NPL::itoa(i+1)+"_CALOR_POS";
    family = "COMPTONTELESCOPE/PHY/CALOR_POS";
    if (Physics->CalorPosX.size() == Physics->CalorPosY.size()) {
      for (int j = 0; j < Physics->CalorPosX.size(); j++) {
        FillSpectra(family, name, Physics->CalorPosX[j], Physics->CalorPosY[j]);
      }
    } else {
      cout << "Position not treated because size of x and y position vectors differs." << endl;
    }
  }

  // Calorimeters spectra
//  double energy = 0;
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    name = "CT"+NPL::itoa(i+1)+"_CALOR_SPECTRUM";
    family = "COMPTONTELESCOPE/PHY/CALOR_SPECTRUM";
//    energy = 0;
//    for (unsigned int j = 0; j < Physics->Calor_E.size(); j++) {
    for (unsigned int j = 0; j < Physics->Calor_E_calib.size(); j++) {
//      energy += Physics->Calor_E[j];
      FillSpectra(family, name, Physics->GetCalorEnergy(j));
    }
//    FillSpectra(family, name, energy);
  }

  // Sum spectrum
  double energy = 0;
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    name = "CT"+NPL::itoa(i+1)+"_SUM_SPECTRUM";
    family = "COMPTONTELESCOPE/PHY/CALOR_SUM";
    energy = 0;
//    for (unsigned int j = 0; j < Physics->Strip_E.size();j++) {
    for (unsigned int j = 0; j < Physics->GetEventMultiplicity();j++) {
//      energy += Physics->Strip_E[j];
      energy = Physics->GetHalfEnergy(j);
    }
//    for (unsigned int j = 0; j < Physics->Calor_E.size(); j++) {
    for (unsigned int j = 0; j < Physics->Calor_E_calib.size(); j++) {
      energy += Physics->GetCalorEnergy(j);
    }
    FillSpectra(family, name, energy);
  }

  // SUM BIDIM not calibrated
  for (unsigned int i = 0 ; i < Physics->GetEventMultiplicity(); i++) {
    for (unsigned int j = 0; j < Physics->Calor_E.size(); j++) {
      name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_SUM_BIDIM";
      family = "COMPTONTELESCOPE/PHY/SUM_BIDIM";
/*      int sumE = 0;
      for (int j = 0; j<64; j++) {
        sumE += RawData->GetCTCalorimeterEEnergy(j);
      }*/
      FillSpectra(family, name,
          Physics->GetFrontEnergy(i), Physics->Calor_E[j]);
    }
  }

  // SUM BIDIM calibrated
  for (unsigned int i = 0 ; i < Physics->GetEventMultiplicity(); i++) {
    for (unsigned int j = 0; j < Physics->Calor_E_calib.size(); j++) {
      name = "CT"+NPL::itoa(Physics->GetTowerNumber(i))+"_SUM_BIDIM_CAL";
      family = "COMPTONTELESCOPE/PHY/SUM_BIDIM";
//      int sumE = 0;
//      for (int j = 0; j<64; j++) {
//        sumE += RawData->GetCTCalorimeterEEnergy(j);
//      }
      FillSpectra(family, name,
          Physics->GetFrontEnergy(i), Physics->GetCalorEnergy(j));
    }
  }


  // DELTA T
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) {
    name = "CT"+NPL::itoa(i+1)+"_DELTA_T";
    family = "COMPTONTELESCOPE/PHY/DELTA_T";
    for (vector<int>::iterator it = Physics->deltaT.begin(); it != Physics->deltaT.end(); ++it) {
      FillSpectra(family, name, *it);
    }
  }

}

