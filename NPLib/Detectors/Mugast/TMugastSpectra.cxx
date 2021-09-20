/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta         contact address: matta@lpccaen.in2p3.fr *
 * Creation Date  : February 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for Mugast                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + first version                                                        *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// NPL
#include "TMugastSpectra.h"
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
TMugastSpectra::TMugastSpectra(){
  SetName("Mugast");
  fStripX=128;
  fStripY=128;
  fStripSecondLayer=16;
}

////////////////////////////////////////////////////////////////////////////////
TMugastSpectra::TMugastSpectra(map<int,int> TelescopeIndex){
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TMugastSpectra : Initalising control spectra for " 
      << TelescopeIndex.size() << " Telescopes" << endl
      << "************************************************" << endl ;

  SetName("Mugast");
  for(map<int,int>::iterator it=TelescopeIndex.begin() ; it!=TelescopeIndex.end() ; it++){
    fTelescopeToIndex[it->second-1]=it->first; 
    }

  fIndexToTelescope=TelescopeIndex;
  fNumberOfTelescope=fTelescopeToIndex.size();
  fStripX=128;
  fStripY=128;
  fStripSecondLayer=16;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}

////////////////////////////////////////////////////////////////////////////////
TMugastSpectra::~TMugastSpectra(){
}

////////////////////////////////////////////////////////////////////////////////
void TMugastSpectra::InitRawSpectra(){

  static string name;
  for (unsigned int i = 0; i < fTelescopeToIndex.size(); i++) { // loop on number of detectors
    TString nbr = NPL::itoa(fTelescopeToIndex[i]);
    
    // STRX_E_RAW
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRX_E_RAW";
    // AddHisto2D(name, name, fStripX, 1, fStripX+1, 512, 8192, 16384,  "Mugast/RAW/STRXE");
    AddHisto2D(name, name, fStripX, 1, fStripX+1, 10000, 7000, 17000,  "Mugast/RAW/STRXE");

    // STRY_E_RAW
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRY_E_RAW";
    // AddHisto2D(name, name, fStripY, 1, fStripY+1, 512, 0, 8192, "Mugast/RAW/STRYE");
    AddHisto2D(name, name, fStripY, 1, fStripY+1,10000, 0, 10000, "Mugast/RAW/STRYE");

    // STRX_T_RAW
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRX_T_RAW";
    AddHisto2D(name, name, fStripX, 1, fStripX+1, 512, 0, 8192, "Mugast/RAW/STRXT");

    // STRY_T_RAW
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRY_T_RAW";
    AddHisto2D(name, name, fStripY, 1, fStripY+1, 512, 0, 8192, "Mugast/RAW/STRYT");

    // SDLR_E_RAW
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_SDLR_E_RAW";
    AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 512, 0, 8192, "Mugast/RAW/SDLRE");

    // SDLR_T_RAW
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_SDLR_T_RAW";
    AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 512, 0, 8192, "Mugast/RAW/SDLRT");
  }
}

////////////////////////////////////////////////////////////////////////////////
void TMugastSpectra::InitPreTreatedSpectra(){

  string name;
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) { // loop on number of detectors
    // STRX_E_CAL
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRX_E_CAL";
    AddHisto2D(name, name, fStripX, 1, fStripX+1, 10000, 0, 50, "Mugast/CAL/STRXE");

    // STRY_E_CAL
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRY_E_CAL";
    AddHisto2D(name, name, fStripY, 1, fStripY+1, 10000, 0, 50, "Mugast/CAL/STRYE");

    // STR X-Y Correlation
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRXY_CORR_CAL";
    AddHisto2D(name, name, 500, 0, 50, 500, 0, 50, "Mugast/CAL/STRXY");

    // STRX_T_CAL
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRX_T_CAL";
    AddHisto2D(name, name, fStripX, 1, fStripX+1, 1000, 0, 1000, "Mugast/CAL/STRXT");

    // STRY_T_CAL
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_STRY_T_CAL";
    AddHisto2D(name, name, fStripY, 1, fStripY+1, 1000, 0, 1000, "Mugast/CAL/STRYT");

    // SDLR_E_CAL
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_SDLR_E_CAL";
    AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 500, 0, 50, "Mugast/CAL/SDLRE");

    // SDLR_T_CAL
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_SDLR_T_CAL";
    AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 500, 0, 50, "Mugast/CAL/SDLRT");

  }  // end loop on number of detectors

}

////////////////////////////////////////////////////////////////////////////////
void TMugastSpectra::InitPhysicsSpectra(){

  static string name;
  // X-Y Impact Matrix
  name = "MG_IMPACT_MATRIX";
  AddHisto2D(name, name,500,-150,150,500,-150,150, "Mugast/PHY");

  // X-Y Impact Matrix
  name = "MG_THETA_E";
  AddHisto2D(name, name,360,0,180,500,0,50,"Mugast/PHY");

  // ID Plot
  // E-TOF:
  name = "MG_E_TOF";
  AddHisto2D(name, name,500,0,50,500,200,1200,"Mugast/PHY");

  // SDLRE-DE:
  name = "MG_SDLRE_E";
  AddHisto2D(name, name,500,0,200,500,0,50,"Mugast/PHY");

  // Etot-DE:
  name = "MG_Etot_E";
  AddHisto2D(name, name,500,0,500,500,0,50,"Mugast/PHY");


  // ID plot detector by detector
  for (unsigned int i = 0; i < fNumberOfTelescope; i++) { // loop on number of detectors
    // E-TOF:
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_E_TOF";
    AddHisto2D(name, name,500,0,50,500,200,1200,"Mugast/PHY");

    // SDLRE-DE:
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_SDLRE_E";
    AddHisto2D(name, name,500,0,200,500,0,50,"Mugast/PHY");

    // Etot-DE:
    name = "MG"+NPL::itoa(fTelescopeToIndex[i])+"_Etot_E";
    AddHisto2D(name, name,500,0,500,500,0,50,"Mugast/PHY");
  }

}



////////////////////////////////////////////////////////////////////////////////
void TMugastSpectra::FillRawSpectra(TMugastData* RawData){

  static string name;
  static string family;

  // STRX_E 
  unsigned int size = RawData->GetDSSDXEMult();
  for (unsigned int i = 0; i < size ; i++) {
    name = "MG"+NPL::itoa( RawData->GetDSSDXEDetectorNbr(i))+"_STRX_E_RAW";
    family = "Mugast/RAW/STRXE";

    FillSpectra(family,name
        ,RawData->GetDSSDXEStripNbr(i), 
        RawData->GetDSSDXEEnergy(i));
  }

  // STRY_E
  size = RawData->GetDSSDYEMult();
  for (unsigned int i = 0; i < size ; i++) {
    name = "MG"+NPL::itoa( RawData->GetDSSDYEDetectorNbr(i) )+"_STRY_E_RAW";
    family = "Mugast/RAW/STRYE";

    FillSpectra(family,name
        ,RawData->GetDSSDYEStripNbr(i),
        RawData->GetDSSDYEEnergy(i));
  }

  // STRX_T
  size = RawData->GetDSSDXTMult();
  for (unsigned int i = 0; i < size; i++) {
    name = "MG"+NPL::itoa(RawData->GetDSSDXTDetectorNbr(i))+"_STRX_T_RAW";
    family = "Mugast/RAW/STRXT";

    FillSpectra(family,name
        ,RawData->GetDSSDXTStripNbr(i),
        RawData->GetDSSDXTTime(i));
  }

  // STRY_T
  size = RawData->GetDSSDYTMult();
  for (unsigned int i = 0; i < size ; i++) {
    name = "MG"+NPL::itoa( RawData->GetDSSDYTDetectorNbr(i))+"_STRY_T_RAW";
    family = "Mugast/RAW/STRYT";

    FillSpectra(family,name
        ,RawData->GetDSSDYTStripNbr(i),
        RawData->GetDSSDYTTime(i));
  }

  // SDLR_E
  size = RawData->GetSecondLayerEMult();
  for (unsigned int i = 0; i < size ; i++) {
    name = "MG"+NPL::itoa( RawData->GetSecondLayerEDetectorNbr(i))+"_SDLR_E_RAW";
    family = "Mugast/RAW/SDLRE";

    FillSpectra(family,name
        ,RawData->GetSecondLayerEStripNbr(i),
        RawData->GetSecondLayerEEnergy(i));
  }

  // SDLR_T
  size = RawData->GetSecondLayerTMult();
  for (unsigned int i = 0; i < size ; i++) {
    name = "MG"+NPL::itoa(RawData->GetSecondLayerTDetectorNbr(i))+"_SDLR_T_RAW";
    family = "Mugast/RAW/SDLRT";

    FillSpectra(family,name
        ,RawData->GetSecondLayerTStripNbr(i),
        RawData->GetSecondLayerTTime(i));
  }

}

////////////////////////////////////////////////////////////////////////////////
void TMugastSpectra::FillPreTreatedSpectra(TMugastData* PreTreatedData){

  static string name ;
  static string family;
  // STRX_E
  unsigned int sizeX = PreTreatedData->GetDSSDXEMult();
  for (unsigned int i = 0; i < sizeX ; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDXEDetectorNbr(i))+"_STRX_E_CAL";
    family = "Mugast/CAL/STRXE";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDXEStripNbr(i), 
        PreTreatedData->GetDSSDXEEnergy(i));
  }
  // STRY_E
  unsigned int sizeY = PreTreatedData->GetDSSDYEMult();
  for (unsigned int i = 0; i < sizeY; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDYEDetectorNbr(i))+"_STRY_E_CAL";
    family = "Mugast/CAL/STRYE";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDYEStripNbr(i), 
        PreTreatedData->GetDSSDYEEnergy(i));
  }

  // STR XY Correlation
  for (unsigned int i = 0; i < sizeX; i++) {
    for (unsigned int j = 0; j < sizeY; j++) {
      if(PreTreatedData->GetDSSDXEDetectorNbr(i)==PreTreatedData->GetDSSDYEDetectorNbr(j))
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDXEDetectorNbr(i) )+"_STRXY_CORR_CAL";
    family = "Mugast/CAL/STRXY";
    FillSpectra(family,name
        ,PreTreatedData->GetDSSDXEEnergy(i),
        PreTreatedData->GetDSSDYEEnergy(j));
    }
  }


  // STRX_T
  unsigned int size = PreTreatedData->GetDSSDXTMult();
  for (unsigned int i = 0; i < size; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDXTDetectorNbr(i))+"_STRX_T_CAL";
    family = "Mugast/CAL/STRXT";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDXTStripNbr(i), 
        PreTreatedData->GetDSSDXTTime(i));
  }
  // STRY_T
  size = PreTreatedData->GetDSSDYTMult();
  for (unsigned int i = 0; i < size; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDYTDetectorNbr(i))+"_STRY_T_CAL";
    family = "Mugast/CAL/STRYT";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDYTStripNbr(i), 
        PreTreatedData->GetDSSDYTTime(i));
  }
  // SDLR_E
  size = PreTreatedData->GetSecondLayerEMult();
  for (unsigned int i = 0; i < size; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetSecondLayerEDetectorNbr(i))+"_SDLR_E_CAL";
    family = "Mugast/CAL/SDLRE";

    FillSpectra(family,name
        ,PreTreatedData->GetSecondLayerEStripNbr(i), 
        PreTreatedData->GetSecondLayerEEnergy(i));
  }
  // SDLR_T
  size = PreTreatedData->GetSecondLayerTMult();
  for (unsigned int i = 0; i < size; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetSecondLayerTDetectorNbr(i))+"_SDLR_T_CAL";
    family = "Mugast/CAL/SDLRT";

    FillSpectra(family,name
        ,PreTreatedData->GetSecondLayerTStripNbr(i), 
        PreTreatedData->GetSecondLayerTTime(i));
  }
 
  //E-SDLR ID
  family = "Mugast/CAL/ID";
  size = PreTreatedData->GetSecondLayerEMult();
  for (unsigned int i = 0; i < sizeX; i++) {
    for (unsigned int j = 0; j < size; j++) {

      if(PreTreatedData->GetDSSDXEDetectorNbr(i) == PreTreatedData->GetSecondLayerEDetectorNbr(j)){ 
        name = "MG"+NPL::itoa( PreTreatedData->GetDSSDXEDetectorNbr(i))+"_SDLR"+NPL::itoa(PreTreatedData->GetSecondLayerEStripNbr(j))+"_CAL_ID";

        FillSpectra(family,name
            ,PreTreatedData->GetSecondLayerEEnergy(j), 
            PreTreatedData->GetDSSDXEEnergy(i));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void TMugastSpectra::FillPhysicsSpectra(TMugastPhysics* Physics){

  static string name;
  static string family= "Mugast/PHY";
  // X-Y Impact Matrix

  for(unsigned int i = 0 ; i < Physics->DSSD_E.size(); i++){
    name = "MG_IMPACT_MATRIX";
    double x = Physics->GetPositionOfInteraction(i).x();
    double y = Physics->GetPositionOfInteraction(i).y();
    FillSpectra(family,name,x,y);

    name = "MG_THETA_E";
    double Theta = Physics->GetPositionOfInteraction(i).Angle(TVector3(0,0,1));
    Theta = Theta/deg;
    FillSpectra(family,name,Theta,Physics->DSSD_E[i]);

    // Fill only for particle stopped in the first stage
    if(Physics->SecondLayer_E[i]<0 ){
      // E-TOF:
      name = "MG_E_TOF";
      FillSpectra(family,name,Physics->DSSD_E[i],Physics->DSSD_T[i]);

      name = "MG"+NPL::itoa(Physics->TelescopeNumber[i])+"_E_TOF";
      FillSpectra(family,name,Physics->DSSD_E[i],Physics->DSSD_T[i]);
    }

    double Etot=0;
    if(Physics->SecondLayer_E[i]>0){
      name = "MG_SDLRE_E";
      Etot = Physics->SecondLayer_E[i];
      FillSpectra(family,name,Physics->SecondLayer_E[i],Physics->DSSD_E[i]);

      name = "MG"+NPL::itoa(Physics->TelescopeNumber[i])+"_SDLRE_E";
      FillSpectra(family,name,Physics->SecondLayer_E[i],Physics->DSSD_E[i]);
    }

    if(Etot>0){
      name = "MG_Etot_E";
      FillSpectra(family,name,Etot,Physics->DSSD_E[i]);
      name = "MG"+NPL::itoa(Physics->TelescopeNumber[i])+"_Etot_E";
      FillSpectra(family,name,Etot,Physics->DSSD_E[i]);
    }
  }

}

