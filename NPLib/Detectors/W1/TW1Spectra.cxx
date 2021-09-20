/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville         address: deserevi@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : November 2015                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for W1                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// class header
#include "TW1Spectra.h"

// C++ headers
#include <iostream>  
#include <string>
using namespace std;

// NPTool headers
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TW1Spectra::TW1Spectra()
   : fNumberOfDetectors(0),
     fNumberOfStripsFront(16),
     fNumberOfStripsBack(16),
     fNumberOfCounters(10)
{
   SetName("W1");
}



////////////////////////////////////////////////////////////////////////////////
TW1Spectra::TW1Spectra(unsigned int NumberOfDetectors)
{
   if (NPOptionManager::getInstance()->GetVerboseLevel() > 0)
      cout << "*********************************************************" << endl
           << "TW1Spectra : Initializing control spectra for " 
           << NumberOfDetectors << " detectors" << endl
           << "*********************************************************" << endl ;
   SetName("W1");
   fNumberOfDetectors = NumberOfDetectors;
   fNumberOfStripsFront = 16;
   fNumberOfStripsBack  = 16;
   fNumberOfCounters    = 10;

   InitRawSpectra();
   InitPreTreatedSpectra();
   InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TW1Spectra::~TW1Spectra()
{
}



////////////////////////////////////////////////////////////////////////////////
void TW1Spectra::InitRawSpectra()
{
   static string name;
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
      // STR_FRONT_E_RAW
      name = "W1_D"+NPL::itoa(i+1)+"_STR_FRONT_E_RAW";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 4096, 0, 4096, "W1/RAW/STR_FRONT_E");

      // STR_BACK_E_RAW
      name = "W1_D"+NPL::itoa(i+1)+"_STR_BACK_E_RAW";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 4096, 0, 4096, "W1/RAW/STR_BACK_E");

      // STR_FRONT_T_RAW
      name = "W1_D"+NPL::itoa(i+1)+"_STR_FRONT_T_RAW";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 4096, 0, 4096, "W1/RAW/STR_FRONT_T");

      // STR_BACK_T_RAW
      name = "W1_D"+NPL::itoa(i+1)+"_STR_BACK_T_RAW";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 4096, 0, 4096, "W1/RAW/STR_BACK_T");

      // STR_FRONT_EMAX_RAW
      name = "W1_D"+NPL::itoa(i+1)+"_STR_FRONT_EMAX_RAW";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 4096, 0, 4096, "W1/RAW/STR_FRONT_EMAX");

      // STR_BACK_EMAX_Raw
      name = "W1_D"+NPL::itoa(i+1)+"_STR_BACK_EMAX_RAW";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 4096, 0, 4096, "W1/RAW/STR_BACK_EMAX");

      // STR_FRONT_RAW_MULT
      name = "W1_D"+NPL::itoa(i+1)+"_STR_FRONT_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 1, fNumberOfStripsFront+1, "W1/RAW/MULT");

      // STR_BACK_RAW_MULT
      name = "W1_D"+NPL::itoa(i+1)+"_STR_BACK_RAW_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 1, fNumberOfStripsFront+1, "W1/RAW/MULT");
   } // end loop on number of detectors

   // One single STR_E_RAW spectrum (E vs absolute strip number)
   name = "W1_STR_E_RAW";
   Int_t ntot = (fNumberOfStripsFront+fNumberOfStripsBack) * fNumberOfDetectors;
   AddHisto2D(name, name, ntot, 0, ntot, 4096, 0, 4096, "W1/RAW");

   // One single STR_T_RAW spectrum (T vs absolute strip number)
   name = "W1_STR_T_RAW";
   ntot = (fNumberOfStripsFront+fNumberOfStripsBack) * fNumberOfDetectors;
   AddHisto2D(name, name, ntot, 0, ntot, 4096, 0, 4096, "W1/RAW");

   // Energy vs time correlation plot
   name = "W1_ET_COR_RAW";
   ntot = (fNumberOfStripsFront+fNumberOfStripsBack) * fNumberOfDetectors;
   AddHisto2D(name, name, ntot, 0, ntot, ntot, 0, ntot, "W1/RAW");

   // One single Energy hit pattern
   name = "W1_HIT_E_RAW";
   Int_t ntotF = fNumberOfStripsFront * fNumberOfDetectors;
   Int_t ntotB = fNumberOfStripsBack  * fNumberOfDetectors;
   AddHisto2D(name, name, ntotF, 0, ntotF, ntotB, 0, ntotB, "W1/RAW");

   // One single Time hit pattern
   name = "W1_HIT_T_RAW";
   ntotF = fNumberOfStripsFront * fNumberOfDetectors;
   ntotB = fNumberOfStripsBack  * fNumberOfDetectors;
   AddHisto2D(name, name, ntotF, 0, ntotF, ntotB, 0, ntotB, "W1/RAW");
}



////////////////////////////////////////////////////////////////////////////////
void TW1Spectra::InitPreTreatedSpectra()
{
   static string name;
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
      // STR_FRONT_E_CAL
      name = "W1_D"+NPL::itoa(i+1)+"_STR_FRONT_E_CAL";
      AddHisto2D(name, name, fNumberOfStripsFront, 0, fNumberOfStripsFront, 4000, 0, 8, "W1/CAL/STR_FRONT_E");

      // STR_BACK_E_CAL
      name = "W1_D"+NPL::itoa(i+1)+"_STR_BACK_E_CAL";
      AddHisto2D(name, name, fNumberOfStripsBack, 0, fNumberOfStripsBack, 4000, 0, 8, "W1/CAL/STR_BACK_E");

      // STR_FRONT_CAL_MULT
      name = "W1_D"+NPL::itoa(i+1)+"_STR_FRONT_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 1, fNumberOfStripsFront+1, "W1/CAL/MULT");

      // STR_BACK_CAL_MULT
      name = "W1_D"+NPL::itoa(i+1)+"_STR_BACK_CAL_MULT";
      AddHisto1D(name, name, fNumberOfStripsFront, 1, fNumberOfStripsFront+1, "W1/CAL/MULT");

      // Front-Back Energy Correlation
      name = "W1_D"+NPL::itoa(i+1)+"_FB_COR";
      AddHisto2D(name, name, 2000, 0, 8, 2000, 0, 8, "W1/CAL/FB"); 
   }  // end loop on number of detectors

   // One single STR_E_CAL spectrum
   name = "W1_STR_E_CAL";
   Int_t ntot = (fNumberOfStripsFront+fNumberOfStripsBack) * fNumberOfDetectors;
   AddHisto2D(name, name, ntot, 0, ntot, 2000, 0, 8, "W1/CAL");
 
   // One single STR_T_CAL spectrum
   name = "W1_STR_T_CAL";
   ntot = (fNumberOfStripsFront+fNumberOfStripsBack) * fNumberOfDetectors;
   AddHisto2D(name, name, ntot, 0, ntot, 4096, 0, 4096, "W1/CAL");
 
   // Energy vs time correlation plot
   name = "W1_ET_COR_CAL";
   ntot = (fNumberOfStripsFront+fNumberOfStripsBack) * fNumberOfDetectors;
   AddHisto2D(name, name, ntot, 0, ntot, ntot, 0, ntot, "W1/CAL");

   // One single Energy hit pattern
   name = "W1_HIT_E_CAL";
   Int_t ntotF = fNumberOfStripsFront * fNumberOfDetectors;
   Int_t ntotB = fNumberOfStripsBack  * fNumberOfDetectors;
   AddHisto2D(name, name, ntotF, 0, ntotF, ntotB, 0, ntotB, "W1/CAL");

   // One single Time hit pattern
   name = "W1_HIT_T_CAL";
   ntotF = fNumberOfStripsFront * fNumberOfDetectors;
   ntotB = fNumberOfStripsBack  * fNumberOfDetectors;
   AddHisto2D(name, name, ntotF, 0, ntotF, ntotB, 0, ntotB, "W1/CAL");
}



////////////////////////////////////////////////////////////////////////////////
void TW1Spectra::InitPhysicsSpectra()
{
   static string name;
   name = "W1_CUTS";
   AddHisto1D(name, name, fNumberOfCounters, 1, fNumberOfCounters+1, "W1/PHY");

   // Kinematic Plot 
   name = "W1_THETA_E";
   AddHisto2D(name, name,360,0,180,500,0,50,"W1/PHY");
 
   // Energy vs time correlation plot
   name = "W1_ET_COR_PHY";
   Int_t ntot = (fNumberOfStripsFront+fNumberOfStripsBack) * fNumberOfDetectors;
   AddHisto2D(name, name, ntot, 0, ntot, ntot, 0, ntot, "W1/PHY");

   // One single Energy hit pattern
   name = "W1_HIT_PHY";
   Int_t ntotF = fNumberOfStripsFront * fNumberOfDetectors;
   Int_t ntotB = fNumberOfStripsBack  * fNumberOfDetectors;
   AddHisto2D(name, name, ntotF, 0, ntotF, ntotB, 0, ntotB, "W1/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TW1Spectra::FillRawSpectra(TW1Data* RawData)
{
   static string index;

   // STR_FRONT_E 
   unsigned int mysize = RawData->GetFrontEMult();
   double EFMAX = 0 ;
   int SFMAX = 0;
   int DFMAX = 0 ;
   for (unsigned int i = 0; i < mysize; i++) {
      index = "W1/RAW/STR_FRONT_E/W1_D"+NPL::itoa(RawData->GetFrontEDetectorNbr(i))+"_STR_FRONT_E_RAW";
      if(RawData->GetFrontEEnergy(i) > EFMAX){
         EFMAX = RawData->GetFrontEEnergy(i);
         SFMAX = RawData->GetFrontEStripNbr(i);
         DFMAX = RawData->GetFrontEDetectorNbr(i);
      }
      FillSpectra(index,RawData->GetFrontEStripNbr(i), RawData->GetFrontEEnergy(i));
   }

   if (DFMAX != 0) {
      index = "W1/RAW/STR_FRONT_EMAX/W1_D"+NPL::itoa(DFMAX)+"_STR_FRONT_EMAX_RAW";
      FillSpectra(index,SFMAX, EFMAX);
   }

   // STR_BACK_E
   mysize = RawData->GetBackEMult();
   double EBMAX = 0 ;
   int SBMAX = 0;
   int DBMAX = 0 ;

   for (unsigned int i = 0; i < mysize; i++) {
      index = "W1/RAW/STR_BACK_E/W1_D"+NPL::itoa( RawData->GetBackEDetectorNbr(i) )+"_STR_BACK_E_RAW";
      if(RawData->GetBackEEnergy(i) > EBMAX){
         EBMAX = RawData->GetBackEEnergy(i);
         SBMAX = RawData->GetBackEStripNbr(i);
         DBMAX = RawData->GetBackEDetectorNbr(i);
      }
      FillSpectra(index,RawData->GetBackEStripNbr(i), RawData->GetBackEEnergy(i));
   }

   if (DBMAX != 0) {
      index = "W1/RAW/STR_BACK_EMAX/W1_D"+NPL::itoa(DBMAX)+"_STR_BACK_EMAX_RAW";
      FillSpectra(index,SBMAX, EBMAX);
   }

   // STR_FRONT_T 
   mysize = RawData->GetFrontTMult();
   for (unsigned int i = 0; i < mysize; i++) {
      index = "W1/RAW/STR_FRONT_T/W1_D"+NPL::itoa(RawData->GetFrontTDetectorNbr(i))+"_STR_FRONT_T_RAW";
      FillSpectra(index,RawData->GetFrontTStripNbr(i), RawData->GetFrontTTime(i));
   }

   // STR_BACK_T 
   mysize = RawData->GetBackTMult();
   for (unsigned int i = 0; i < mysize; i++) {
      index = "W1/RAW/STR_BACK_T/W1_D"+NPL::itoa(RawData->GetBackTDetectorNbr(i))+"_STR_BACK_T_RAW";
      FillSpectra(index,RawData->GetBackTStripNbr(i), RawData->GetBackTTime(i));
   }

   // STR_E_RAW
   index = "W1/RAW/W1_STR_E_RAW";
   UShort_t multFront = RawData->GetFrontEMult();
   for (UShort_t i = 0; i < multFront; ++i) {   // loop on front strips
      FillSpectra(index,RawData->GetFrontEStripNbr(i) + 2*fNumberOfStripsFront*(RawData->GetFrontEDetectorNbr(i)-1), RawData->GetFrontEEnergy(i));
   } // end loop on front strips
   UShort_t multBack = RawData->GetBackEMult();
   for (UShort_t i = 0; i < multBack; ++i) {   // loop on front strips
      FillSpectra(index,RawData->GetBackEStripNbr(i) + (2*(RawData->GetBackEDetectorNbr(i)-1)+1)*fNumberOfStripsBack, RawData->GetBackEEnergy(i));
   } // end loop on front strips

   // STR_T_RAW
   index = "W1/RAW/W1_STR_T_RAW";
   multFront = RawData->GetFrontTMult();
   for (UShort_t i = 0; i < multFront; ++i) {   // loop on front strips
      FillSpectra(index,RawData->GetFrontTStripNbr(i) + 2*fNumberOfStripsFront*(RawData->GetFrontTDetectorNbr(i)-1), RawData->GetFrontTTime(i));
   } // end loop on front strips
   multBack = RawData->GetBackTMult();
   for (UShort_t i = 0; i < multBack; ++i) {   // loop on front strips
      FillSpectra(index,RawData->GetBackTStripNbr(i) + (2*(RawData->GetBackTDetectorNbr(i)-1)+1)*fNumberOfStripsBack, RawData->GetBackTTime(i));
   } // end loop on front strips

   // W1_ET_COR_RAW
   index = "W1/RAW/W1_ET_COR_RAW";
   UShort_t multFrontE = RawData->GetFrontEMult();
   for (UShort_t i = 0; i < multFrontE; ++i) {   // loop on front strips (E)
      UShort_t multFrontT = RawData->GetFrontTMult();
      for (UShort_t j = 0; j < multFrontT; ++j) {   // loop on front strips (T)
         FillSpectra(index,RawData->GetFrontEStripNbr(i) + 2*fNumberOfStripsFront*(RawData->GetFrontEDetectorNbr(i)-1), RawData->GetFrontTStripNbr(i) + 2*fNumberOfStripsFront*(RawData->GetFrontTDetectorNbr(i)-1));
      } // end loop on front strips (T)
   } // end loop on front strips (E)
   UShort_t multBackE = RawData->GetBackEMult();
   for (UShort_t i = 0; i < multBackE; ++i) {   // loop on front strips (E)
      UShort_t multBackT = RawData->GetBackTMult();
      for (UShort_t j = 0; j < multBackT; ++j) {   // loop on front strips (T)
         FillSpectra(index,RawData->GetBackEStripNbr(i) + (2*(RawData->GetBackEDetectorNbr(i)-1)+1)*fNumberOfStripsBack, RawData->GetBackTStripNbr(i) + (2*(RawData->GetBackTDetectorNbr(i)-1)+1)*fNumberOfStripsBack);
      } // end loop on front strips (T)
   } // end loop on front strips (E)
   
   // W1_HIT_E_RAW
   index = "W1/RAW/W1_HIT_E_RAW";
   for (UShort_t i = 0; i < RawData->GetFrontEMult(); ++i) {   // loop on front strips
      for (UShort_t j = 0; j < RawData->GetBackEMult(); ++j) { // loop on back strips
         FillSpectra(index,RawData->GetFrontEStripNbr(i) + fNumberOfStripsFront*(RawData->GetFrontEDetectorNbr(i)-1), RawData->GetBackEStripNbr(i) + fNumberOfStripsBack*(RawData->GetBackEDetectorNbr(i)-1));
      }
   } // end loop on front strips

   // W1_HIT_T_RAW
   index = "W1/RAW/W1_HIT_T_RAW";
   for (UShort_t i = 0; i < RawData->GetFrontTMult(); ++i) {   // loop on front strips
      for (UShort_t j = 0; j < RawData->GetBackTMult(); ++j) { // loop on back strips
         FillSpectra(index,RawData->GetFrontTStripNbr(i) + fNumberOfStripsFront*(RawData->GetFrontTDetectorNbr(i)-1), RawData->GetBackTStripNbr(i) + fNumberOfStripsBack*(RawData->GetBackTDetectorNbr(i)-1));
      }
   } // end loop on front strips

   // STR_FRONT MULT
   int myMULT[fNumberOfDetectors];
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) myMULT[i] = 0;

   for (unsigned int i = 0; i < RawData->GetFrontEMult();i++) {
      myMULT[RawData->GetFrontEDetectorNbr(i)-1] += 1;  
   }

   for (unsigned int i = 0; i < fNumberOfDetectors; i++) {
      index = "W1/RAW/MULT/W1_D"+NPL::itoa(i+1)+"_STR_FRONT_RAW_MULT";
      FillSpectra(index,myMULT[i]);
   }

   // STR_BACK MULT
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) myMULT[i] = 0;

   mysize = RawData->GetBackEMult();
   for (unsigned int i = 0; i < mysize;i++) {
      myMULT[RawData->GetBackEDetectorNbr(i)-1] += 1;  
   }

   for (unsigned int i = 0; i < fNumberOfDetectors; i++) {
      index= "W1/RAW/MULT/W1_D"+NPL::itoa(i+1)+"_STR_BACK_RAW_MULT";
      FillSpectra(index,myMULT[i]);
   }
}



////////////////////////////////////////////////////////////////////////////////
void TW1Spectra::FillPreTreatedSpectra(TW1Data* PreTreatedData)
{
   static string index;

   // Front-Back
   unsigned int mysizeF = PreTreatedData->GetFrontEMult();
   unsigned int mysizeB = PreTreatedData->GetBackEMult();

   for (unsigned int i = 0; i < mysizeF; i++) {
      for (unsigned int j = 0; j < mysizeB; j++) {
         if(PreTreatedData->GetFrontEDetectorNbr(i)==PreTreatedData->GetBackEDetectorNbr(j)){
            index="W1/CAL/FB/W1_D"+NPL::itoa(PreTreatedData->GetFrontEDetectorNbr(i))+"_FB_COR";
            FillSpectra(index,PreTreatedData->GetFrontEEnergy(i), PreTreatedData->GetBackEEnergy(j));
         }
      }
   } 

   // STR_FRONT_E
   unsigned int mysize = PreTreatedData->GetFrontEMult();
   for (unsigned int i = 0; i < mysize; i++) {
      index = "W1/CAL/STR_FRONT_E/W1_D"+NPL::itoa(PreTreatedData->GetFrontEDetectorNbr(i))+"_STR_FRONT_E_CAL";
      FillSpectra(index,PreTreatedData->GetFrontEStripNbr(i), PreTreatedData->GetFrontEEnergy(i));
   }
   // STR_BACK_E
   mysize = PreTreatedData->GetBackEMult();
   for (unsigned int i = 0; i < mysize; i++) {
      index = "W1/CAL/STR_BACK_E/W1_D"+NPL::itoa( PreTreatedData->GetBackEDetectorNbr(i))+"_STR_BACK_E_CAL";
      FillSpectra(index,PreTreatedData->GetBackEStripNbr(i), PreTreatedData->GetBackEEnergy(i));
   }


   // STR_E_CAL
   index = "W1/CAL/W1_STR_E_CAL";
   UShort_t multFront = PreTreatedData->GetFrontEMult();
   for (UShort_t i = 0; i < multFront; ++i) {   // loop on front strips
      FillSpectra(index,PreTreatedData->GetFrontEStripNbr(i) + 2*fNumberOfStripsFront*(PreTreatedData->GetFrontEDetectorNbr(i)-1), PreTreatedData->GetFrontEEnergy(i));
   } // end loop on front strips
   UShort_t multBack = PreTreatedData->GetBackEMult();
   for (UShort_t i = 0; i < multBack; ++i) {   // loop on front strips
      FillSpectra(index,PreTreatedData->GetBackEStripNbr(i) + (2*(PreTreatedData->GetBackEDetectorNbr(i)-1)+1)*fNumberOfStripsBack, PreTreatedData->GetBackEEnergy(i));
   } // end loop on front strips
 

   // STR_T_CAL
   index = "W1/CAL/W1_STR_T_CAL";
   multFront = PreTreatedData->GetFrontTMult();
   for (UShort_t i = 0; i < multFront; ++i) {   // loop on front strips
      FillSpectra(index,PreTreatedData->GetFrontTStripNbr(i) + 2*fNumberOfStripsFront*(PreTreatedData->GetFrontTDetectorNbr(i)-1), PreTreatedData->GetFrontTTime(i));
   } // end loop on front strips
   multBack = PreTreatedData->GetBackTMult();
   for (UShort_t i = 0; i < multBack; ++i) {   // loop on front strips
      FillSpectra(index,PreTreatedData->GetBackTStripNbr(i) + (2*(PreTreatedData->GetBackTDetectorNbr(i)-1)+1)*fNumberOfStripsBack, PreTreatedData->GetBackTTime(i));
   } // end loop on front strips
 
   // W1_ET_COR_CAL
   index = "W1/CAL/W1_ET_COR_CAL";
   UShort_t multFrontE = PreTreatedData->GetFrontEMult();
   for (UShort_t i = 0; i < multFrontE; ++i) {   // loop on front strips (E)
      UShort_t multFrontT = PreTreatedData->GetFrontTMult();
      for (UShort_t j = 0; j < multFrontT; ++j) {   // loop on front strips (T)
         FillSpectra(index,PreTreatedData->GetFrontEStripNbr(i) + 2*fNumberOfStripsFront*(PreTreatedData->GetFrontEDetectorNbr(i)-1), PreTreatedData->GetFrontTStripNbr(i) + 2*fNumberOfStripsFront*(PreTreatedData->GetFrontTDetectorNbr(i)-1));
      } // end loop on front strips (T)
   } // end loop on front strips (E)
   UShort_t multBackE = PreTreatedData->GetBackEMult();
   for (UShort_t i = 0; i < multBackE; ++i) {   // loop on front strips (E)
      UShort_t multBackT = PreTreatedData->GetBackTMult();
      for (UShort_t j = 0; j < multBackT; ++j) {   // loop on front strips (T)
         FillSpectra(index,PreTreatedData->GetBackEStripNbr(i) + (2*(PreTreatedData->GetBackEDetectorNbr(i)-1)+1)*fNumberOfStripsBack, PreTreatedData->GetBackTStripNbr(i) + (2*(PreTreatedData->GetBackTDetectorNbr(i)-1)+1)*fNumberOfStripsBack);
      } // end loop on front strips (T)
   } // end loop on front strips (E)
   

   // W1_HIT_E_CAL
   index = "W1/CAL/W1_HIT_E_CAL";
   for (UShort_t i = 0; i < PreTreatedData->GetFrontEMult(); ++i) {   // loop on front strips
      for (UShort_t j = 0; j < PreTreatedData->GetBackEMult(); ++j) { // loop on back strips
         FillSpectra(index,PreTreatedData->GetFrontEStripNbr(i) + fNumberOfStripsFront*(PreTreatedData->GetFrontEDetectorNbr(i)-1), PreTreatedData->GetBackEStripNbr(i) + fNumberOfStripsBack*(PreTreatedData->GetBackEDetectorNbr(i)-1));
      }
   } // end loop on front strips

   // W1_HIT_T_CAL
   index = "W1/CAL/W1_HIT_T_CAL";
   for (UShort_t i = 0; i < PreTreatedData->GetFrontTMult(); ++i) {   // loop on front strips
      for (UShort_t j = 0; j < PreTreatedData->GetBackTMult(); ++j) { // loop on back strips
         FillSpectra(index,PreTreatedData->GetFrontTStripNbr(i) + fNumberOfStripsFront*(PreTreatedData->GetFrontTDetectorNbr(i)-1), PreTreatedData->GetBackTStripNbr(i) + fNumberOfStripsBack*(PreTreatedData->GetBackTDetectorNbr(i)-1));
      }
   } // end loop on front strips


   // STR_FRONT MULT
   int myMULT[fNumberOfDetectors];
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) myMULT[i] = 0;

   mysize = PreTreatedData->GetFrontEMult(); 
   for (unsigned int i = 0 ; i < mysize ;i++) {
      myMULT[PreTreatedData->GetFrontEDetectorNbr(i)-1] += 1;  
   }

   for (unsigned int i = 0; i < fNumberOfDetectors; i++) {
      index= "W1/CAL/MULT/W1_D"+NPL::itoa(i+1)+"_STR_FRONT_CAL_MULT";
      FillSpectra(index,myMULT[i]);
   }

   // STR_BACK MULT
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) myMULT[i] = 0; 

   mysize = PreTreatedData->GetBackEMult();
   for (unsigned int i = 0 ; i < mysize ;i++) {
      myMULT[PreTreatedData->GetBackEDetectorNbr(i)-1] += 1;  
   }

   for (unsigned int i = 0; i < fNumberOfDetectors; i++) {
      index= "W1/CAL/MULT/W1_D"+NPL::itoa(i+1)+"_STR_BACK_CAL_MULT";
      FillSpectra(index,myMULT[i]);
   }
}



////////////////////////////////////////////////////////////////////////////////
void TW1Spectra::FillPhysicsSpectra(TW1Physics* Physics)
{
   static string index;

   // counters
   index = "W1/PHY/W1_CUTS";
   for (Int_t i = 0; i < fNumberOfCounters; ++i) {   // loop on counters
      FillSpectra(index,i+1, Physics->m_Counter[i]);
   } // end loop on counters


   // W1_ET_COR_PHY
   index = "W1/PHY/W1_ET_COR_PHY";
   for (UShort_t i = 0; i < Physics->GetEventMultiplicity(); ++i) {   // loop on multiplicity
      Int_t stripF = Physics->GetFrontStrip(i) + 2*fNumberOfStripsFront*(Physics->GetDetectorNumber(i)-1);
      if (Physics->GetFrontEnergy(i) > 0 && Physics->GetFrontTime(i) > 0) FillSpectra(index,stripF, stripF);
      Int_t stripB = Physics->GetBackStrip(i)  + (2*(Physics->GetDetectorNumber(i)-1)+1)*fNumberOfStripsBack;
      if (Physics->GetBackEnergy(i) > 0  && Physics->GetBackTime(i) > 0)  FillSpectra(index,stripB, stripB);
   } // end loop on multiplicity 
 

   // W1_HIT_PHY
   index = "W1/PHY/W1_HIT_PHY";
   for (UShort_t i = 0; i < Physics->GetEventMultiplicity(); ++i) {   // loop on multiplicity 
      FillSpectra(index,Physics->GetFrontStrip(i) + fNumberOfStripsFront*(Physics->GetDetectorNumber(i)-1), Physics->GetBackStrip(i) + fNumberOfStripsBack*(Physics->GetDetectorNumber(i)-1));
   } // end loop on multiplicity 


/*   static string index;

   // Kine plot
   unsigned int mysize = Physics->Strip_E.size();
   for(unsigned int i = 0 ; i < mysize ; i++){
      double Theta = Physics->GetPositionOfInteraction(i).Angle(TVector3(0,0,1));
      Theta = Theta/deg;
      double Etot=Physics->Strip_E[i];

      if(Physics->PAD_E[i]>0){
         index = "W1/PHY/W1_PAD_E_E";
         Etot += Physics->PAD_E[i];
         FillSpectra(index,Physics->PAD_E[i],Physics->Strip_E[i]);
         index = "W1/PHY/W1"+NPL::itoa(Physics->DetectorNumber[i])+"_PAD_E_E";
         FillSpectra(index,Physics->PAD_E[i],Physics->Strip_E[i]);

      }
      index = "W1/PHY/W1_THETA_E";
      FillSpectra(index,Theta,Etot);
   }*/
}

