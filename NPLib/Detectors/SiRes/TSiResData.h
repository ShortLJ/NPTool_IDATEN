#ifndef __SIRESDATA__
#define __SIRESDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author:    contact address:                                      *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include <vector>
#include <iostream>

#include "TObject.h"


class TSiResData : public TObject {
 private:
   // Energy 
   std::vector<int>     fSiRes_E_Number;
   std::vector<int>     fSiRes_E_Channel;
   std::vector<double>  fSiRes_E_Energy;   
   // Time
   std::vector<int>     fSiRes_T_Number;
   std::vector<double>  fSiRes_T_Time;
   std::vector<int>     fSiRes_E_EnergyBack_Number;
   std::vector<double>  fSiRes_E_EnergyBack;
   
 public:
   TSiResData();
   virtual ~TSiResData();

   void   Clear();
   void   Clear(const Option_t*) {};
   void   Dump() const;

   /////////////////////           GETTERS           ////////////////////////
   // Energy
   unsigned int   GetEnergyMult()   {return fSiRes_E_Number.size();}
   int            GetEChannelNumber(int i) {return fSiRes_E_Channel[i];}
   int            GetEDetectorNumber(int i) {return fSiRes_E_Number[i];}
   double         GetEEnergy(int i) {return fSiRes_E_Energy[i];}
   // Time 
   unsigned int   GetTimeMult()     {return fSiRes_T_Number.size();}
   int            GetTDetectorNumber(int i) {return fSiRes_T_Number[i];}
   double         GetTTime(int i)   {return fSiRes_T_Time[i];}
   double         GetEEnergyBack(int i)   	{return fSiRes_E_EnergyBack[i];}
   int            GetEEnergyBackDetectorNumber(int i)   {return fSiRes_E_EnergyBack_Number[i];}
   double         GetEEnergyBackMult()   	{return fSiRes_E_EnergyBack.size();}

   /////////////////////           SETTERS           ////////////////////////
   // Energy
   void     SetEDetectorNumber(int N)    	{fSiRes_E_Number.push_back(N);}
   void     SetEChannelNumber(int channel)    	{fSiRes_E_Channel.push_back(channel);}
   void     SetEEnergy(double E) 		{fSiRes_E_Energy.push_back(E);}
   // time
   void     SetTDetectorNumber(int N)    	{fSiRes_T_Number.push_back(N);}
   void     SetTTime(double T)   		{fSiRes_T_Time.push_back(T);}
   void     SetEEnergyBack(double E)   		{fSiRes_E_EnergyBack.push_back(E);}
   void     SetEEnergyBackDetectorNumber(int N) {fSiRes_E_EnergyBack_Number.push_back(N);}

   ClassDef(TSiResData,1)  // SiResData structure
};

#endif
