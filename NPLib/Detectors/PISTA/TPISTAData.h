#ifndef __PISTADATA__
#define __PISTADATA__
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
 *  This class hold PISTA Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TPISTAData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // First Stage Front Energy
    vector<unsigned short>  fFirstStage_XE_DetectorNbr;
    vector<unsigned short>  fFirstStage_XE_StripNbr;
    vector<double>          fFirstStage_XE_Energy;
    // First Stage Front Time
    vector<unsigned short>  fFirstStage_XT_DetectorNbr;
    vector<unsigned short>  fFirstStage_XT_StripNbr;
    vector<double>          fFirstStage_XT_Time;
    // First Stage Back Energy
    vector<unsigned short>  fFirstStage_YE_DetectorNbr;
    vector<unsigned short>  fFirstStage_YE_StripNbr;
    vector<double>          fFirstStage_YE_Energy;
    // First Stage Back Time
    vector<unsigned short>  fFirstStage_YT_DetectorNbr;
    vector<unsigned short>  fFirstStage_YT_StripNbr;
    vector<double>          fFirstStage_YT_Time;

    // Second Stage Front Energy
    vector<unsigned short>  fSecondStage_XE_DetectorNbr;
    vector<unsigned short>  fSecondStage_XE_StripNbr;
    vector<double>          fSecondStage_XE_Energy;
    // Second Stage Front Time
    vector<unsigned short>  fSecondStage_XT_DetectorNbr;
    vector<unsigned short>  fSecondStage_XT_StripNbr;
    vector<double>          fSecondStage_XT_Time;
    // Second Stage Back Energy
    vector<unsigned short>  fSecondStage_YE_DetectorNbr;
    vector<unsigned short>  fSecondStage_YE_StripNbr;
    vector<double>          fSecondStage_YE_Energy;
    // Second Stage Back Time
    vector<unsigned short>  fSecondStage_YT_DetectorNbr;
    vector<unsigned short>  fSecondStage_YT_StripNbr;
    vector<double>          fSecondStage_YT_Time;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TPISTAData();
    ~TPISTAData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};
    void Dump() const;


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    // First Stage Energy Front
    inline void SetFirstStageXE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fFirstStage_XE_DetectorNbr.push_back(DetNbr);
      fFirstStage_XE_StripNbr.push_back(StripNbr);
      fFirstStage_XE_Energy.push_back(Energy);
    };//!
    // First Stage Energy Back
    inline void SetFirstStageYE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fFirstStage_YE_DetectorNbr.push_back(DetNbr);
      fFirstStage_YE_StripNbr.push_back(StripNbr);
      fFirstStage_YE_Energy.push_back(Energy);
    };//!
    // First Stage Time Front
    inline void SetFirstStageXT(const UShort_t& DetNbr, const UShort_t StripNbr, const Double_t& Time)	{
      fFirstStage_XT_DetectorNbr.push_back(DetNbr);     
      fFirstStage_XT_StripNbr.push_back(StripNbr);     
      fFirstStage_XT_Time.push_back(Time);
    };//!
    // First Stage Time Back
    inline void SetFirstStageYT(const UShort_t& DetNbr, const UShort_t StripNbr, const Double_t& Time)	{
      fFirstStage_YT_DetectorNbr.push_back(DetNbr);     
      fFirstStage_YT_StripNbr.push_back(StripNbr);     
      fFirstStage_YT_Time.push_back(Time);
    };//!
    
    //////
    // Second Stage Energy Front
    inline void SetSecondStageXE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fSecondStage_XE_DetectorNbr.push_back(DetNbr);
      fSecondStage_XE_StripNbr.push_back(StripNbr);
      fSecondStage_XE_Energy.push_back(Energy);
    };//!
    // Second Stage Energy Back
    inline void SetSecondStageYE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
      fSecondStage_YE_DetectorNbr.push_back(DetNbr);
      fSecondStage_YE_StripNbr.push_back(StripNbr);
      fSecondStage_YE_Energy.push_back(Energy);
    };//!
    // Second Stage Time Front
    inline void SetSecondStageXT(const UShort_t& DetNbr, const UShort_t StripNbr, const Double_t& Time)	{
      fSecondStage_XT_DetectorNbr.push_back(DetNbr);     
      fSecondStage_XT_StripNbr.push_back(StripNbr);     
      fSecondStage_XT_Time.push_back(Time);
    };//!
    // Second Stage Time Back
    inline void SetSecondStageYT(const UShort_t& DetNbr, const UShort_t StripNbr, const Double_t& Time)	{
      fSecondStage_YT_DetectorNbr.push_back(DetNbr);     
      fSecondStage_YT_StripNbr.push_back(StripNbr);     
      fSecondStage_YT_Time.push_back(Time);
    };//!

    //////////////////////    GETTERS    ////////////////////////
    // First Stage Energy X
    inline UShort_t GetFirstStageMultXEnergy() const
      {return fFirstStage_XE_DetectorNbr.size();}
    inline UShort_t GetFirstStage_XE_DetectorNbr(const unsigned int &i) const 
      {return fFirstStage_XE_DetectorNbr[i];}//!
    inline UShort_t GetFirstStage_XE_StripNbr(const unsigned int &i) const 
      {return fFirstStage_XE_StripNbr[i];}//!
    inline Double_t GetFirstStage_XE_Energy(const unsigned int &i) const 
      {return fFirstStage_XE_Energy[i];}//!
    // First Stage Energy Y
    inline UShort_t GetFirstStageMultYEnergy() const
      {return fFirstStage_YE_DetectorNbr.size();}
    inline UShort_t GetFirstStage_YE_DetectorNbr(const unsigned int &i) const 
      {return fFirstStage_YE_DetectorNbr[i];}//!
    inline UShort_t GetFirstStage_YE_StripNbr(const unsigned int &i) const 
      {return fFirstStage_YE_StripNbr[i];}//!
    inline Double_t GetFirstStage_YE_Energy(const unsigned int &i) const 
      {return fFirstStage_YE_Energy[i];}//!
    // First Stage Time X
    inline UShort_t GetFirstStageMultXTime() const
      {return fFirstStage_XT_DetectorNbr.size();}
    inline UShort_t GetFirstStage_XT_DetectorNbr(const unsigned int &i) const 
      {return fFirstStage_XT_DetectorNbr[i];}//!
    inline UShort_t GetFirstStage_XT_StripNbr(const unsigned int &i) const 
      {return fFirstStage_XT_StripNbr[i];}//!
    inline Double_t GetFirstStage_XT_Time(const unsigned int &i) const 
      {return fFirstStage_XT_Time[i];}//!
    // First Stage Time Y
    inline UShort_t GetFirstStageMultYTime() const
      {return fFirstStage_YT_DetectorNbr.size();}
    inline UShort_t GetFirstStage_YT_DetectorNbr(const unsigned int &i) const 
      {return fFirstStage_YT_DetectorNbr[i];}//!
    inline UShort_t GetFirstStage_YT_StripNbr(const unsigned int &i) const 
      {return fFirstStage_YT_StripNbr[i];}//!
    inline Double_t GetFirstStage_YT_Time(const unsigned int &i) const 
      {return fFirstStage_YT_Time[i];}//!
    
    //////
    // Second Stage Energy X
    inline UShort_t GetSecondStageMultXEnergy() const
      {return fSecondStage_XE_DetectorNbr.size();}
    inline UShort_t GetSecondStage_XE_DetectorNbr(const unsigned int &i) const 
      {return fSecondStage_XE_DetectorNbr[i];}//!
    inline UShort_t GetSecondStage_XE_StripNbr(const unsigned int &i) const 
      {return fSecondStage_XE_StripNbr[i];}//!
    inline Double_t GetSecondStage_XE_Energy(const unsigned int &i) const 
      {return fSecondStage_XE_Energy[i];}//!
    // Second Stage Energy Y
    inline UShort_t GetSecondStageMultYEnergy() const
      {return fSecondStage_YE_DetectorNbr.size();}
    inline UShort_t GetSecondStage_YE_DetectorNbr(const unsigned int &i) const 
      {return fSecondStage_YE_DetectorNbr[i];}//!
    inline UShort_t GetSecondStage_YE_StripNbr(const unsigned int &i) const 
      {return fSecondStage_YE_StripNbr[i];}//!
    inline Double_t GetSecondStage_YE_Energy(const unsigned int &i) const 
      {return fSecondStage_YE_Energy[i];}//!
    // Second Stage Time X
    inline UShort_t GetSecondStageMultXTime() const
      {return fSecondStage_XT_DetectorNbr.size();}
    inline UShort_t GetSecondStage_XT_DetectorNbr(const unsigned int &i) const 
      {return fSecondStage_XT_DetectorNbr[i];}//!
    inline UShort_t GetSecondStage_XT_StripNbr(const unsigned int &i) const 
      {return fSecondStage_XT_StripNbr[i];}//!
    inline Double_t GetSecondStage_XT_Time(const unsigned int &i) const 
      {return fSecondStage_XT_Time[i];}//!
    // Second Stage Time Y
    inline UShort_t GetSecondStageMultYTime() const
      {return fSecondStage_YT_DetectorNbr.size();}
    inline UShort_t GetSecondStage_YT_DetectorNbr(const unsigned int &i) const 
      {return fSecondStage_YT_DetectorNbr[i];}//!
    inline UShort_t GetSecondStage_YT_StripNbr(const unsigned int &i) const 
      {return fSecondStage_YT_StripNbr[i];}//!
    inline Double_t GetSecondStage_YT_Time(const unsigned int &i) const 
      {return fSecondStage_YT_Time[i];}//!


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TPISTAData,1)  // PISTAData structure
};

#endif
