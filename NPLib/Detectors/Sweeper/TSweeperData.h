#ifndef __SweeperDATA__
#define __SweeperDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: B. Monteagudo  contact address: monteagu@frib.msu.edu                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Sweeper Raw data                                    *
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

class TSweeperData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // Energy (Ion Chamber)
    vector<UShort_t>   fSweeper_E_DetectorNbr;
    vector<Double_t>   fSweeper_Energy;

    // Time (Thin Scint)
    vector<UShort_t>   fSweeper_T_DetectorNbr;
    vector<Double_t>   fSweeper_Time;

    // Position (Drift Chambers)
    vector<UShort_t>   fSweeper_DC_DetectorNbr;
    vector<Double_t>   fSweeper_X;
    vector<Double_t>   fSweeper_Y;
    vector<Double_t>   fSweeper_DriftTime;
    
    /* vector<UShort_t>   fSweeper_CRDC2_DetectorNbr; */
    /* vector<Double_t>   fSweeper_CRDC2_X; */
    /* vector<Double_t>   fSweeper_CRDC2_DriftTime; */
    

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSweeperData();
    ~TSweeperData();
    

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
    // Energy
    inline void SetEnergy(const UShort_t& DetNbr,const Double_t& Energy){
      fSweeper_E_DetectorNbr.push_back(DetNbr);
      fSweeper_Energy.push_back(Energy);
    };//!

    // Time
    inline void SetTime(const UShort_t& DetNbr,const Double_t& Time)	{
      fSweeper_T_DetectorNbr.push_back(DetNbr);     
      fSweeper_Time.push_back(Time);
    };//!

    // Position (CRDCs)
    inline void SetDriftPosition(const UShort_t& DetNbr,const Double_t& DriftTime, const Double_t& X)  {
      fSweeper_DC_DetectorNbr.push_back(DetNbr);     
      fSweeper_X.push_back(X);
      fSweeper_DriftTime.push_back(DriftTime); 
    };//!
    inline void SetPosition(const UShort_t& DetNbr,const Double_t& X, const Double_t& Y)  {
      fSweeper_DC_DetectorNbr.push_back(DetNbr);     
      fSweeper_X.push_back(X);
      fSweeper_Y.push_back(Y); 
    };//!
    
    
    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy() const
      {return fSweeper_E_DetectorNbr.size();}
    inline UShort_t GetE_DetectorNbr(const unsigned int &i) const 
      {return fSweeper_E_DetectorNbr[i];}//!
    inline Double_t Get_Energy(const unsigned int &i) const 
      {return fSweeper_Energy[i];}//!

    // Time
    inline UShort_t GetMultTime() const
      {return fSweeper_T_DetectorNbr.size();}
    inline UShort_t GetT_DetectorNbr(const unsigned int &i) const 
      {return fSweeper_T_DetectorNbr[i];}//!
    inline Double_t Get_Time(const unsigned int &i) const 
      {return fSweeper_Time[i];}//!

    // Position
    inline UShort_t GetMultPosition() const
      {return fSweeper_DC_DetectorNbr.size();}
    inline UShort_t GetDC_DetectorNbr(const unsigned int &i) const 
      {return fSweeper_DC_DetectorNbr[i];}//!
    inline Double_t Get_X(const unsigned int &i) const 
      {return fSweeper_X[i];}//!
    inline Double_t Get_Y(const unsigned int &i) const 
      {return fSweeper_Y[i];}//!
    inline Double_t Get_DriftTime(const unsigned int &i) const 
      {return fSweeper_DriftTime[i];}//!
    
    
  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSweeperData,1)  // SweeperData structure
};

#endif
