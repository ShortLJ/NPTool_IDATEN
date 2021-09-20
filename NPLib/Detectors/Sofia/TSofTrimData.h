#ifndef __SofTrimDATA__
#define __SofTrimDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SofTrim Raw data                                    *
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

class TSofTrimData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int>      fTrim_SectionNbr;
    vector<int>      fTrim_AnodeNbr;
    vector<double>   fTrim_Energy;
    vector<double>   fTrim_DriftTime;
    vector<bool>     fTrim_PileUp;
    vector<bool>     fTrim_Overflow;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofTrimData();
    ~TSofTrimData();
    

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
    inline void SetSectionNbr(int sec){fTrim_SectionNbr.push_back(sec);};//!
    inline void SetAnodeNbr(int det){fTrim_AnodeNbr.push_back(det);};//!
    inline void SetEnergy(double Energy){fTrim_Energy.push_back(Energy);};//!
    inline void SetDriftTime(double Time){fTrim_DriftTime.push_back(Time);};//!
    inline void SetPileUp(bool ispileup){fTrim_PileUp.push_back(ispileup);};//!
    inline void SetOverflow(bool isoverflow){fTrim_Overflow.push_back(isoverflow);};//!

    //////////////////////    GETTERS    ////////////////////////
    inline int GetMultiplicity() const {return fTrim_AnodeNbr.size();}//!
    inline int GetSectionNbr(const unsigned int &i) const {return fTrim_SectionNbr[i];}//! 
    inline int GetAnodeNbr(const unsigned int &i) const {return fTrim_AnodeNbr[i];}//! 
    inline double GetEnergy(const unsigned int &i) const {return fTrim_Energy[i];}//!     
    inline double GetDriftTime(const unsigned int &i) const {return fTrim_DriftTime[i];}//!     
    inline bool GetPileUp(const unsigned int &i) const {return fTrim_PileUp[i];}//!     
    inline bool GetOverflow(const unsigned int &i) const {return fTrim_Overflow[i];}//!     

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofTrimData,1)  // SofTrimData structure
};

#endif
