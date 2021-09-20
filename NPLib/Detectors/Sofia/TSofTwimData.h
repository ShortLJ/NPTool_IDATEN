#ifndef __SofTwimDATA__
#define __SofTwimDATA__
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
 *  This class hold SofTwim Raw data                                    *
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

class TSofTwimData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int>      fTwim_SectionNbr;
    vector<int>      fTwim_AnodeNbr;
    vector<double>   fTwim_Energy;
    vector<double>   fTwim_DriftTime;
    vector<bool>     fTwim_PileUp;
    vector<bool>     fTwim_Overflow;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofTwimData();
    ~TSofTwimData();
    

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
    inline void SetSectionNbr(int sec){fTwim_SectionNbr.push_back(sec);};//!
    inline void SetAnodeNbr(int det){fTwim_AnodeNbr.push_back(det);};//!
    inline void SetEnergy(double Energy){fTwim_Energy.push_back(Energy);};//!
    inline void SetDriftTime(double Time){fTwim_DriftTime.push_back(Time);};//!
    inline void SetPileUp(bool ispileup){fTwim_PileUp.push_back(ispileup);};//!
    inline void SetOverflow(bool isoverflow){fTwim_Overflow.push_back(isoverflow);};//!

    //////////////////////    GETTERS    ////////////////////////
    inline int GetMultiplicity() const {return fTwim_AnodeNbr.size();}//!
    inline int GetSectionNbr(const unsigned int &i) const {return fTwim_SectionNbr[i];}//! 
    inline int GetAnodeNbr(const unsigned int &i) const {return fTwim_AnodeNbr[i];}//! 
    inline double GetEnergy(const unsigned int &i) const {return fTwim_Energy[i];}//!     
    inline double GetDriftTime(const unsigned int &i) const {return fTwim_DriftTime[i];}//!     
    inline bool GetPileUp(const unsigned int &i) const {return fTwim_PileUp[i];}//!     
    inline bool GetOverflow(const unsigned int &i) const {return fTwim_Overflow[i];}//!     

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofTwimData,1)  // SofTwimData structure
};

#endif
