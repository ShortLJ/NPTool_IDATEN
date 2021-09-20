#ifndef __SofAtDATA__
#define __SofAtDATA__
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
 *  This class hold SofAt Raw data                                    *
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

class TSofAtData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int>      fAT_AnodeNbr;
    vector<double>   fAT_Energy;
    vector<double>   fAT_Time;
    vector<bool>     fAT_PileUp;
    vector<bool>     fAT_Overflow;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofAtData();
    ~TSofAtData();
    

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
    inline void SetAnodeNbr(int det){fAT_AnodeNbr.push_back(det);};//!
    inline void SetEnergy(double Energy){fAT_Energy.push_back(Energy);};//!
    inline void SetTime(double Time){fAT_Time.push_back(Time);};//!
    inline void SetPileUp(bool ispileup){fAT_PileUp.push_back(ispileup);};//!
    inline void SetOverflow(bool isoverflow){fAT_Overflow.push_back(isoverflow);};//!

    //////////////////////    GETTERS    ////////////////////////
    inline int GetMultiplicity() const {return fAT_AnodeNbr.size();}//!
    inline int GetAnodeNbr(const unsigned int &i) const {return fAT_AnodeNbr[i];}//! 
    inline double GetEnergy(const unsigned int &i) const {return fAT_Energy[i];}//!     
    inline double GetTime(const unsigned int &i) const {return fAT_Time[i];}//!     
    inline bool GetPileUp(const unsigned int &i) const {return fAT_PileUp[i];}//!     
    inline bool GetOverflow(const unsigned int &i) const {return fAT_Overflow[i];}//!     

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofAtData,1)  // SofAtData structure
};

#endif
