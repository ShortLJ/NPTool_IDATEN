#ifndef __StrasseDATA__
#define __StrasseDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny    contact : flavigny@lpccaen.in2p3.fr       *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Strasse Raw data                                         *
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

class TStrasseData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  // private: 
  public: 
    // First Stage Front Energy
    vector<unsigned short>  fInner_TE_DetectorNbr;
    vector<unsigned short>  fInner_TE_StripNbr;
    vector<double>          fInner_TE_Energy;
    // First Stage Back Energy
    vector<unsigned short>  fInner_LE_DetectorNbr;
    vector<unsigned short>  fInner_LE_StripNbr;
    vector<double>          fInner_LE_Energy;

    // Second Stage Front Energy
    vector<unsigned short>  fOuter_TE_DetectorNbr;
    vector<unsigned short>  fOuter_TE_StripNbr;
    vector<double>          fOuter_TE_Energy;
    // Second Stage Back Energy
    vector<unsigned short>  fOuter_LE_DetectorNbr;
    vector<unsigned short>  fOuter_LE_StripNbr;
    vector<double>          fOuter_LE_Energy;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TStrasseData();
    ~TStrasseData();
    

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
    inline void SetInnerTE(const unsigned short& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fInner_TE_DetectorNbr.push_back(DetNbr);
      fInner_TE_StripNbr.push_back(StripNbr);
      fInner_TE_Energy.push_back(Energy);
    };//!
    // First Stage Energy Back
    inline void SetInnerLE(const unsigned short& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fInner_LE_DetectorNbr.push_back(DetNbr);
      fInner_LE_StripNbr.push_back(StripNbr);
      fInner_LE_Energy.push_back(Energy);
    };//!
   
    //////
    // Second Stage Energy Front
    inline void SetOuterTE(const unsigned short& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fOuter_TE_DetectorNbr.push_back(DetNbr);
      fOuter_TE_StripNbr.push_back(StripNbr);
      fOuter_TE_Energy.push_back(Energy);
    };//!
    // Second Stage Energy Back
    inline void SetOuterLE(const unsigned short& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fOuter_LE_DetectorNbr.push_back(DetNbr);
      fOuter_LE_StripNbr.push_back(StripNbr);
      fOuter_LE_Energy.push_back(Energy);
    };//!

    //////////////////////    GETTERS    ////////////////////////
    // First Stage Energy T
    inline unsigned short GetInnerMultTEnergy() const
      {return fInner_TE_DetectorNbr.size();}
    inline unsigned short GetInner_TE_DetectorNbr(const unsigned int &i) const 
      {return fInner_TE_DetectorNbr[i];}//!
    inline unsigned short GetInner_TE_StripNbr(const unsigned int &i) const 
      {return fInner_TE_StripNbr[i];}//!
    inline Double_t GetInner_TE_Energy(const unsigned int &i) const 
      {return fInner_TE_Energy[i];}//!
    // First Stage Energy L
    inline unsigned short GetInnerMultLEnergy() const
      {return fInner_LE_DetectorNbr.size();}
    inline unsigned short GetInner_LE_DetectorNbr(const unsigned int &i) const 
      {return fInner_LE_DetectorNbr[i];}//!
    inline unsigned short GetInner_LE_StripNbr(const unsigned int &i) const 
      {return fInner_LE_StripNbr[i];}//!
    inline Double_t GetInner_LE_Energy(const unsigned int &i) const 
      {return fInner_LE_Energy[i];}//!
   
    //////
    // Second Stage Energy T
    inline unsigned short GetOuterMultTEnergy() const
      {return fOuter_TE_DetectorNbr.size();}
    inline unsigned short GetOuter_TE_DetectorNbr(const unsigned int &i) const 
      {return fOuter_TE_DetectorNbr[i];}//!
    inline unsigned short GetOuter_TE_StripNbr(const unsigned int &i) const 
      {return fOuter_TE_StripNbr[i];}//!
    inline Double_t GetOuter_TE_Energy(const unsigned int &i) const 
      {return fOuter_TE_Energy[i];}//!
    // Second Stage Energy L
    inline unsigned short GetOuterMultLEnergy() const
      {return fOuter_LE_DetectorNbr.size();}
    inline unsigned short GetOuter_LE_DetectorNbr(const unsigned int &i) const 
      {return fOuter_LE_DetectorNbr[i];}//!
    inline unsigned short GetOuter_LE_StripNbr(const unsigned int &i) const 
      {return fOuter_LE_StripNbr[i];}//!
    inline Double_t GetOuter_LE_Energy(const unsigned int &i) const 
      {return fOuter_LE_Energy[i];}//!

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TStrasseData,1)  // StrasseData structure
};

#endif
