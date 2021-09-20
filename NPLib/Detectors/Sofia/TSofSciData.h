#ifndef __SofSciDATA__
#define __SofSciDATA__
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
 *  This class hold SofSci Raw data                                    *
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

class TSofSciData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int> fSofSci_DetNbr;
    vector<int> fSofSci_Pmt;
    vector<int> fSofSci_CT;
    vector<int> fSofSci_FT;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofSciData();
    ~TSofSciData();
    

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
    inline void SetDetectorNbr(int det){fSofSci_DetNbr.push_back(det);};//!
    inline void SetPmt(int pmt){fSofSci_Pmt.push_back(pmt);};//!
    inline void SetCoarseTime(int Time){fSofSci_CT.push_back(Time);};//!
    inline void SetFineTime(int Time){fSofSci_FT.push_back(Time);};//!

    //////////////////////    GETTERS    ////////////////////////
    inline int GetMultiplicity() const {return fSofSci_DetNbr.size();}//!
    inline int GetDetectorNbr(const unsigned int &i) const {return fSofSci_DetNbr[i];}//! 
    inline int GetPmt(const unsigned int &i) const {return fSofSci_Pmt[i];}//! 
    inline int GetCoarseTime(const unsigned int &i) const {return fSofSci_CT[i];}//!     
    inline int GetFineTime(const unsigned int &i) const {return fSofSci_FT[i];}//!     

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofSciData,1)  // SofSciData structure
};

#endif
