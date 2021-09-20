#ifndef __SofMwpcDATA__
#define __SofMwpcDATA__
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
 *  This class hold SofMwpc Raw data                                    *
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

class TSofMwpcData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int> fMwpc_DetNbr;
    vector<int> fMwpc_Plane;
    vector<int> fMwpc_Pad;
    vector<int> fMwpc_Charge;

    //////////////////////////////////////////////////////////////
    // Constructor and destructor
  public: 
    TSofMwpcData();
    ~TSofMwpcData();


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
    inline void SetDetectorNbr(int det){fMwpc_DetNbr.push_back(det);};//!
    inline void SetPlane(int p){fMwpc_Plane.push_back(p);};//!
    inline void SetPad(int p){fMwpc_Pad.push_back(p);};//!
    inline void SetCharge(int q){fMwpc_Charge.push_back(q);};//!

    //////////////////////    GETTERS    ////////////////////////
    inline int GetMultiplicity() const {return fMwpc_DetNbr.size();}//!
    inline int GetDetectorNbr(const unsigned int &i) const {return fMwpc_DetNbr[i];}//! 
    inline int GetPlane(const unsigned int &i) const {return fMwpc_Plane[i];}//! 
    inline int GetPad(const unsigned int &i) const {return fMwpc_Pad[i];}//! 
    inline int GetCharge(const unsigned int &i) const {return fMwpc_Charge[i];}//! 

    //////////////////////////////////////////////////////////////
    // Required for ROOT dictionnary
    ClassDef(TSofMwpcData,1)  // SofMwpcData structure
};

#endif
