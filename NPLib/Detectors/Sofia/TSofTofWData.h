#ifndef __SofTofWDATA__
#define __SofTofWDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : November 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SofTofW Raw data                                    *
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

class TSofTofWData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int>      fTOF_PlasticNbr;
    vector<int>      fTOF_Pmt;
    vector<double>   fTOF_Energy;
    vector<double>   fTOF_CT;
    vector<double>   fTOF_FT;
    vector<bool>     fTOF_WhichFlag;


    //////////////////////////////////////////////////////////////
    // Constructor and destructor
  public: 
    TSofTofWData();
    ~TSofTofWData();


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
    inline void SetPlasticNbr(int plastic){fTOF_PlasticNbr.push_back(plastic);};//!
    inline void SetPmt(int pmt){fTOF_Pmt.push_back(pmt);};//!
    inline void SetEnergy(double Energy){fTOF_Energy.push_back(Energy);};//!
    inline void SetCoarseTime(double Time){fTOF_CT.push_back(Time);};//!
    inline void SetFineTime(double Time){fTOF_FT.push_back(Time);};//!
    inline void SetWhichFlag(bool flag){fTOF_WhichFlag.push_back(flag);};//!

    //////////////////////    GETTERS    ////////////////////////
    inline int GetMultiplicity() const {return fTOF_PlasticNbr.size();}//!
    inline int GetPlasticNbr(const unsigned int &i) const {return fTOF_PlasticNbr[i];}//!     
    inline int GetPmt(const unsigned int &i) const {return fTOF_Pmt[i];}//!     
    inline double GetEnergy(const unsigned int &i) const {return fTOF_Energy[i];}//!     
    inline double GetCoarseTime(const unsigned int &i) const {return fTOF_CT[i];}//!     
    inline double GetFineTime(const unsigned int &i) const {return fTOF_FT[i];}//!     
    inline bool GetWhichFlag(const unsigned int &i) const {return fTOF_WhichFlag[i];}//!     

    //////////////////////////////////////////////////////////////
    // Required for ROOT dictionnary
    ClassDef(TSofTofWData,1)  // SofTofWData structure
};

#endif
