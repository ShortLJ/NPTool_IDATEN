#ifndef __FissionFragmentDATA__
#define __FissionFragmentDATA__
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
 *  This class hold FissionFragment Raw data                                    *
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

class TSofFissionFragment : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private:
    vector<double> fFF_Z; 
    vector<double> fFF_Qmax;
    vector<double> fFF_AoQ;
    vector<double> fFF_A;
    vector<double> fFF_Beta;
    vector<double> fFF_TOF;
    vector<double> fFF_Gamma;
    vector<double> fFF_Brho;
    double fFF_Zsum;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofFissionFragment();
    ~TSofFissionFragment();
    

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
    inline void SetZsum(double val){fFF_Zsum = val;};//!
    inline void SetZ(double val){fFF_Z.push_back(val);};//!
    inline void SetQmax(double val){fFF_Qmax.push_back(val);};//!
    inline void SetAoQ(double val){fFF_AoQ.push_back(val);};//!
    inline void SetA(double val){fFF_A.push_back(val);};//!
    inline void SetBeta(double val){fFF_Beta.push_back(val);};//!
    inline void SetTOF(double val){fFF_TOF.push_back(val);};//!
    inline void SetGamma(double val){fFF_Gamma.push_back(val);};//!
    inline void SetBrho(double val){fFF_Brho.push_back(val);};//!

    //////////////////////    GETTERS    ////////////////////////
    int GetMult() {return fFF_Z.size();}//!
    inline double GetZsum() const {return fFF_Zsum;}//! 
    inline double GetZ(int i) const {return fFF_Z[i];}//! 
    inline double GetQmax(int i) const {return fFF_Qmax[i];}//! 
    inline double GetAoQ(int i) const {return fFF_AoQ[i];}//! 
    inline double GetA(int i) const {return fFF_A[i];}//! 
    inline double GetBeta(int i) const {return fFF_Beta[i];}//! 
    inline double GetTOF(int i) const {return fFF_TOF[i];}//! 
    inline double GetGamma(int i) const {return fFF_Gamma[i];}//! 
    inline double GetBrho(int i) const {return fFF_Brho[i];}//! 

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofFissionFragment,1)  // FissionFragment structure
};

#endif
