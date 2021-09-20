#ifndef __SofBeamIDDATA__
#define __SofBeamIDDATA__
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
 *  This class hold SofBeamID Raw data                                    *
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

class TSofBeamID : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private:
    double fBeam_Z; 
    double fBeam_Qmax;
    double fBeam_AoQ;
    double fBeam_A;
    double fBeam_Beta;
    double fBeam_Gamma;
    double fBeam_Brho;
    double fBeam_XS2;
    double fBeam_XCC;
    double fBeam_YCC;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofBeamID();
    ~TSofBeamID();
    

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
    inline void SetZbeam(double val){fBeam_Z = val;};//!
    inline void SetQmax(double val){fBeam_Qmax = val;};//!
    inline void SetAoQ(double val){fBeam_AoQ = val;};//!
    inline void SetAbeam(double val){fBeam_A = val;};//!
    inline void SetBeta(double val){fBeam_Beta = val;};//!
    inline void SetGamma(double val){fBeam_Gamma = val;};//!
    inline void SetBrho(double val){fBeam_Brho = val;};//!
    inline void SetXS2(double val){fBeam_XS2 = val;};//!
    inline void SetXCC(double val){fBeam_XCC = val;};//!
    inline void SetYCC(double val){fBeam_YCC = val;};//!

    //////////////////////    GETTERS    ////////////////////////
    inline double GetZbeam() const {return fBeam_Z;}//! 
    inline double GetQmax() const {return fBeam_Qmax;}//! 
    inline double GetAoQ() const {return fBeam_AoQ;}//! 
    inline double GetAbeam() const {return fBeam_A;}//! 
    inline double GetBeta() const {return fBeam_Beta;}//! 
    inline double GetGamma() const {return fBeam_Gamma;}//! 
    inline double GetBrho() const {return fBeam_Brho;}//! 
    inline double GetXS2() const {return fBeam_XS2;}//! 
    inline double GetXCC() const {return fBeam_XCC;}//! 
    inline double GetYCC() const {return fBeam_YCC;}//! 

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofBeamID,1)  // SofBeamID structure
};

#endif
