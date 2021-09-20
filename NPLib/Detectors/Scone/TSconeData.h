#ifndef __SconeDATA__
#define __SconeDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : March 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Scone Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>

// ROOT
#include "TObject.h"

class TSconeData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into std::vectors in order 
  // to allow multiplicity treatment
  private: 
    // Energy
    std::vector<UShort_t>   fScone_E_DetectorNbr;
    std::vector<UShort_t>   fScone_E_PlasticNbr;
    std::vector<Double_t>   fScone_Energy;

    // Time
    std::vector<UShort_t>   fScone_T_DetectorNbr;
    std::vector<UShort_t>   fScone_T_PlasticNbr;
    std::vector<Double_t>   fScone_Time;

    // Flag for simulation
    std::vector<UShort_t>   fScone_HasCaptured; //=1 for neutron capture on H, 2 on Gd 0 otherwise
    std::vector<double>     fScone_CaptureTime; 
    std::vector<double>     fScone_GammaEnergy;
    std::vector<double>     fScone_ProtonEnergy;
    std::vector<double>     fScone_ProtonTime;
    std::vector<int>        fScone_FCProcess;
    //////////////////////////////////////////////////////////////
    // Constructor and destructor
  public: 
    TSconeData();
    ~TSconeData();


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
    inline void SetEnergy(const UShort_t& DetNbr,const UShort_t& PlasticNbr, const Double_t& Energy){
      fScone_E_DetectorNbr.push_back(DetNbr);
      fScone_E_PlasticNbr.push_back(PlasticNbr);
      fScone_Energy.push_back(Energy);
    };//!

    // Time
    inline void SetTime(const UShort_t& DetNbr,const UShort_t& PlasticNbr, const Double_t& Time)	{
      fScone_T_DetectorNbr.push_back(DetNbr);     
      fScone_T_PlasticNbr.push_back(PlasticNbr);     
      fScone_Time.push_back(Time);
    };//!


    // Flag for simulation
    inline void SetCapture(const UShort_t& capture){
      fScone_HasCaptured.push_back(capture);
    };//
    inline void SetCaptureTime(double capture_time){
      fScone_CaptureTime.push_back(capture_time);
    };//
    inline void SetGammaEnergy(double energy){
      fScone_GammaEnergy.push_back(energy);
    };//
    inline void SetProtonEnergy(double energy){
      fScone_ProtonEnergy.push_back(energy);
    };//
    inline void SetProtonTime(double time){
      fScone_ProtonTime.push_back(time);
    };//
    inline void SetFCProcess(int val){
      fScone_FCProcess.push_back(val);
    };//


    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy() const
    {return fScone_E_DetectorNbr.size();}
    inline UShort_t GetE_DetectorNbr(const unsigned int &i) const 
    {return fScone_E_DetectorNbr[i];}//!
    inline UShort_t GetE_PlasticNbr(const unsigned int &i) const 
    {return fScone_E_PlasticNbr[i];}//!
    inline Double_t Get_Energy(const unsigned int &i) const 
    {return fScone_Energy[i];}//!

    // Time
    inline UShort_t GetMultTime() const
    {return fScone_T_DetectorNbr.size();}
    inline UShort_t GetT_DetectorNbr(const unsigned int &i) const 
    {return fScone_T_DetectorNbr[i];}//!
    inline UShort_t GetT_PlasticNbr(const unsigned int &i) const 
    {return fScone_T_PlasticNbr[i];}//!
    inline Double_t Get_Time(const unsigned int &i) const 
    {return fScone_Time[i];}//!

    // Flag for simulation
    inline UShort_t GetCapture(const unsigned int &i) const
    {return fScone_HasCaptured[i];}//!
    inline double GetCaptureTime(const unsigned int &i) const
    {return fScone_CaptureTime[i];}//!
    inline UShort_t GetGammaMult() const
    {return fScone_GammaEnergy.size();}
    inline double GetGammaEnergy(const unsigned int &i) const
    {return fScone_GammaEnergy[i];}//!
    
    inline UShort_t GetProtonMult() const
    {return fScone_ProtonEnergy.size();}
    inline double GetProtonEnergy(const unsigned int &i) const
    {return fScone_ProtonEnergy[i];}//!
    inline double GetProtonTime(const unsigned int &i) const
    {return fScone_ProtonTime[i];}//!

    inline UShort_t GetFCProcessMult() const
    {return fScone_FCProcess.size();}//!
    inline int GetFCProcess(const unsigned int &i) const
    {return fScone_FCProcess[i];}//!
 
    //////////////////////////////////////////////////////////////
    // Required for ROOT dictionnary
    ClassDef(TSconeData,1)  // SconeData structure
};

#endif
