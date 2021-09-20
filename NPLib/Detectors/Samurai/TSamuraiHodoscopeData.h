#ifndef __SamuraiHodoscopeDATA__
#define __SamuraiHodoscopeDATA__
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : April 2021                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiHodoscope Raw data                                *
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

class TSamuraiHodoscopeData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // UP // 
    // Charge 
    vector<UShort_t>   fSamuraiHodoscope_Qu_ID;
    vector<Double_t>   fSamuraiHodoscope_Qu_Charge;
    
    // Time
    vector<UShort_t>   fSamuraiHodoscope_Tu_ID;
    vector<Double_t>   fSamuraiHodoscope_Tu_Time;
    
    // DOWN // 
    // Charge 
    vector<UShort_t>   fSamuraiHodoscope_Qd_ID;
    vector<Double_t>   fSamuraiHodoscope_Qd_Charge;
    
    // Time
    vector<UShort_t>   fSamuraiHodoscope_Td_ID;
    vector<Double_t>   fSamuraiHodoscope_Td_Time;
 

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSamuraiHodoscopeData();
    ~TSamuraiHodoscopeData();
    

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
    // UP // 
    // Charge
    inline void SetChargeUp(const UShort_t& ID, const Double_t& Charge){
    fSamuraiHodoscope_Qu_ID.push_back(ID);
    fSamuraiHodoscope_Qu_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeUp(const UShort_t& ID, const Double_t& Time){
    fSamuraiHodoscope_Tu_ID.push_back(ID);
    fSamuraiHodoscope_Tu_Time.push_back(Time);
    };//!

    // DOWN // 
    // Charge
    inline void SetChargeDown(const UShort_t& ID, const Double_t& Charge){
    fSamuraiHodoscope_Qd_ID.push_back(ID);
    fSamuraiHodoscope_Qd_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeDown(const UShort_t& ID, const Double_t& Time){
    fSamuraiHodoscope_Td_ID.push_back(ID);
    fSamuraiHodoscope_Td_Time.push_back(Time);
    };//!

    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultQUp() const
      {return fSamuraiHodoscope_Qu_ID.size();}
    inline UShort_t GetQUp_ID(const unsigned int &i) const 
      {return fSamuraiHodoscope_Qu_ID[i];}//!
    inline UShort_t GetQUp_Charge(const unsigned int &i) const 
      {return fSamuraiHodoscope_Qu_Charge[i];}//!

    // Time
    inline UShort_t GetMultTUp() const
      {return fSamuraiHodoscope_Tu_ID.size();}
    inline UShort_t GetTUp_ID(const unsigned int &i) const 
      {return fSamuraiHodoscope_Tu_ID[i];}//!
    inline Double_t GetTUp_Time(const unsigned int &i) const 
      {return fSamuraiHodoscope_Tu_Time[i];}//!
    // Energy
    inline UShort_t GetMultQDown() const
      {return fSamuraiHodoscope_Qd_ID.size();}
    inline UShort_t GetQDown_ID(const unsigned int &i) const 
      {return fSamuraiHodoscope_Qd_ID[i];}//!
    inline Double_t GetQDown_Charge(const unsigned int &i) const 
      {return fSamuraiHodoscope_Qd_Charge[i];}//!

    // Time
    inline UShort_t GetMultTDown() const
      {return fSamuraiHodoscope_Td_ID.size();}
    inline UShort_t GetTDown_ID(const unsigned int &i) const 
      {return fSamuraiHodoscope_Td_ID[i];}//!
    inline Double_t GetTDown_Time(const unsigned int &i) const 
      {return fSamuraiHodoscope_Td_Time[i];}//!

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSamuraiHodoscopeData,1)  // SamuraiHodoscopeData structure
};

#endif
