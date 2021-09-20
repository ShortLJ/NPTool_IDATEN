#ifndef __NebulaDATA__
#define __NebulaDATA__
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Nebula Raw data                                          *
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

class TNebulaData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // UP // 
    // Charge 
    vector<UShort_t>   fNebula_Qu_ID;
    vector<Double_t>   fNebula_Qu_Charge;
    
    // Time
    vector<UShort_t>   fNebula_Tu_ID;
    vector<Double_t>   fNebula_Tu_Time;
    
    // DOWN // 
    // Charge 
    vector<UShort_t>   fNebula_Qd_ID;
    vector<Double_t>   fNebula_Qd_Charge;
    
    // Time
    vector<UShort_t>   fNebula_Td_ID;
    vector<Double_t>   fNebula_Td_Time;
 

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TNebulaData();
    ~TNebulaData();
    

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
    inline void SetChargeUp(const Double_t& ID, const Double_t& Charge){
    fNebula_Qu_ID.push_back(ID);
    fNebula_Qu_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeUp(const Double_t& ID, const Double_t& Time){
    fNebula_Tu_ID.push_back(ID);
    fNebula_Tu_Time.push_back(Time);
    };//!

    // DOWN // 
    // Charge
    inline void SetChargeDown(const Double_t& ID, const Double_t& Charge){
    fNebula_Qd_ID.push_back(ID);
    fNebula_Qd_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeDown(const Double_t& ID, const Double_t& Time){
    fNebula_Td_ID.push_back(ID);
    fNebula_Td_Time.push_back(Time);
    };//!

    //////////////////////    GETTERS    ////////////////////////
    // MULT //
    // Charge 
    inline unsigned int GetChargeUpMult() const
      {return fNebula_Qu_ID.size();};
    // Time
    inline unsigned int GetTimeUpMult() const
      {return fNebula_Tu_ID.size();};
    // Charge
    inline unsigned int GetChargeDownMult() const
      {return fNebula_Qd_ID.size();};
    // Time
    inline unsigned int GetTimeDownMult() const
      {return fNebula_Td_ID.size();};

    // Value // 
    // Charge 
    inline UShort_t GetChargeUpID(unsigned int& i) const
      {return fNebula_Qu_ID[i];};
    inline double GetChargeUp(unsigned int& i) const
      {return fNebula_Qu_Charge[i];};
    // Time 
    inline UShort_t GetTimeUpID(unsigned int& i) const
      {return fNebula_Tu_ID[i];};
    inline double GetTimeUp(unsigned int& i) const
      {return fNebula_Tu_Time[i];};
    // Charge 
    inline UShort_t GetChargeDownID(unsigned int& i) const
      {return fNebula_Qd_ID[i];};
    inline double GetChargeDown(unsigned int& i) const
      {return fNebula_Qd_Charge[i];};
    // Time 
    inline UShort_t GetTimeDownID(unsigned int& i) const
      {return fNebula_Td_ID[i];};
    inline double GetTimeDown(unsigned int& i) const
      {return fNebula_Td_Time[i];};


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TNebulaData,1)  // NebulaData structure
};

#endif
