#ifndef TMinosDATA_H
#define TMinosDATA_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin                                        *
 * Maintainer : Adrien Matta                                                 *
 * contact address:  matta@lpccaen.in2p3.fr                                  *
 *                                                                           *
 *                                                                           *
 * Creation Date  : October 2018                                             *
 * Last update    : April 2021                                               *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Raw data                                           *
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

using namespace std;

class TMinosData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
  
    // Pads
    vector<unsigned short>  fMinos_PadNumber;
    vector<unsigned short>  fMinos_HitCount;
    vector< vector<unsigned short> >   fMinos_Charge;
    vector< vector<unsigned short> >   fMinos_Time;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TMinosData();
    ~TMinosData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};

  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    inline void SetPad(const unsigned short& Pad, const unsigned short& HitCount, const vector<int>* time , const vector<int>* charge){
      fMinos_HitCount.push_back(HitCount);
      fMinos_PadNumber.push_back(Pad);
      unsigned int size=time->size();
      static vector<unsigned short> Charge,Time; 
      Charge.clear();Time.clear();
      for(unsigned int i = 0 ; i < size ; i++){
        Time.push_back((*time)[i]);    
        Charge.push_back((*charge)[i]);    
      }
      fMinos_Time.push_back(Time);    
      fMinos_Charge.push_back(Charge);    
    }//!

    //////////////////////    GETTERS    ////////////////////////

      inline unsigned int GetPadMult() const
        {return fMinos_PadNumber.size() ;}//!
      inline unsigned short GetPadNumber(const unsigned int& i) const
        {return fMinos_PadNumber[i] ;}//!
      inline vector<unsigned short> GetCharge(const unsigned int& i)const
        {return fMinos_Charge[i] ;}//!
      inline vector<unsigned short> GetTime(const unsigned int& i)const 
        {return fMinos_Time[i] ;}//!
      inline vector<unsigned short>* GetChargePtr(const unsigned int& i)
        {return &fMinos_Charge[i];}//!
      inline vector<unsigned short>* GetTimePtr(const unsigned int& i)
        {return &fMinos_Time[i];}//!


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TMinosData,2)  // MinosData structure
};

#endif
