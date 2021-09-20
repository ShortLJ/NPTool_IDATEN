#ifndef __KhalaDATA__
#define __KhalaDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: M. Labiche    contact address: marc.labiche@stfc.ac.uk   *
 *                                                                           *
 * Creation Date  : 04/12/2009                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class described the raw data of the Khala detector                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <vector>

#include "TObject.h"
using namespace std ;


class TKhalaData : public TObject {

  protected:
    // LaBr3
    // Energy
    vector<UShort_t>  fKHALA_LaBr3_E_DetectorNbr;
    vector<Double_t>  fKHALA_LaBr3_E_Energy;
    // Time
    vector<UShort_t>  fKHALA_LaBr3_T_DetectorNbr;
    vector<Double_t>  fKHALA_LaBr3_T_Time;

  public:
    TKhalaData();
    virtual ~TKhalaData();

    void   Clear();
    void   Clear(const Option_t*) {};
    void   Dump() const;

    /////////////////////           GETTERS           ////////////////////////
    // (E)
    UShort_t GetKhalaLaBr3EMult()               {
      return fKHALA_LaBr3_E_DetectorNbr.size();
    }
    
    UShort_t GetKhalaLaBr3EDetectorNbr(Int_t i)    {
      return fKHALA_LaBr3_E_DetectorNbr.at(i);
    }
    Double_t GetKhalaLaBr3EEnergy(Int_t i)      {
      return fKHALA_LaBr3_E_Energy.at(i);
    }


    // (T)
    UShort_t GetKhalaLaBr3TMult()               {
      return fKHALA_LaBr3_E_DetectorNbr.size();
    }
    
    UShort_t GetKhalaLaBr3TDetectorNbr(Int_t i)    {
      return fKHALA_LaBr3_T_DetectorNbr.at(i);
    }
    Double_t GetKhalaLaBr3TTime(Int_t i)      {
      return fKHALA_LaBr3_T_Time.at(i);
    }

    /////////////////////           SETTERS           ////////////////////////
    // (E)
    void SetKhalaLaBr3E(UShort_t DetectorNbr, Double_t E){
      fKHALA_LaBr3_E_DetectorNbr.push_back(DetectorNbr); 
      fKHALA_LaBr3_E_Energy.push_back(E); 
    }

    void SetKhalaLaBr3EDetectorNbr(UShort_t DetectorNbr)    {
      fKHALA_LaBr3_E_DetectorNbr.push_back(DetectorNbr);
    }
    void SetKhalaLaBr3EEnergy(Double_t Energy)      {
      fKHALA_LaBr3_E_Energy.push_back(Energy);
    }

    // (T)
    void SetKhalaLaBr3T(UShort_t DetectorNbr, Double_t T){
      fKHALA_LaBr3_T_DetectorNbr.push_back(DetectorNbr); 
      fKHALA_LaBr3_T_Time.push_back(T); 
    }

    void SetKhalaLaBr3TDetectorNbr(UShort_t DetectorNbr)    {
      fKHALA_LaBr3_T_DetectorNbr.push_back(DetectorNbr);
    }
    void SetKhalaLaBr3TTime(Double_t Time)      {
      fKHALA_LaBr3_T_Time.push_back(Time);
    }

    ClassDef(TKhalaData,1)  // KhalaData structure
};

#endif
