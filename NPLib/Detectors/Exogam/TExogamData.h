#ifndef __EXOGAMDATA__
#define __EXOGAMDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : march 2009                                               *
 * Last update    : july 2019                                                         *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Exogam Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: Added vectors for real energy/time values (double) (T.Goigoux CEA)                                                                 *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/


// STL
#include <vector>
using namespace std;
#include "TObject.h"



class TExogamData : public TObject {
 private:
  // real energy value (for npsimulation)
  vector<UShort_t>	fEXO_Clover;
  vector<UShort_t>	fEXO_Cristal;
  vector<double> 	fEXO_Energy;
  vector<double> 	fEXO_Time;
  // ECC / Energy
  vector<UShort_t>	fEXO_ECC_E_Clover;
  vector<UShort_t>	fEXO_ECC_E_Cristal;
  vector<UShort_t>	fEXO_ECC_E_Energy;
  // ECC / Time
  vector<UShort_t>	fEXO_ECC_T_Clover;
  vector<UShort_t>	fEXO_ECC_T_Cristal;
  vector<UShort_t>	fEXO_ECC_T_Time;
  // GOCCE / Energy
  vector<UShort_t>	fEXO_GOCCE_E_Clover;
  vector<UShort_t>	fEXO_GOCCE_E_Cristal;
  vector<UShort_t>	fEXO_GOCCE_E_Segment;
  vector<UShort_t>	fEXO_GOCCE_E_Energy;
  // GOCCE / Time
  vector<UShort_t>	fEXO_GOCCE_T_Clover;
  vector<UShort_t>	fEXO_GOCCE_T_Cristal;
  vector<UShort_t>	fEXO_GOCCE_T_Segment;
  vector<UShort_t>	fEXO_GOCCE_T_Time;
  // GeFill
  UShort_t             fEXO_Fill;
 public:
  TExogamData();
  virtual ~TExogamData();

  void Clear();
  void Clear(const Option_t*) {};
  void Dump() const;


  /////////////////////           SETTERS           ////////////////////////
    void	SetClover(UShort_t clo)	{ fEXO_Clover.push_back(clo);}
    void	SetCristal(UShort_t cris)	{ fEXO_Cristal.push_back(cris);}
    void	SetEnergy(double ener)	{ fEXO_Energy.push_back(ener);}
    void	SetTime(double time)	{ fEXO_Time.push_back(time);}
   // ECC / Energy
    void	SetECCEClover(UShort_t clov)	{ fEXO_ECC_E_Clover.push_back(clov);}
    void	SetECCECristal(UShort_t cris)	{ fEXO_ECC_E_Cristal.push_back(cris);}
    void	SetECCEEnergy(UShort_t ener)	{ fEXO_ECC_E_Energy.push_back(ener);}
    // ECC / Time
    void	SetECCTClover(UShort_t clov)	{ fEXO_ECC_T_Clover.push_back(clov);}
    void	SetECCTCristal(UShort_t cris)	{ fEXO_ECC_T_Cristal.push_back(cris);}
    void	SetECCTTime(UShort_t time)	{ fEXO_ECC_T_Time.push_back(time);}
    // GOCCE / Energy
    void	SetGOCCEEClover(UShort_t clov)	{ fEXO_GOCCE_E_Clover.push_back(clov);}
    void	SetGOCCEECristal(UShort_t cris)	{ fEXO_GOCCE_E_Cristal.push_back(cris);}
    void	SetGOCCEESegment(UShort_t seg)	{ fEXO_GOCCE_E_Segment.push_back(seg);}
    void	SetGOCCEEEnergy(UShort_t ener)	{ fEXO_GOCCE_E_Energy.push_back(ener);}
    // GOCCE / Time
    void	SetGOCCETClover(UShort_t clov)	{ fEXO_GOCCE_T_Clover.push_back(clov);}
    void	SetGOCCETCristal(UShort_t cris)	{ fEXO_GOCCE_T_Cristal.push_back(cris);}
    void	SetGOCCETSegment(UShort_t seg)	{ fEXO_GOCCE_T_Segment.push_back(seg);}
    void	SetGOCCETTime(UShort_t time)	{ fEXO_GOCCE_T_Time.push_back(time);}
    //GeFill
    void SetGeFill(UShort_t Fill)        {fEXO_Fill = Fill;}

    /////////////////////           GETTERS           ////////////////////////
      UShort_t	GetClover(Int_t i)	{return fEXO_Clover[i];}
      UShort_t	GetCristal(Int_t i)	{return fEXO_Cristal[i];}
      UShort_t	GetEnergy(Int_t i)	{return fEXO_Energy[i];}
      UShort_t	GetTime(Int_t i)	{return fEXO_Time[i];}
      // ECC / Energy
      // UShort_t	GetCloverMult()		{return fEXO_ECC_E_Clover.size();}       
      UShort_t	GetECCEMult()		{return fEXO_ECC_E_Clover.size();}             
      UShort_t	GetECCEClover(Int_t i)	{return fEXO_ECC_E_Clover[i];}
      UShort_t	GetECCECristal(Int_t i)	{return fEXO_ECC_E_Cristal[i];}
      UShort_t	GetECCEEnergy(Int_t i)	{return fEXO_ECC_E_Energy[i];}
      // ECC / Time
      UShort_t	GetECCTMult()		{return fEXO_ECC_T_Clover.size();}
      UShort_t	GetECCTClover(Int_t i)	{return fEXO_ECC_T_Clover[i];}
      UShort_t	GetECCTCristal(Int_t i)	{return fEXO_ECC_T_Cristal[i];}
      UShort_t	GetECCTTime(Int_t i)	{return fEXO_ECC_T_Time[i];}
      // GOCCE / Energy
      UShort_t	GetGOCCEEMult()			{return fEXO_GOCCE_E_Clover.size();}    // multiplicity of segments hit in one clover
      UShort_t	GetGOCCEEClover(Int_t i)	{return fEXO_GOCCE_E_Clover[i];}
      UShort_t	GetGOCCEECristal(Int_t i)	{return fEXO_GOCCE_E_Cristal[i];}
      UShort_t	GetGOCCEESegment(Int_t i)	{return fEXO_GOCCE_E_Segment[i];}
      UShort_t	GetGOCCEEEnergy(Int_t i)	{return fEXO_GOCCE_E_Energy[i];}
      // GOCCE / Time
      UShort_t	GetGOCCETMult()			{return fEXO_GOCCE_T_Clover.size();}
      UShort_t	GetGOCCETClover(Int_t i)	{return fEXO_GOCCE_T_Clover[i];}
      UShort_t	GetGOCCETCristal(Int_t i)	{return fEXO_GOCCE_T_Cristal[i];}
      UShort_t	GetGOCCETSegment(Int_t i)	{return fEXO_GOCCE_T_Segment[i];}
      UShort_t	GetGOCCETTime(Int_t i)		{return fEXO_GOCCE_T_Time[i];}
      //GeFill
      UShort_t     GetGeFill()               {return fEXO_Fill;}

      ClassDef(TExogamData,1)  // ExogamData structure
	};

#endif
