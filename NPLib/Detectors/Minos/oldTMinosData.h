#ifndef __MinosDATA__
#define __MinosDATA__
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: tronchin@lpccaen.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Raw data                                    *
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
#include "TGraph.h"

using namespace std;

class TMinosData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
  
    // Pads
    vector<UShort_t>   fMinos_PadNumber;
    vector<Double_t>   fMinos_PadX;
    vector<Double_t>   fMinos_PadY;
    /* vector<Double_t>   fMinos_PadCharge; */
    /* vector<Double_t>   fMinos_Charge; */
    /* vector<Double_t>   fMinos_Time; */
    vector< vector<Double_t> >   fMinos_Charge;
    vector< vector<Double_t> >   fMinos_Time;
    /* vector<Double_t>   fMinos_DriftTime; */
    /* vector<UShort_t>   fMinos_EventNumber; */
/* //From Santamaria:*/


/* //from simulation*/
/*   vector<double> x_tpc,y_tpc,z_tpc,e_tpc;*/
/*   vector<double> x_trigger,y_trigger,z_trigger,e_trigger;*/
/* vector<double> x_tar,y_tar,z_tar,e_tar;*/
/* vector<double> x_ch,y_ch,z_ch,e_ch;*/
/* vector<double> x_win,y_win,z_win,e_win;*/
/* vector<double> x_InRoh,y_InRoh,z_InRoh,e_InRoh;*/
/* vector<double> x_OutRoh,y_OutRoh,z_OutRoh,e_OutRoh;*/
/* vector<double> x_Kap,y_Kap,z_Kap,e_Kap;*/
/* double Et_tpc_tot;*/
/* vector<double> Et_tar,Et_ch,Et_tpc,Et_trigger,Et_win,Et_InnerRohacell, Et_OuterRohacell, Et_Kapton;*/
/* vector<int> A, Z;*/
/* vector<int> trackID, parentID;*/

/* */     //unuseful, cause nptool should make already that*/
/* //initial conditions*/
/* //double x0,y0,z0,theta0,phi0,energy0;*/
/* vector<double> x0, y0, z0, theta0, phi0, energy0;*/
/* vector<bool> detection;*/
/* int event;*/
/* */
 vector<Double_t> MINOSx_0; //!
 vector<Double_t> MINOSy_0; //!
 vector<Double_t> MINOSz_0; //! 
 vector<Double_t> MINOS_D_min; //! 
 vector<Double_t> MINOS_Radius; //! 
 vector<Double_t> MINOS_NumberTracks; //! 
 vector<Double_t> theta0; //! 
 vector<Double_t> phi0; //! 
 vector<Double_t> energy0; //! 
 
 //For take fitpar values

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
    /* void Dump() const; */


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////

     // Minos Pads
    inline void SetCharge(const UShort_t& Pad,/*const vector<Double_t>& Charge,const vector<Double_t>& Time,*/  const Double_t& X,const Double_t& Y/*,const Double_t& PadCharge*/){
      fMinos_PadNumber.push_back(Pad);
      fMinos_PadX.push_back(X);
      fMinos_PadY.push_back(Y);
    };//!

    inline void AddChargePoint(const vector< Double_t >& Q,const vector<Double_t>& T){
      fMinos_Charge.push_back(Q);
      fMinos_Time.push_back(T);
    };//!

    /* inline void AddChargePoint(const double& Q,const double& T){ */
    /*   fMinos_Charge.push_back(Q); */
    /*   fMinos_Time.push_back(); */
    /* };//! */
    
   //

    //Setters for position vertex and obsv in experiment analysis

    // Position
    inline void SetVertexPos(const Double_t& x,const Double_t& y,const Double_t& z)	{
      MINOSx_0.push_back(x);
      MINOSy_0.push_back(y);
      MINOSz_0.push_back(z);     
    };//!

    // Min Distance
    inline void SetD_min(const Double_t& dmin)     {
      MINOS_D_min.push_back(dmin);
    };//!

    //////////////////////    GETTERS    ////////////////////////

      inline int GetPadMult()
        {return fMinos_PadNumber.size() ;}//!

      /* inline int GetEventNumberMult() */
      /*   {return fMinos_EventNumber.size() ;}//! */

      inline double GetPadX(const unsigned int&i) const
        {return fMinos_PadX[i] ;}//! 
      inline double GetPadY(const unsigned int&i) const
        {return fMinos_PadY[i] ;}//!
      inline UShort_t GetPadNbr(const unsigned int&i) const
        {return fMinos_PadNumber[i] ;}//!

       /* inline double GetPadCharge(const unsigned int&i) const */
       /*  {return fMinos_PadCharge[i] ;}//! */
       
       /* inline double GetPadTime(const unsigned int&i) const */
        /* {return fMinos_DriftTime[i] ;}//! */
       
      inline vector<double>& GetCharge(const unsigned int&i)
        {return fMinos_Charge[i] ;}//!
      inline vector<double>& GetTime(const unsigned int&i)
        {return fMinos_Time[i] ;}//!

      /* TGraph* GetChargeAsGraph(const unsigned int&i) ; //! */      

      /* inline vector<double> GetCharge() */
      /*   {return fMinos_PadChar   TGraph* GetEnergyAsGraph();ge() ;}//! */



    // Position
    inline Double_t GetVertexPos() const
    {return MINOSz_0[0] ;}//!
    inline Double_t GetVertexPosX() const
    {return MINOSx_0[0] ;}//!
    inline Double_t GetVertexPosY() const
    {return MINOSy_0[0] ;}//!
    inline Double_t GetVertexPosZ() const
    {return MINOSz_0[0] ;}//!

    // Min Distance
    inline Double_t GetD_min() const
    {return MINOS_D_min[0] ; }//!

    // Charge


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TMinosData,1)  // MinosData structure
};

#endif
