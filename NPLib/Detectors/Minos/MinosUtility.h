#ifndef MINOSUTILITY_H
#define MINOSUTILITY_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta   contact address: matta@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : october 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Minos class utility. Hold various usefull fonction to analyse MINOS      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include<vector>
#include "TMinosData.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"

namespace NPL{

  class MinosUtility{
    public:
      MinosUtility();
      ~MinosUtility(){};

    public:
      // take the vector describing the charge, perform a fit and return the charge and time of the trace
      double Calibrate(const std::vector<unsigned short>* T,const std::vector<unsigned short>* Q, const unsigned int& i , double& time, double& charge);
      // This function describe the shape of signal 
      double signal_shape(const double& x, const double* p);
      double signal_chi2(const double* p);
      void   sample_signal(); // take one point every sampling
      void   find_maxQ();

    private:

      std::vector<std::pair<unsigned short,unsigned short>> m_vTQ;
      ROOT::Math::Minimizer* m_signal_min; // minimiser used for the signal fitting
      ROOT::Math::Functor    m_signal_func; // functor passed to the miniser
      const std::vector<unsigned short>* m_fitSignalT;
      const std::vector<unsigned short>* m_fitSignalQ;
      unsigned int m_signal_size,m_minbin,m_maxbin,m_qmax,m_qmaxbin,m_guess_t0_bin;
      unsigned int m_counter;
      double m_TimeBin; // time between two sample in the waveform (10 to 30 ns depending on the FEM clock)
      double m_ShapingTime; // Setting of the preamp on the AGET
      unsigned int m_Baseline; // Waveform offset ~250
      double m_Tau;// Time constant of the pad signal expressed in bin (depending on TimeBin and ShapingTime)
      unsigned int m_Sampling;// Allow to take less point in the wave form to perform the fit to improve speed. 10 gives good results and is twice faster than 1 
    
    public:
      void SetParameters(double TimeBin,double ShapingTime, unsigned int Baseline, unsigned int Sampling){
        m_TimeBin=TimeBin;
        m_ShapingTime=ShapingTime;
        m_Tau=m_ShapingTime/(log(2)*m_TimeBin);
        m_Baseline = Baseline;
        m_Sampling = Sampling;};
  };



}  

#endif
