#ifndef NPDCRECONSTRUCTIONMT_H
#define NPDCRECONSTRUCTIONMT_H
#if __cplusplus > 199711L // require c++11 
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta   contact address: matta@lpcaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : October 2020                                            *
 * Last update     : October 2020                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class have all the method needed to analyse Drift Chambers          *
 *  Notation:                                                                *
 *    - The Wire are along the Y axis, there fore providing a measure of the *
 *    X positions                                                            *
 *    - First we find a 2D track in the X-Z plane based on the measurement of* 
 *    several drift radius by minimisation:                                  *
 *      - Let be Di the distance between a trial track and one of the drift  *
 *      circle of coordinate Xi Zi and radius Ri.                            *
 *      - We compute for the track P= sum_i(Di^2/Ri)                         *
 *      - We minimise P to find the correct track                            *
 *      - Warning : this algo assume only one track                          *
 *    - Once the track found, we return the X position at Z=0 and Z=100 to   *
 *    the user.                                                              *
 *    - User can rotate found position depending on the real position of the *
 *    Wires.                                                                 *
 *    - Once the position found for all wire plan, one can look for the      *
 *    intersaction in any given plane. This is done using ResolvePlane.      *
 *    - Resolving plane for two Z plane will provide a point and a direction *
 *****************************************************************************/

// stl
#include <map>
#include <vector>
#include <thread>
#include <mutex>
class TVector3;
class TGraph;
namespace NPL{

  class DCReconstructionMT{
    public:
      DCReconstructionMT(){DCReconstructionMT(1);};
      DCReconstructionMT(unsigned int number_thread);
      ~DCReconstructionMT();

    public:
      // Set number of thread
      // require to stop and reinit the thread for this to be taken into account
      void SetNumberOfThread(unsigned int number_thread){m_nbr_thread=number_thread;};
      // Add a plan to the job list
      void AddPlan(unsigned int ID,const std::vector<double>& X, const std::vector<double>& Z, const std::vector<double>& R);
      // Get results back for the plan with corresponding ID
      double GetResults(unsigned int ID,double& X0,double& X100,double& a, double& b);
      // Build a track in 2D based on drift circle of Radius R and position X,Z
      // return X0(X100) the X position at Z=0 (Z=100)
      // return a and b the coeff of the 2D line
      // when all thread are done
      void BuildTrack2D();

      // Compute X and Y crossing coordinate of 2 plane of Wire
      void ResolvePlane(const TVector3& L,const double& ThetaU ,const TVector3& H, const double& ThetaV, TVector3& PosXY);

      // Function used by the minimizer in BuildTrack2D
      double SumD(const double* parameter );

      // For debugging/optimisation
      // Scan Sumd versus parameter a or b (tovary =0 for a, 1 for b)
      // return a TGraph for display
      TGraph* Scan(double a, double b, int tovary, double minV, double maxV);

    private: // private member used by SumD
      // data to minize index by thread ID
      std::map<unsigned int,unsigned int> sizeX;
      std::map<unsigned int,unsigned int> m_uid; // match thread id and user id
      std::map<unsigned int,const std::vector<double>*> fitX;
      std::map<unsigned int,const std::vector<double>*> fitZ;
      std::map<unsigned int,const std::vector<double>*> fitR;
      // Computed value indexed by user ID
      std::map<unsigned int,double> m_minimum;
      std::map<unsigned int,double> m_X0;
      std::map<unsigned int,double> m_X100;
      std::map<unsigned int,double> m_a;
      std::map<unsigned int,double> m_b;

      // used by resolve plane
      long double av,bv,au,bu;
      double xM,yM;

    private: // Thread Pool defined if C++11 is available
      unsigned int m_nbr_thread;
      std::vector<std::thread> m_ThreadPool;
      std::vector<bool> m_Ready;
      bool m_stop;
      std::mutex m_mtx;

    public: // Init the Thread Pool

      void InitThreadPool(); 
      void StartThread(unsigned int);
      void StopThread();
      bool IsDone();

  };
}
#endif//c++11
#endif//ndef
