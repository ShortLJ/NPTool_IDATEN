#ifndef NPDCRECONSTRUCTION_H
#define NPDCRECONSTRUCTION_H
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
#include <vector>
#include "TVector3.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
class TGraph;
namespace NPL{

  class DCReconstruction{
    public:
      DCReconstruction();
      ~DCReconstruction();

    public:
      // Build a track in 2D based on drift circle of Radius R and position X,Z
      // return X0(X100) the X position at Z=0 (Z=100)
      // return a and b the coeff of the 2D line
      double BuildTrack2D(const std::vector<double>& X,const std::vector<double>& Z,const std::vector<double>& R,double& X0,double& X100,double& a, double& b );

      // Compute X and Y crossing coordinate of 2 plane of Wire
      void ResolvePlane(const TVector3& L,const double& ThetaU ,const TVector3& H, const double& ThetaV, TVector3& PosXY);

      // Function used by the minimizer in BuildTrack2D
      double SumD(const double* parameter );

      
      // For debugging/optimisation
      // Scan Sumd versus parameter a or b (tovary =0 for a, 1 for b)
      // return a TGraph for display
      TGraph* Scan(double a, double b, int tovary, double minV, double maxV);


    private: // private member used by SumD
      ROOT::Math::Minimizer* m_min;
      ROOT::Math::Functor    m_func;
      const std::vector<double>* fitX;
      const std::vector<double>* fitZ;
      const std::vector<double>* fitR;

      // used by SumD
      unsigned int sizeX ;
      double P,p,a,b,ab,a2,c,d,x,z,r;
      // used by BuildTrack
      double ai,bi;
      double parameter[2];
      // used by resolve plane
      long double av,bv,au,bu;
      double xM,yM;

  };
}

#endif
