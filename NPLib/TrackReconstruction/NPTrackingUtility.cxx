/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta  contact address: matta@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : July 2020                                               *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with finding all the track event by event                *
 *****************************************************************************/

#include "NPTrackingUtility.h"

////////////////////////////////////////////////////////////////////////////////
double NPL::MinimumDistanceTwoLines(const TVector3& v1,const TVector3& v2, const TVector3& w1, const TVector3& w2, TVector3& BestPosition, TVector3& delta){
  TVector3 v = v2-v1;
  TVector3 w = w2-w1;
  // Finding best position
  TVector3 e = v1-w1;
  double A = -(v.Mag2()*w.Mag2()-(v.Dot(w)*v.Dot(w)));
  double s = (-v.Mag2()*(w.Dot(e))+(v.Dot(e))*(w.Dot(v)))/A;
  double t = (w.Mag2()*(v.Dot(e))-(w.Dot(e)*w.Dot(v)))/A;
  double d = sqrt((e+v*t-w*s).Mag2());
 
  BestPosition = 0.5*(v1+t*v+w1+s*w);
  delta = (v1+t*v-w1-s*w);
  return d;
  }
////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // return the minimum distance between the line defines by v1,v2 and the point
  // in space x
  // demo is here: https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
double NPL::MinimumDistancePointLine(const TVector3& v1, const TVector3& v2, const TVector3& x){
    TVector3 w1 = x-v1;
    TVector3 w2 = x-v2;
    TVector3 w = w1.Cross(w2);
    return w.Mag()/(v2-v1).Mag();
    }

