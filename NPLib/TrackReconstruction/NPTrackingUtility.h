#ifndef NPTrackingUtility_H
#define NPTrackingUtility_H
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
#include "TVector3.h"
namespace NPL{
  //////////////////////////////////////////////////////////////////////////////
  // return the minimum distance between v and w defined respectively by points 
  // v1,v2 and w1 w1
  // Also compute the best crossing position BestPosition, i.e. average position
  // at the minimum distance.
  double MinimumDistanceTwoLines(const TVector3& v1,const TVector3& v2, const TVector3& w1, const TVector3& w2, TVector3& BestPosition, TVector3& delta);

  //////////////////////////////////////////////////////////////////////////////
  // return the minimum distance between the line defines by v1,v2 and the point
  // in space x
  // demo is here: https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  double MinimumDistancePointLine(const TVector3& v1, const TVector3& v2, const TVector3& x);
}


#endif
