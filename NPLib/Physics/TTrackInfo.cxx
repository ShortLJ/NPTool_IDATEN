/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 10/06/09                                                 *
 * Last update    : 04/09/09                                                 *
 *---------------------------------------------------------------------------*
 * Decription: This class records all the information concerning the event   *
 *             generators, e.g. vertex of interaction, angles of emitted     *
 *             particles...                                                  *
 *             This class derives from TObject (ROOT) and its aim is to be   *
 *             stored in the output TTree of the G4 simulation               *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + 04/09/09: Add private members for emittance  (N. de Sereville)       *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
using namespace std;

#include "TTrackInfo.h"

ClassImp(TTrackInfo)
////////////////////////////////////////////////////////////////////////////////
TTrackInfo::TTrackInfo(){
}
////////////////////////////////////////////////////////////////////////////////
TTrackInfo::~TTrackInfo(){
}
////////////////////////////////////////////////////////////////////////////////
void TTrackInfo::Clear(){
  
  // emmitted particle
  fTI_Particle_Name.clear();
  fTI_Volume_Name.clear();
  fTI_Kinetic_Energy.clear();
  fTI_Mass.clear();
  fTI_Charge.clear();
  fTI_Z.clear();
  fTI_A.clear();
  fTI_Brho.clear();
  fTI_Momentum.clear();
  fTI_Momentum_X.clear();
  fTI_Momentum_Y.clear();
  fTI_Momentum_Z.clear();
  fTI_PositionX.clear();
  fTI_PositionY.clear();
  fTI_PositionZ.clear();
  fTI_Index.clear();
  fTI_Theta.clear();
  fTI_Phi.clear();
}
////////////////////////////////////////////////////////////////////////////////
void TTrackInfo::Dump() const{
  cout << "--------- Initial Condition Dump ---------" << endl ;
  
  // emmitted particle
  unsigned int size = fTI_Particle_Name.size();
  for(unsigned int i = 0 ; i < size; i ++){
    cout << "\t ---- Particle " << i << " ---- " << endl;
    cout << "\t Particle Name" <<   fTI_Particle_Name[i] << endl;
    cout << "\t Energy" <<   fTI_Kinetic_Energy[i] << endl;
    cout << "\t Momentum Direction: ( "
    << fTI_Momentum_X[i] << " ; "
    << fTI_Momentum_Y[i] << " ; "
    << fTI_Momentum_Z[i] << ")" << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
TVector3 TTrackInfo::GetParticleDirection (const int &i) const {
  return TVector3(  fTI_Momentum_X[i],
                    fTI_Momentum_Y[i],
                    fTI_Momentum_Z[i]).Unit();
}

