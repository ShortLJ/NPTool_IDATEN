#ifndef __TRACKINFO__
#define __TRACKINFO__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
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

// STL Header
#include <cmath>
#include <string>
#include <vector>
using namespace std;

// Root Header
#include "TObject.h"
#include "TVector3.h"

// NPTOOL headers
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

class TTrackInfo : public TObject {
private:
  // Track info

public:
  // Particles info
  vector<string> fTI_Particle_Name;
  vector<string> fTI_Volume_Name;
  vector<double> fTI_Kinetic_Energy;
  vector<double> fTI_Theta;
  vector<double> fTI_Phi;

  vector<double> fTI_Mass;
  vector<double> fTI_Charge;
  vector<double> fTI_Z;
  vector<double> fTI_A;
  vector<double> fTI_Brho;

  // vector<double> fTI_Spin;

  vector<TVector3> fTI_Momentum;
  vector<double>   fTI_Momentum_X;
  vector<double>   fTI_Momentum_Y;
  vector<double>   fTI_Momentum_Z;

  vector<double> fTI_PositionX;
  vector<double> fTI_PositionY;
  vector<double> fTI_PositionZ;

  vector<int> fTI_Index;

  // vector<double> fTI_;

public:
  TTrackInfo();
  virtual ~TTrackInfo();

  void Clear();
  void Clear(const Option_t*) { Clear(); };
  void Dump() const;

  // emmitted particle
  void SetParticleName(const string& Particle_Name) {
    fTI_Particle_Name.push_back(Particle_Name);
  }
  void SetVolumeName(const string& Volume_Name) {
    fTI_Volume_Name.push_back(Volume_Name);
  }
  void SetKineticEnergy(const double& Kinetic_Energy) {
    fTI_Kinetic_Energy.push_back(Kinetic_Energy);
  }
  void SetTheta(const double& Theta) {
    fTI_Theta.push_back(Theta);
  }
  void SetPhi(const double& Phi) {
    fTI_Phi.push_back(Phi);
  }
  void SetMass(const double& Mass) { fTI_Mass.push_back(Mass); }
  void SetCharge(const double& Charge) { fTI_Charge.push_back(Charge); }
  void SetA(const double& A) { fTI_A.push_back(A); }
  void SetZ(const double& Z) { fTI_Z.push_back(Z); }
  void SetBrho(const double& Brho) { fTI_Brho.push_back(Brho); }
  void SetMomentum(const TVector3& Momentum) {
    fTI_Momentum.push_back(Momentum);
  }
  void SetMomentumX(const double& Momentum_X) {
    fTI_Momentum_X.push_back(Momentum_X);
  }
  void SetMomentumY(const double& Momentum_Y) {
    fTI_Momentum_Y.push_back(Momentum_Y);
  }
  void SetMomentumZ(const double& Momentum_Z) {
    fTI_Momentum_Z.push_back(Momentum_Z);
  }
  void SetPositionX(const double& X) { fTI_PositionX.push_back(X); }
  void SetPositionY(const double& Y) { fTI_PositionY.push_back(Y); }
  void SetPositionZ(const double& Z) { fTI_PositionZ.push_back(Z); }
  void SetIndex(const int& Index) { fTI_Index.push_back(Index); }

  // emmitted particle
  int    GetParticleMultiplicity() const { return fTI_Kinetic_Energy.size(); }
  string GetParticleName(const int& i) const { return fTI_Particle_Name[i]; }
  double GetKineticEnergy(const int& i) const { return fTI_Kinetic_Energy[i]; }
  TVector3 GetMomentum(const int& i) const { return fTI_Momentum[i]; }
  double   GetMomentumX(const int& i) const { return fTI_Momentum_X[i]; }
  double   GetMomentumY(const int& i) const { return fTI_Momentum_Y[i]; }
  double   GetMomentumZ(const int& i) const { return fTI_Momentum_Z[i]; }

  TVector3 GetParticleDirection(const int& i) const;

  double GetThetaLab_WorldFrame(const int& i) const {
    return (GetParticleDirection(i).Angle(TVector3(0, 0, 1))) / deg;
  }

  unsigned int GetEmittedMult() const { return fTI_Particle_Name.size(); }

  ClassDef(TTrackInfo, 1) // TrackInfo structure
};

#endif
