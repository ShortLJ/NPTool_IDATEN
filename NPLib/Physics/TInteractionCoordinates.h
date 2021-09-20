#ifndef __INTERACTIONCOORDINATES__
#define __INTERACTIONCOORDINATES__
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
 * Last update    : 01/09/2021 Valerian Alcindor adding particle name to     * 
 *                  interaction coordinate for easier g4 analysis simulation *
 *---------------------------------------------------------------------------*
 * Decription: This class mainly records the coordinates of interaction      *
 *             between a particle and a detector.                            *
 *             This class derives from TObject (ROOT) and its aim is to be   *
 *             stored in the output TTree of the G4 simulation               *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/



#include "NPFunction.h"
#include <vector>
#include "TObject.h"
#include <iostream>

using namespace std ;


class TInteractionCoordinates : public TObject{
  private:
    // TrackID or index for correlations
    vector<int> fDetected_Index;
    // Detected particle properties (before interactions in the target)
    // Energy and Time
    vector<double> fDetected_Energy;
    vector<double> fDetected_Time;
    // Vertex of interaction
    vector<double>  fDetected_Position_X;
    vector<double>  fDetected_Position_Y;
    vector<double>  fDetected_Position_Z;
    // Particle angles
    vector<double>  fDetected_Angle_Theta;
    vector<double>  fDetected_Angle_Phi;
    // Particle characteristics
    vector<std::string>  fDetected_Particle_Name;
    vector<int>          fDetected_A;
    vector<int>          fDetected_Z;
    vector<double>       fDetected_Mass;
    vector<int>          fDetected_Charge;
    vector<double>       fDetected_Brho;

  public:
    TInteractionCoordinates();
    virtual ~TInteractionCoordinates();

    void  Clear();
    void  Clear(const Option_t*) {};
    void  Dump() const;


    void SetInteraction(const double& Energy, const double&Time, const double& PositionX, const double& PositionY, const double& PositionZ,const double& Theta, const double& Phi){
      fDetected_Energy.push_back(Energy);
      fDetected_Time.push_back(Time);
      fDetected_Position_X.push_back(PositionX);
      fDetected_Position_Y.push_back(PositionY);
      fDetected_Position_Z.push_back(PositionZ);
      fDetected_Angle_Theta.push_back(Theta);
      fDetected_Angle_Phi.push_back(Phi);
    }

    void SetInteraction(const int& Index, const double& Energy, const double&Time, const double& PositionX, const double& PositionY, const double& PositionZ,const double& Theta, const double& Phi){
      fDetected_Index.push_back(Index);
      fDetected_Energy.push_back(Energy);
      fDetected_Time.push_back(Time);
      fDetected_Position_X.push_back(PositionX);
      fDetected_Position_Y.push_back(PositionY);
      fDetected_Position_Z.push_back(PositionZ);
      fDetected_Angle_Theta.push_back(Theta);
      fDetected_Angle_Phi.push_back(Phi);
    }

    void SetInteraction(const int& Index, const double& Energy, const double&Time, const double& PositionX, const double& PositionY, const double& PositionZ,const double& Theta, const double& Phi, const std::string &ParticleName, const int &A, const int &Z, const double &Mass, const int &Charge, const double &Brho){
      fDetected_Index.push_back(Index);
      fDetected_Energy.push_back(Energy);
      fDetected_Time.push_back(Time);
      fDetected_Position_X.push_back(PositionX);
      fDetected_Position_Y.push_back(PositionY);
      fDetected_Position_Z.push_back(PositionZ);
      fDetected_Angle_Theta.push_back(Theta);
      fDetected_Angle_Phi.push_back(Phi);
      if(ParticleName != "e-" && ParticleName != "e+")
        fDetected_Particle_Name.push_back(NPL::ChangeNameFromG4Standard(ParticleName));
      else
        fDetected_Particle_Name.push_back(ParticleName);
      fDetected_A.push_back(A);
      fDetected_Z.push_back(Z);
      fDetected_Mass.push_back(Mass);
      fDetected_Charge.push_back(Charge);
      fDetected_Brho.push_back(Brho);
    }

    /////////////////////           SETTERS           ////////////////////////
    // Incident particle properties (before interactions in the target)
    // Vertex of interaction
    void SetDetectedPositionX(const double& PositionX)      {fDetected_Position_X.push_back(PositionX);}//!
    void SetDetectedPositionY(const double& PositionY)      {fDetected_Position_Y.push_back(PositionY);}//!
    void SetDetectedPositionZ(const double& PositionZ)      {fDetected_Position_Z.push_back(PositionZ);}//!
    // Incident particle angles
    void SetDetectedAngleTheta(const double& AngleTheta)  {fDetected_Angle_Theta.push_back(AngleTheta);}//!
    void SetDetectedAnglePhi(const double& AnglePhi)      {fDetected_Angle_Phi.push_back(AnglePhi);}//!

    void SetDetectedParticleName(const std::string& ParticleName)      {fDetected_Particle_Name.push_back(ParticleName);}//!
    void SetDetectedA(const int& A)      {fDetected_A.push_back(A);}//!
    void SetDetectedZ(const int& Z)      {fDetected_Z.push_back(Z);}//!
    void SetDetectedMass(const double& Mass)      {fDetected_Mass.push_back(Mass);}//!
    void SetDetectedCharge(const int& Charge)      {fDetected_Charge.push_back(Charge);}//!
    void SetDetectedBrho(const double& Brho)      {fDetected_Brho.push_back(Brho);}//!

    /////////////////////           GETTERS           ////////////////////////
    // Number of interactions (multiplicity)
    int    GetDetectedMultiplicity() const     {return fDetected_Position_X.size();}
    // Incident particle properties (before interactions in the target)
    // Energy and Time
    double GetEnergy(const int& i) const {return fDetected_Energy[i];}//!
    double GetTime(const int& i) const   {return fDetected_Time[i];}//!
    // Vertex of interaction
    double GetDetectedPositionX(const int& i) const   {return fDetected_Position_X[i];}//!
    double GetDetectedPositionY(const int& i) const   {return fDetected_Position_Y[i];}//!
    double GetDetectedPositionZ(const int& i) const   {return fDetected_Position_Z[i];}//!
    // Incident particle angles
    double GetDetectedAngleTheta(const int& i) const {return fDetected_Angle_Theta[i];}//!
    double GetDetectedAnglePhi(const int& i) const   {return fDetected_Angle_Phi[i];}//!

    std::string GetParticleName(const int& i) const   {return fDetected_Particle_Name[i];}//!

    int GetA(const int& i) const   {return fDetected_A[i];}//!
    int GetZ(const int& i) const   {return fDetected_Z[i];}//!
    double GetMass(const int& i) const   {return fDetected_Mass[i];}//!
    int GetCharge(const int& i) const   {return fDetected_Charge[i];}//!
    double GetBrho(const int& i) const   {return fDetected_Brho[i];}//!

    ClassDef(TInteractionCoordinates, 2) // InteractionCoordinates structure
};

#endif
