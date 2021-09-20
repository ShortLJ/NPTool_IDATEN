/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  pierre.morfouace2@cea.fr               *
 *                                                                           *
 * Creation Date  : Octobre 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
*                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#ifndef FissionDecay_h
#define FissionDecay_h

#include "G4VFastSimulationModel.hh"
#include "PhysicsList.hh"
#include "NPFissionDecay.h"
#include "TFissionConditions.h"
class G4VPhysicalVolume;

////////////////////////////////////////////////////////////////////////////////
namespace NPS{

class FissionDecay : public G4VFastSimulationModel{
  public:
    FissionDecay (G4String, G4Region*);
    FissionDecay (G4String);
    ~FissionDecay();

  public:
    void ReadConfiguration();
    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    virtual G4bool ModelTrigger(const G4FastTrack &);
    virtual void DoIt(const G4FastTrack&, G4FastStep&);

  private:
    NPL::FissionDecay m_FissionDecay;
    NPL::Particle m_CompoundParticle;
    std::string m_CompoundName;
    std::string m_CurrentName;
    double m_ExcitationEnergy;
    double m_PreviousEnergy;
    double m_PreviousLength;

  private:
    TFissionConditions* m_FissionConditions;

  public:
    void AttachFissionConditions();
};
}

#endif 
