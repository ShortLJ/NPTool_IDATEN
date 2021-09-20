/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : Octobre 2017                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Use to kill the beam track and replace it with the reaction product       *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#ifndef BeamReaction_h
#define BeamReaction_h

#include "G4VFastSimulationModel.hh"
#include "G4Abla.hh"
#include "G4AblaInterface.hh"
#include "G4Fragment.hh"
#include "PhysicsList.hh"
#include "NPReaction.h"
#include "NPQFS.h"
#include "TReactionConditions.h"
class G4VPhysicalVolume;
namespace NPS{
  enum ReactionType{
    TwoBody,
    QFS,
    Fusion
    };

  class BeamReaction : public G4VFastSimulationModel{
    public:
      BeamReaction (G4String, G4Region*);
      BeamReaction (G4String);
      ~BeamReaction ();

    public:
      void ReadConfiguration();
      G4bool IsApplicable(const G4ParticleDefinition&);
      G4bool ModelTrigger(const G4FastTrack &);
      void DoIt(const G4FastTrack&, G4FastStep&);
 
    private:
      NPL::Reaction m_Reaction;
      NPL::QFS m_QFS;
      string m_BeamName;
      int m_ReactionType;
      G4AblaInterface* ABLA;

      bool   m_active;// is the process active
      bool   m_shoot;
      double m_StepSize;
      double m_Z;
      double m_S;
      double m_rand;
      double m_length;
      int    m_Parent_ID;
      double SlowDownBeam(const G4ParticleDefinition* Beam, double IncidentEnergy, double Thickness,G4Material* Material);
    
    private:// specific for the simple case of fusion
      string m_TargetNuclei;
      string m_FusionProduct;
      double m_FusionExcitation;
  
   private:
     TReactionConditions* m_ReactionConditions;
 
   public:
    void AttachReactionConditions();
    void SetStepSize(double step){m_StepSize=step;};
  };
}


#endif 
