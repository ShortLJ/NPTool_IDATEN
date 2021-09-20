/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: ValerianAlcindor  contact address: valcindor@@ikp.tu-darmstadt.de
 *                                                                           *
 * Creation Date  : September 2021                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                  []                                             *
 *  A quite Standard Geant4 EventAction class.                               *
 *  Call the Fill method of the output tree.                                 *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "SteppingAction.hh"
#include "G4UnitsTable.hh"
#include "NPFunction.h"
#include "NPOptionManager.h"
#include "RootOutput.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction() {
  m_record_track = false;
  m_tree         = RootOutput::getInstance()->GetTree();

  m_First       = true;
  ParticleName  = "";
  KineticEnergy = -1000;
  Theta         = -1000;
  Phi           = -1000;
  Mass          = -1000;
  Charge        = -1000;
  Z             = -1000;
  A             = -1000;
  Momentum.setX(-1000);
  Momentum.setY(-1000);
  Momentum.setZ(-1000);
  Position.setX(-1000);
  Position.setY(-1000);
  Position.setZ(-1000);
  LastStepNumber = -1000;
  VolumeName     = "";
  nClear         = 0;
  TrackID        = -1000;

  m_TrackInfo = new TTrackInfo();
  AttachTrackInfo();
  if (!RootOutput::getInstance()->GetTree()->FindBranch("TrackInfo"))
    RootOutput::getInstance()->GetTree()->Branch("TrackInfo", "TTrackInfo",
                                                 &m_TrackInfo);

  // Attach track info class to m_tree
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::AttachTrackInfo() {
  // Reasssigned the branch address
  if (RootOutput::getInstance()->GetTree()->FindBranch("TrackInfo"))
    RootOutput::getInstance()->GetTree()->SetBranchAddress("TrackInfo",
                                                           &m_TrackInfo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step) {
  if (m_record_track) {
    TrackRecording(step);
  }
}

// FIXME Still underdeveloppement
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::TrackRecording(const G4Step* step) {
  if (m_First)
    m_TrackInfo->Clear();
  m_First = false;

  G4Track* track      = step->GetTrack();
  int      StepNumber = track->GetCurrentStepNumber();

  if (eventID
      < G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()) {
    m_TrackInfo->Clear();
    nClear++;
  }

  if (StepNumber == 1) {
    ParticleName  = track->GetParticleDefinition()->GetParticleName();
    KineticEnergy = track->GetDynamicParticle()->GetKineticEnergy();
    Mass          = track->GetDynamicParticle()->GetMass();
    Charge        = track->GetDynamicParticle()->GetCharge();
    Z             = track->GetParticleDefinition()->GetAtomicNumber();
    A             = track->GetParticleDefinition()->GetAtomicMass();
    Momentum      = track->GetDynamicParticle()->GetMomentum();
    Theta         = Momentum.theta() * 180. / M_PI;
    Phi           = Momentum.phi() * 180. / M_PI;
    Position      = track->GetPosition();
    G4VPhysicalVolume* volume = track->GetVolume();
    VolumeName                = volume->GetName();
    TrackID                   = track->GetTrackID();

    double c_light = 299.792458; // To go from T.m to MeV/e
    double Brho = sqrt(KineticEnergy * KineticEnergy + 2 * KineticEnergy * Mass)
                  / (c_light * Charge);

    m_TrackInfo->SetKineticEnergy(KineticEnergy);
    m_TrackInfo->SetTheta(Theta);
    m_TrackInfo->SetPhi(Phi);
    m_TrackInfo->SetMass(Mass);
    m_TrackInfo->SetCharge(Charge);
    m_TrackInfo->SetZ(Z);
    m_TrackInfo->SetA(A);
    m_TrackInfo->SetBrho(Brho);
    TVector3 Mom;
    Mom.SetX(Momentum.x());
    Mom.SetY(Momentum.y());
    Mom.SetZ(Momentum.z());

    m_TrackInfo->SetMomentum(Mom);

    m_TrackInfo->SetMomentumX(Momentum.x());
    m_TrackInfo->SetMomentumY(Momentum.y());
    m_TrackInfo->SetMomentumZ(Momentum.z());

    m_TrackInfo->SetPositionX(Position.x());
    m_TrackInfo->SetPositionY(Position.y());
    m_TrackInfo->SetPositionZ(Position.z());

    m_TrackInfo->SetVolumeName(VolumeName);
    m_TrackInfo->SetIndex(eventID + TrackID);

    if (ParticleName != "e-" && ParticleName != "e+")
      m_TrackInfo->SetParticleName(NPL::ChangeNameFromG4Standard(ParticleName));
    else
      m_TrackInfo->SetParticleName(ParticleName);
  }
  eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
}
