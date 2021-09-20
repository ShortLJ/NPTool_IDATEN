/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace address: pierre.morfouace2@cea.fr       *
 *                                                                           *
 * Creation Date  : Octobre 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Use to kill the beam track and replace it with the fission products       *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <string>
#include <set>
#include "NPFunction.h"
#include "FissionDecay.hh"
#include "NPOptionManager.h"
#include "NPInputParser.h"
#include "G4VPhysicalVolume.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "RootOutput.h"

using namespace NPS;
////////////////////////////////////////////////////////////////////////////////
FissionDecay::FissionDecay(G4String modelName,G4Region* envelope) :
  G4VFastSimulationModel(modelName, envelope) {
    m_PreviousEnergy=0 ;
    m_PreviousLength=0 ;
    m_FissionConditions = new TFissionConditions();
    
    ReadConfiguration();
    //AttachFissionConditions();
    //if(!RootOutput::getInstance()->GetTree()->FindBranch("FissionConditions"))
    //  RootOutput::getInstance()->GetTree()->Branch("FissionConditions", "TFissionConditions", &m_FissionConditions);
  }


////////////////////////////////////////////////////////////////////////////////
FissionDecay::FissionDecay(G4String modelName) :
  G4VFastSimulationModel(modelName) {
  }

////////////////////////////////////////////////////////////////////////////////
FissionDecay::~FissionDecay() {
}

////////////////////////////////////////////////////////////////////////////////
void FissionDecay::AttachFissionConditions(){
  if(RootOutput::getInstance()->GetTree()->FindBranch("FissionConditions"))
    RootOutput::getInstance()->GetTree()->SetBranchAddress("FissionConditions", &m_FissionConditions);
}

////////////////////////////////////////////////////////////////////////////////
void FissionDecay::ReadConfiguration(){
  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());
  m_FissionDecay.ReadConfiguration(input);

  if(m_FissionDecay.GetFissionToken()>0){
    std::string Mother = m_FissionDecay.GetCompoundName(); 
    m_CompoundParticle = NPL::Particle(Mother);
    m_CompoundName = NPL::ChangeNameToG4Standard(Mother,true);
    AttachFissionConditions();
    if(!RootOutput::getInstance()->GetTree()->FindBranch("FissionConditions")){
      RootOutput::getInstance()->GetTree()->Branch("FissionConditions", "TFissionConditions", &m_FissionConditions);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
G4bool FissionDecay::IsApplicable( const G4ParticleDefinition& particleType) {
  
  m_CurrentName = particleType.GetParticleName();
  // Extract Ex from name
  if(m_CurrentName.find("[")!=std::string::npos)
    m_ExcitationEnergy = atof(m_CurrentName.substr(m_CurrentName.find("[")+1,m_CurrentName.find("]")-1).c_str())*keV;
  else
    m_ExcitationEnergy=0;

  // Strip name from excitation energy
  m_CurrentName = m_CurrentName.substr(0,m_CurrentName.find("["));
  if (m_CompoundName==m_CurrentName) {
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
G4bool FissionDecay::ModelTrigger(const G4FastTrack& fastTrack) {
  bool Trigger = true;

  m_FissionConditions->Clear();

  m_PreviousEnergy=fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  // To be coded;

  return Trigger;
}

////////////////////////////////////////////////////////////////////////////////
void FissionDecay::DoIt(const G4FastTrack& fastTrack,G4FastStep& fastStep){
  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  G4ThreeVector pdirection = PrimaryTrack->GetMomentum().unit();
  G4ThreeVector localdir = fastTrack.GetPrimaryTrackLocalDirection();

  G4ThreeVector worldPosition = PrimaryTrack->GetPosition();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();

  double energy = PrimaryTrack->GetKineticEnergy();
  double time = PrimaryTrack->GetGlobalTime();

  // Randomize within the step
  // Assume energy loss is linear within the step
  // Assume no scattering
  double rand =  G4RandFlat::shoot(); 
  double length = rand*(m_PreviousLength); 
  energy += (1-rand)*(m_PreviousEnergy-energy); 
  G4ThreeVector ldir = pdirection;
  ldir*=length;
  localPosition = localPosition - ldir; 
  //////////////////////////////////////////////////
  //////Define the kind of particle to shoot////////
  //////////////////////////////////////////////////
  std::vector<NPL::Particle> FissionFragment;
  std::vector<double> Ex;
  std::vector<double> DEK;
  std::vector<double> DPx;
  std::vector<double> DPy;
  std::vector<double> DPz;
  double TKE, KE1, KE2;


  m_FissionDecay.GenerateEvent(NPL::ChangeNameFromG4Standard(m_CurrentName),m_ExcitationEnergy,energy,
      pdirection.x(),pdirection.y(),pdirection.z(),
      FissionFragment, Ex,DEK,DPx,DPy,DPz,
      TKE, KE1, KE2);

  /////////////////////////////////////////////////
  // Fillion the attached Fission condition Tree //
  /////////////////////////////////////////////////
  // Fissionning system
  m_FissionConditions->SetZ_CN(m_CompoundParticle.GetZ());
  m_FissionConditions->SetA_CN(m_CompoundParticle.GetA());
  m_FissionConditions->SetEx_CN(m_ExcitationEnergy);
  m_FissionConditions->SetELab_CN(energy);
  m_FissionConditions->SetThetaLab_CN(pdirection.theta()*180./3.1415);

  // Fission Process
  m_FissionConditions->Set_TKE(TKE);
  m_FissionConditions->Set_KE1(KE1);
  m_FissionConditions->Set_KE2(KE2);

  G4ParticleDefinition* FissionFragmentDef; 
  unsigned int size = FissionFragment.size();

  if(size == 0)
    return;

  // Get Neutron Multiplicity
  int Zsum = FissionFragment[0].GetZ() + FissionFragment[1].GetZ();
  int Asum = FissionFragment[0].GetA() + FissionFragment[1].GetA();
  if(Zsum == m_CompoundParticle.GetZ())
    m_FissionConditions->SetNeutronMultiplicity(m_CompoundParticle.GetA()-Asum);

  for(unsigned int i = 0 ; i < size ; i++){
    // Get the decaying particle
    int FFZ = FissionFragment[i].GetZ();
    int FFA = FissionFragment[i].GetA();
    FissionFragmentDef=NULL;

    // Set the momentum direction
    G4ThreeVector Momentum (DPx[i],DPy[i],DPz[i]);
    Momentum=Momentum.unit();

    double Brho = FissionFragment[i].GetBrho();
    double KineticEnergy = FissionFragment[i].GetEnergy();

    m_FissionConditions->SetFragmentZ(FFZ);
    m_FissionConditions->SetFragmentA(FFA);
    //m_FissionConditions->SetFragmentKineticEnergy(DEK[i]);
    m_FissionConditions->SetFragmentKineticEnergy(KineticEnergy);
    m_FissionConditions->SetFragmentBrho(Brho);
    m_FissionConditions->SetFragmentTheta(Momentum.theta()/deg);
    m_FissionConditions->SetFragmentPhi(Momentum.phi()/deg);
    m_FissionConditions->SetFragmentMomentumX(DPx[i]);
    m_FissionConditions->SetFragmentMomentumX(DPy[i]);
    m_FissionConditions->SetFragmentMomentumX(DPz[i]);

    // neutral particle
    if(FFZ==0){
      if(FFA==1)
        FissionFragmentDef=G4ParticleTable::GetParticleTable()->FindParticle("neutron");

      else if(FFA==0){
        FissionFragmentDef=G4ParticleTable::GetParticleTable()->FindParticle("gamma");
      }

    }
    // proton
    else if (FFZ==1 && FFA==1 )
      FissionFragmentDef=G4ParticleTable::GetParticleTable()->FindParticle("proton");
    // the rest
    else
      FissionFragmentDef=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(FFZ, FFA, Ex[i]);


    G4DynamicParticle DynamicFissionFragment(FissionFragmentDef,Momentum,DEK[i]);
    fastStep.CreateSecondaryTrack(DynamicFissionFragment, localPosition, time);
  }

  if(size){
    // Set the end of the step conditions
    fastStep.SetPrimaryTrackFinalKineticEnergyAndDirection(0,pdirection);
    fastStep.SetPrimaryTrackFinalPosition(worldPosition);  
    fastStep.SetTotalEnergyDeposited(0);
    fastStep.SetPrimaryTrackFinalTime (time);
    fastStep.KillPrimaryTrack();
    fastStep.SetPrimaryTrackPathLength(0.0);
  }

}
