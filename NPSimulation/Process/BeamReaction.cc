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

#include "BeamReaction.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4IonTable.hh"
#include "G4UserLimits.hh"
#include "NPFunction.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "RootOutput.h"
#include "TLorentzVector.h"
#include "NPSFunction.hh"
#include <Randomize.hh>
#include <iostream>
#include <string>
#include "G4UserLimits.hh"

////////////////////////////////////////////////////////////////////////////////
NPS::BeamReaction::BeamReaction(G4String modelName, G4Region* envelope)
  : G4VFastSimulationModel(modelName, envelope) {
    ReadConfiguration();
    m_shoot=false;
    m_rand=0;
    m_Z=0;

    ABLA = new G4AblaInterface();
  }

////////////////////////////////////////////////////////////////////////////////
NPS::BeamReaction::BeamReaction(G4String modelName)
  : G4VFastSimulationModel(modelName) {}

  ////////////////////////////////////////////////////////////////////////////////
  NPS::BeamReaction::~BeamReaction() {}

  ////////////////////////////////////////////////////////////////////////////////
  void NPS::BeamReaction::AttachReactionConditions() {
    // Reasssigned the branch address
    if (RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
      RootOutput::getInstance()->GetTree()->SetBranchAddress(
          "ReactionConditions", &m_ReactionConditions);
  }

////////////////////////////////////////////////////////////////////////////////
void NPS::BeamReaction::ReadConfiguration() {
  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());

  if(input.GetAllBlocksWithToken("TwoBodyReaction").size()>0) m_ReactionType =TwoBody;
  else if(input.GetAllBlocksWithToken("QFSReaction").size()>0) m_ReactionType =QFS;
  else if(input.GetAllBlocksWithToken("FusionReaction").size()>0) m_ReactionType =Fusion;

  // Two body
  if (m_ReactionType==TwoBody ) {
    m_Reaction.ReadConfigurationFile(input);
    m_BeamName = NPL::ChangeNameToG4Standard(m_Reaction.GetParticle1()->GetName());
    if(m_Reaction.GetParticle3()->GetName() != ""){
      m_active = true;
      m_ReactionConditions = new TReactionConditions();
      AttachReactionConditions();
      if (!RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
        RootOutput::getInstance()->GetTree()->Branch(
            "ReactionConditions", "TReactionConditions", &m_ReactionConditions);
    }
  } 

  // QFS
  else  if (m_ReactionType==QFS) {
    m_QFS.ReadConfigurationFile(input);
    m_BeamName = NPL::ChangeNameToG4Standard(m_QFS.GetParticleA()->GetName());
    m_active = true;
    m_ReactionConditions = new TReactionConditions();
    AttachReactionConditions();
    if (!RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
      RootOutput::getInstance()->GetTree()->Branch(
          "ReactionConditions", "TReactionConditions", &m_ReactionConditions);
  }

  // Fusion
  else  if (m_ReactionType==Fusion) {
    vector<InputBlock*> blocks=input.GetAllBlocksWithToken("FusionReaction");
    m_BeamName = NPL::ChangeNameToG4Standard(blocks[0]->GetString("Beam"));
    m_TargetNuclei = blocks[0]->GetString("Target");
    m_FusionProduct = blocks[0]->GetString("Product");
    m_FusionExcitation = blocks[0]->GetDouble("ExcitationEnergy","MeV");
    m_active = true;
    // not used
    m_ReactionConditions = new TReactionConditions();
  }
  else {
    m_active = false;
  }
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::BeamReaction::IsApplicable(const G4ParticleDefinition& particleType) {
  if (!m_active)
    return false;

  static std::string particleName;
  particleName = particleType.GetParticleName();
  if(particleName=="neutron") particleName="n1";
  else if(particleName=="e-") particleName="electron";
  if (particleName.find(m_BeamName) != std::string::npos) {
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::BeamReaction::ModelTrigger(const G4FastTrack& fastTrack) {
  //cout<< "--------- MODEL TRIG ---------"<<endl;
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  m_Parent_ID = PrimaryTrack->GetParentID();
  // process reserved to the beam
  if(m_Parent_ID!=0)
    return false;

  G4ThreeVector  V            = PrimaryTrack->GetMomentum().unit();
  G4ThreeVector  P            = fastTrack.GetPrimaryTrackLocalPosition();
  G4VSolid*      solid        = fastTrack.GetPrimaryTrack()
    ->GetVolume()
    ->GetLogicalVolume()
    ->GetSolid();
  double  to_exit      = solid->DistanceToOut(P, V);
  double  to_entrance  = solid->DistanceToOut(P, -V);
  bool is_first = (to_entrance==0);

  if(is_first && m_shoot){
    /* Does occur rarely when event is tangent to the target surface and scatters out
       std::cout << "Something went wrong in beam reaction, m_shoot and is_first variables cannot be true simultaneously" << std::endl;
       std::cout << "m_shoot: " << m_shoot << std::endl;
       std::cout << "rand: " << m_rand << std::endl;
       std::cout << "to_entrance: " << to_entrance << std::endl;
       std::cout << "to_exit: " << to_exit << std::endl;
       std::cout << "length: " << m_length << std::endl;
       std::cout << "step: " << m_StepSize << std::endl;
       std::cout << "Z: " << m_Z << std::endl;
       std::cout << "S: " << m_S << std::endl;
       */
    m_shoot = false;
  }

  if(is_first){
    m_rand = G4RandFlat::shoot(); 
    // random Z in the Volume
    m_Z = m_rand*(to_exit+to_entrance)-0.5*(to_exit+to_entrance);
    // Clear Previous Event
    m_ReactionConditions->Clear();
    m_shoot=true;
  }

  // curviligne coordinate along beam path
  m_S = to_entrance - 0.5*(to_exit+to_entrance); 
  m_length = m_Z-m_S;

  m_StepSize = PrimaryTrack->GetVolume()->GetLogicalVolume()->GetUserLimits()->GetMaxAllowedStep(*PrimaryTrack);

  // If the condition is met, the event is generated
  if (m_shoot && m_length < m_StepSize) {
    if(m_ReactionType==QFS){
      if ( m_QFS.IsAllowed() ) {
        return true;
      }
      else{
        m_shoot=false;
        std::cout << "QFS not allowed" << std::endl;
      }
    }

    else if(m_ReactionType==TwoBody){
      if ( m_Reaction.IsAllowed(PrimaryTrack->GetKineticEnergy()) ) {
        return true;
      } 
      else{
        m_shoot=false;
        std::cout << "Two body reaction not allowed" << std::endl;
      }
    }
    else if(m_ReactionType==Fusion){
      return true;
    }
  }

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void NPS::BeamReaction::DoIt(const G4FastTrack& fastTrack,
    G4FastStep&        fastStep) {

  // std::cout << "DOIT" << std::endl; 
  m_shoot=false;
  m_length=abs(m_length);
  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  G4ThreeVector  pdirection   = PrimaryTrack->GetMomentum().unit();
  G4ThreeVector  localdir     = fastTrack.GetPrimaryTrackLocalDirection();

  G4ThreeVector worldPosition = PrimaryTrack->GetPosition();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4Material*   material = fastTrack.GetPrimaryTrack()
    ->GetVolume()
    ->GetLogicalVolume()
    ->GetMaterial();

  double energy = PrimaryTrack->GetKineticEnergy();
  double speed  = PrimaryTrack->GetVelocity();
  double time   = PrimaryTrack->GetGlobalTime()+m_length/speed;


  double reac_energy = SlowDownBeam (
      PrimaryTrack->GetParticleDefinition(),
      energy,
      m_length,
      material);


  G4ThreeVector ldir = pdirection;
  ldir *= m_length;
  localPosition = localPosition + ldir;

  // Set the end of the step conditions
  fastStep.SetPrimaryTrackFinalKineticEnergyAndDirection(0, pdirection);
  fastStep.SetPrimaryTrackFinalPosition(worldPosition);
  fastStep.SetTotalEnergyDeposited(0);
  fastStep.SetPrimaryTrackFinalTime(time);// FIXME
  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);

  ///////////////////////////////
  // Two-Body Reaction Case /////
  ///////////////////////////////
  if (m_ReactionType==TwoBody) {

    static G4IonTable* IonTable
      = G4ParticleTable::GetParticleTable()->GetIonTable();

    //////Define the kind of particle to shoot////////
    // Particle 3
    G4ParticleDefinition* LightName;
    if(m_Reaction.GetParticle3()->GetName()=="electron"){
      LightName=G4Electron::Definition();
    }
    else{
      int LightZ = m_Reaction.GetParticle3()->GetZ();
      int LightA = m_Reaction.GetParticle3()->GetA();
      if (LightZ == 0 && LightA == 1){
        LightName = G4Neutron::Definition();
      } 
      else {
        if (m_Reaction.GetUseExInGeant4())
          LightName
            = IonTable->GetIon(LightZ, LightA, m_Reaction.GetExcitation3() * MeV);
        else
          LightName = IonTable->GetIon(LightZ, LightA);
      }

    }


    // Particle 4
    G4int HeavyZ = m_Reaction.GetParticle4()->GetZ();
    G4int HeavyA = m_Reaction.GetParticle4()->GetA();

    // Generate the excitation energy if a distribution is given
    m_Reaction.ShootRandomExcitationEnergy();

    // Use to clean up the IonTable in case of the Ex changing at every event
    G4ParticleDefinition* HeavyName;

    if (m_Reaction.GetUseExInGeant4())
      HeavyName
        = IonTable->GetIon(HeavyZ, HeavyA, m_Reaction.GetExcitation4() * MeV);
    else
      HeavyName = IonTable->GetIon(HeavyZ, HeavyA);

    // Set the Energy of the reaction
    m_Reaction.SetBeamEnergy(reac_energy);

    double Beam_theta = pdirection.theta();
    double Beam_phi   = pdirection.phi();

    ///////////////////////////
    ///// Beam Parameters /////
    ///////////////////////////
    m_ReactionConditions->SetBeamParticleName(
        PrimaryTrack->GetParticleDefinition()->GetParticleName());

    m_ReactionConditions->SetBeamReactionEnergy(reac_energy);
    m_ReactionConditions->SetVertexPositionX(localPosition.x());
    m_ReactionConditions->SetVertexPositionY(localPosition.y());
    m_ReactionConditions->SetVertexPositionZ(localPosition.z());

    G4ThreeVector U(1, 0, 0);
    G4ThreeVector V(0, 1, 0);
    G4ThreeVector ZZ(0, 0, 1);
    m_ReactionConditions->SetBeamEmittanceTheta(
        PrimaryTrack->GetMomentumDirection().theta() / deg);
    m_ReactionConditions->SetBeamEmittancePhi(
        PrimaryTrack->GetMomentumDirection().phi() / deg);
    m_ReactionConditions->SetBeamEmittanceThetaX(
        PrimaryTrack->GetMomentumDirection().angle(U) / deg);
    m_ReactionConditions->SetBeamEmittancePhiY(
        PrimaryTrack->GetMomentumDirection().angle(V) / deg);

    //////////////////////////////////////////////////////////
    ///// Build rotation matrix to go from the incident //////
    ///// beam frame to the "world" frame               //////
    //////////////////////////////////////////////////////////


    //   G4ThreeVector col1(cos(Beam_theta) * cos(Beam_phi),
    //   cos(Beam_theta) * sin(Beam_phi),
    //   -sin(Beam_theta));
    //   G4ThreeVector col2(-sin(Beam_phi),
    //   cos(Beam_phi),
    //   0);
    //   G4ThreeVector col3(sin(Beam_theta) * cos(Beam_phi),
    //   sin(Beam_theta) * sin(Beam_phi),
    //   cos(Beam_theta));
    //   G4RotationMatrix BeamToWorld(col1, col2, col3);

    /////////////////////////////////////////////////////////////////
    ///// Angles for emitted particles following Cross Section //////
    ///// Angles are in the beam frame                         //////
    /////////////////////////////////////////////////////////////////

    // Angles
    // Shoot and Set a Random ThetaCM
    m_Reaction.ShootRandomThetaCM();
    double phi = G4RandFlat::shoot() * 2. * pi;

    //////////////////////////////////////////////////
    /////  Momentum and angles from  kinematics  /////
    /////  Angles are in the beam frame          /////
    //////////////////////////////////////////////////
    // Variable where to store results
    double Theta3, Energy3, Theta4, Energy4;

    // Compute Kinematic using previously defined ThetaCM
    m_Reaction.KineRelativistic(Theta3, Energy3, Theta4, Energy4);
    // Momentum in beam frame for light particle
    G4ThreeVector momentum_kine3_beam(sin(Theta3) * cos(phi),
        sin(Theta3) * sin(phi), cos(Theta3));
    // Momentum in World frame //to go from the incident beam frame to the "world"
    // frame
    G4ThreeVector momentum_kine3_world = momentum_kine3_beam;
    momentum_kine3_world.rotate(Beam_theta,
        V); // rotation of Beam_theta on Y axis
    momentum_kine3_world.rotate(Beam_phi, ZZ); // rotation of Beam_phi on Z axis

    // Momentum in beam frame for heavy particle
    G4ThreeVector momentum_kine4_beam(sin(Theta4) * cos(phi + pi),
        sin(Theta4) * sin(phi + pi), cos(Theta4));
    // Momentum in World frame
    G4ThreeVector momentum_kine4_world = momentum_kine4_beam;
    momentum_kine4_world.rotate(Beam_theta,
        V); // rotation of Beam_theta on Y axis
    momentum_kine4_world.rotate(Beam_phi, ZZ); // rotation of Beam_phi on Z axis

    // Emitt secondary
    if (m_Reaction.GetShoot3()) {
      G4DynamicParticle particle3(LightName, momentum_kine3_world, Energy3);
      fastStep.CreateSecondaryTrack(particle3, localPosition, time);
    }

    if (m_Reaction.GetShoot4()) {
      G4DynamicParticle particle4(HeavyName, momentum_kine4_world, Energy4);
      fastStep.CreateSecondaryTrack(particle4, localPosition, time);
    }

    ///////////////////////////////////////
    ///// Emitted particle Parameters /////
    ///////////////////////////////////////
    // Names 3 and 4//
    m_ReactionConditions->SetParticleName(LightName->GetParticleName());
    m_ReactionConditions->SetParticleName(HeavyName->GetParticleName());
    // Angle 3 and 4 //
    m_ReactionConditions->SetTheta(Theta3 / deg);
    m_ReactionConditions->SetTheta(Theta4 / deg);

    m_ReactionConditions->SetPhi(phi / deg);
    if ((phi + pi) / deg > 360)
      m_ReactionConditions->SetPhi((phi - pi) / deg);
    else
      m_ReactionConditions->SetPhi((phi + pi) / deg);

    // Energy 3 and 4 //
    m_ReactionConditions->SetKineticEnergy(Energy3);
    m_ReactionConditions->SetKineticEnergy(Energy4);
    // ThetaCM and Ex//
    m_ReactionConditions->SetThetaCM(m_Reaction.GetThetaCM() / deg);
    m_ReactionConditions->SetExcitationEnergy3(m_Reaction.GetExcitation3());
    m_ReactionConditions->SetExcitationEnergy4(m_Reaction.GetExcitation4());
    // Momuntum X 3 and 4 //
    m_ReactionConditions->SetMomentumDirectionX(momentum_kine3_world.x());
    m_ReactionConditions->SetMomentumDirectionX(momentum_kine4_world.x());
    // Momuntum Y 3 and 4 //
    m_ReactionConditions->SetMomentumDirectionY(momentum_kine3_world.y());
    m_ReactionConditions->SetMomentumDirectionY(momentum_kine4_world.y());
    // Momuntum Z 3 and 4 //
    m_ReactionConditions->SetMomentumDirectionZ(momentum_kine3_world.z());
    m_ReactionConditions->SetMomentumDirectionZ(momentum_kine4_world.z());

  }// end if TwoBodyReaction


  // QFS
  else if (m_ReactionType==QFS ) {

    //////Define the kind of particle to shoot////////
    //    A --> T  ==> B + (c -> T) =>  B + 1 + 2      

    int Light1_Z = m_QFS.GetParticle1()->GetZ();
    int Light1_A = m_QFS.GetParticle1()->GetA();
    int Light2_Z = m_QFS.GetParticle2()->GetZ();
    int Light2_A = m_QFS.GetParticle2()->GetA();

    static G4IonTable* IonTable
      = G4ParticleTable::GetParticleTable()->GetIonTable();

    G4ParticleDefinition* Light1Name;
    G4ParticleDefinition* Light2Name;

    if (Light1_Z == 0 && Light1_A == 1) // neutron is special case
    {
      Light1Name = G4Neutron::Definition();
    } else {
      Light1Name = IonTable->GetIon(Light1_Z, Light1_A);
    }

    if (Light2_Z == 0 && Light2_A == 1) // neutron is special case
    {
      Light2Name = G4Neutron::Definition();
    } else {
      Light2Name = IonTable->GetIon(Light2_Z, Light2_A);
    }

    // Particle B
    G4int Heavy_Z = m_QFS.GetParticleB()->GetZ();
    G4int Heavy_A = m_QFS.GetParticleB()->GetA();

    G4ParticleDefinition* HeavyName;
    if (m_QFS.GetUseExInGeant4())
      HeavyName
        = IonTable->GetIon(Heavy_Z, Heavy_A, m_QFS.GetExcitationB() * MeV);
    else
      HeavyName = IonTable->GetIon(Heavy_Z, Heavy_A);

    // Set the Energy of the reaction
    m_QFS.SetBeamEnergy(reac_energy);

    double Beam_theta = pdirection.theta();
    double Beam_phi   = pdirection.phi();


    /////////////////////////////////////////////////////////////////
    ///// Angles for emitted particles following Cross Section //////
    ///// Angles are in the beam frame                         //////
    /////////////////////////////////////////////////////////////////

    // Shoot and Set a Random ThetaCM
    double costheta = G4RandFlat::shoot() * 2 - 1;
    double theta = acos(costheta);
    double phi = G4RandFlat::shoot() * 2. * pi - pi; //rand in [-pi,pi]
    m_QFS.SetThetaCM(theta);
    m_QFS.SetPhiCM(phi);

    // Shoot and set internal momentum for the removed cluster
    // if a momentum Sigma is given then shoot in 3 indep. Gaussians
    // if input files are given for distributions use them instead
    m_QFS.ShootInternalMomentum();

    //////////////////////////////////////////////////
    /////  Momentum and angles from  kinematics  /////
    //////////////////////////////////////////////////
    // Variable where to store results
    double Theta1, Phi1, TKE1, Theta2, Phi2, TKE2, ThetaB, PhiB, TKEB;

    // Compute Kinematic using previously defined ThetaCM
    m_QFS.KineRelativistic(Theta1, Phi1, TKE1, Theta2, Phi2, TKE2);

    G4ThreeVector U(1, 0, 0);
    G4ThreeVector V(0, 1, 0);
    G4ThreeVector ZZ(0, 0, 1);

    // Momentum in beam and world frame for light particle 1
    G4ThreeVector momentum_kine1_beam(sin(Theta1) * cos(Phi1),
        sin(Theta1) * sin(Phi1), cos(Theta1));
    G4ThreeVector momentum_kine1_world = momentum_kine1_beam;
    momentum_kine1_world.rotate(Beam_theta, V); // rotation of Beam_theta on Y axis
    momentum_kine1_world.rotate(Beam_phi, ZZ); // rotation of Beam_phi on Z axis

    // Momentum in beam and world frame for light particle 2
    G4ThreeVector momentum_kine2_beam(sin(Theta2) * cos(Phi2),
        sin(Theta2) * sin(Phi2), cos(Theta2));
    G4ThreeVector momentum_kine2_world = momentum_kine2_beam;
    momentum_kine2_world.rotate(Beam_theta, V); // rotation of Beam_theta on Y axis
    momentum_kine2_world.rotate(Beam_phi, ZZ); // rotation of Beam_phi on Z axis

    // Momentum in beam and world frame for heavy residual
    //
    //G4ThreeVector momentum_kineB_beam(sin(ThetaB) * cos(PhiB + pi),
    //        sin(ThetaB) * sin(PhiB + pi), cos(ThetaB));
    //G4ThreeVector momentum_kineB_world =  momentum_kineB_beam;
    //momentum_kineB_world.rotate(Beam_theta, V); // rotation of Beam_theta on Y axis
    //momentum_kineB_world.rotate(Beam_phi, ZZ); // rotation of Beam_phi on Z axis

    TLorentzVector* P_A = m_QFS.GetEnergyImpulsionLab_A();
    TLorentzVector* P_B = m_QFS.GetEnergyImpulsionLab_B();

    G4ThreeVector momentum_kineB_beam( P_B->Px(), P_B->Py(), P_B->Pz() );
    momentum_kineB_beam = momentum_kineB_beam.unit();
    TKEB = P_B->Energy() - m_QFS.GetParticleB()->Mass();
    G4ThreeVector momentum_kineB_world =  momentum_kineB_beam;
    momentum_kineB_world.rotate(Beam_theta, V); // rotation of Beam_theta on Y axis
    momentum_kineB_world.rotate(Beam_phi, ZZ); // rotation of Beam_phi on Z axis

    ThetaB = P_B->Angle(P_A->Vect());
    if (ThetaB < 0) ThetaB += M_PI;
    PhiB = M_PI + P_B->Vect().Phi(); 
    if (fabs(PhiB) < 1e-6) PhiB = 0;


    // Emitt secondary
    if (m_QFS.GetShoot1()) {
      G4DynamicParticle particle1(Light1Name, momentum_kine1_world, TKE1);
      fastStep.CreateSecondaryTrack(particle1, localPosition, time);
    }

    if (m_QFS.GetShoot2()) {
      G4DynamicParticle particle2(Light2Name, momentum_kine2_world, TKE2);
      fastStep.CreateSecondaryTrack(particle2, localPosition, time);
    }
    if (m_QFS.GetShootB()) {
      if(m_QFS.GetDeexcitation()){
        double P_B2 = P_B->Px()*P_B->Px() + P_B->Py()*P_B->Py() + P_B->Pz()*P_B->Pz();
        double scaling = 1;
        if(P_B2>0.0){
          double P_Bnew2 = TKEB*TKEB + 2 * TKEB * m_QFS.GetParticleB()->Mass(); 
          scaling = sqrt(P_Bnew2)/sqrt(P_B2);
        }
        G4LorentzVector G4HeavyMomentum;
        G4HeavyMomentum.setE(P_B->Energy());
        G4HeavyMomentum.setPx(P_B->X());
        G4HeavyMomentum.setPy(P_B->Y());
        G4HeavyMomentum.setPz(P_B->Z());
        G4Fragment HeavyFragment(HeavyName->GetAtomicMass(), HeavyName->GetAtomicNumber(), G4HeavyMomentum);
        G4ReactionProductVector* HeavyDeexcitation = ABLA->DeExcite(HeavyFragment);

        unsigned int sizeFrag = HeavyDeexcitation->size();
        G4ThreeVector Momsum(0,0,0);
        for(unsigned int i = 0 ; i < sizeFrag ; i++){
          const G4ParticleDefinition* HeavyFragName = HeavyDeexcitation->at(i)->GetDefinition();
          G4ThreeVector HeavyFragMomentum = HeavyDeexcitation->at(i)->GetMomentum();
          double TKEFrag = HeavyDeexcitation->at(i)->GetKineticEnergy();
          double HeavyFragModule = sqrt(TKEFrag*TKEFrag + 2*TKEFrag*HeavyDeexcitation->at(i)->GetMass());
          G4ThreeVector HeavyFragMomentumDirection(HeavyFragMomentum.x()/HeavyFragModule,
              HeavyFragMomentum.y()/HeavyFragModule,
              HeavyFragMomentum.z()/HeavyFragModule); 
          G4DynamicParticle particleFrag(HeavyFragName, HeavyFragMomentumDirection, TKEFrag);
          fastStep.CreateSecondaryTrack(particleFrag, localPosition, time);
        }
      }
      else{
        G4DynamicParticle particleB(HeavyName, momentum_kineB_world, TKEB);
        fastStep.CreateSecondaryTrack(particleB, localPosition, time);
      }
    }

    ///////////////////////////////////
    ///// Reaction Condition Save /////
    ///////////////////////////////////
    m_ReactionConditions->SetBeamParticleName(
        PrimaryTrack->GetParticleDefinition()->GetParticleName());
    m_ReactionConditions->SetBeamReactionEnergy(reac_energy);
    m_ReactionConditions->SetVertexPositionX(localPosition.x());
    m_ReactionConditions->SetVertexPositionY(localPosition.y());
    m_ReactionConditions->SetVertexPositionZ(localPosition.z());
    m_ReactionConditions->SetBeamEmittanceTheta(
        PrimaryTrack->GetMomentumDirection().theta() / deg);
    m_ReactionConditions->SetBeamEmittancePhi(
        PrimaryTrack->GetMomentumDirection().phi() / deg);
    m_ReactionConditions->SetBeamEmittanceThetaX(
        PrimaryTrack->GetMomentumDirection().angle(U) / deg);
    m_ReactionConditions->SetBeamEmittancePhiY(
        PrimaryTrack->GetMomentumDirection().angle(V) / deg);

    // Names 1,2 and B//
    m_ReactionConditions->SetParticleName(Light1Name->GetParticleName());
    m_ReactionConditions->SetParticleName(Light2Name->GetParticleName());
    m_ReactionConditions->SetParticleName(HeavyName->GetParticleName());
    // Angle 1,2 and B //
    m_ReactionConditions->SetTheta(Theta1 / deg);
    m_ReactionConditions->SetTheta(Theta2 / deg);
    m_ReactionConditions->SetTheta(ThetaB / deg);
    m_ReactionConditions->SetPhi(Phi1 / deg);
    m_ReactionConditions->SetPhi(Phi2 / deg);
    m_ReactionConditions->SetPhi(PhiB / deg);
    // Energy 1,2 and B //
    m_ReactionConditions->SetKineticEnergy(TKE1);
    m_ReactionConditions->SetKineticEnergy(TKE2);
    m_ReactionConditions->SetKineticEnergy(TKEB);
    // ThetaCM, Ex and Internal Momentum of removed cluster//
    m_ReactionConditions->SetThetaCM(m_QFS.GetThetaCM() / deg);
    m_ReactionConditions->SetInternalMomentum(m_QFS.GetInternalMomentum());
    //m_ReactionConditions->SetExcitationEnergy3(m_QFS.GetExcitation3());
    m_ReactionConditions->SetExcitationEnergy4(m_QFS.GetExcitationB());
    // Momuntum X 1,2 and B //
    m_ReactionConditions->SetMomentumDirectionX(momentum_kine1_world.x());
    m_ReactionConditions->SetMomentumDirectionX(momentum_kine2_world.x());
    m_ReactionConditions->SetMomentumDirectionX(momentum_kineB_world.x());
    // Momuntum Y 1,2 and B //
    m_ReactionConditions->SetMomentumDirectionY(momentum_kine1_world.y());
    m_ReactionConditions->SetMomentumDirectionY(momentum_kine2_world.y());
    m_ReactionConditions->SetMomentumDirectionY(momentum_kineB_world.y());
    // Momuntum Z 1,2 and B //
    m_ReactionConditions->SetMomentumDirectionZ(momentum_kine1_world.z());
    m_ReactionConditions->SetMomentumDirectionZ(momentum_kine2_world.z());
    m_ReactionConditions->SetMomentumDirectionZ(momentum_kineB_world.z());

  }// end QFS

  ////////////////////
  //  Fusion Case   //
  ////////////////////
  else if(m_ReactionType==Fusion){
    // Set the end of the step conditions
    fastStep.SetPrimaryTrackFinalKineticEnergyAndDirection(0, pdirection);
    fastStep.SetPrimaryTrackFinalPosition(worldPosition);
    fastStep.SetTotalEnergyDeposited(0);
    fastStep.SetPrimaryTrackFinalTime(time);// FIXME
    fastStep.KillPrimaryTrack();
    fastStep.SetPrimaryTrackPathLength(0.0);

    static G4IonTable* IonTable
      = G4ParticleTable::GetParticleTable()->GetIonTable();
    //////Define the kind of particle to shoot////////
    G4ParticleDefinition* Product;
    NPL::Particle N(m_FusionProduct);  
    N.SetExcitationEnergy(m_FusionExcitation);
    NPL::Particle T(m_TargetNuclei);  
    int PZ = N.GetZ();
    int PA = N.GetA();
    Product = IonTable->GetIon(PZ, PA, m_FusionExcitation*MeV);
    // setup the daugter
    /////////////////////////////////////////////////////////////////////////
    TVector3 BeamP = NPS::ConvertVector(PrimaryTrack->GetMomentum()); 

    TLorentzVector BeamLV;
    BeamLV.SetVectM(BeamP,N.Mass()*MeV);
    TLorentzVector TargetLV;
    TargetLV.SetVectM(TVector3(0,0,0),T.Mass()*MeV);
    TLorentzVector TotalLV=BeamLV+TargetLV;

    // energy lost in the fusion process to be removed to the total energy 
    // Total Available Ek = Initial Ek + (InitialMass-FinalMass)
    double KineAvailable= TotalLV.Et()+ (TotalLV.Mag()-N.Mass());
    G4ThreeVector momentum_dir = NPS::ConvertVector(TotalLV.Vect().Unit());

    //////FIXME Unsure of this part
    // Randomize Phi after the reaction
    double Phi = CLHEP::RandFlat::shoot()*2*M_PI ;
    momentum_dir.rotate(Phi,PrimaryTrack->GetMomentum());

    G4DynamicParticle particle(Product, momentum_dir, KineAvailable);
    fastStep.CreateSecondaryTrack(particle, localPosition, time);
  }// end fusion

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Return the slow down beam in given thickness of material
double NPS::BeamReaction::SlowDownBeam(const G4ParticleDefinition* Beam, 
    double IncidentEnergy, 
    double Thickness,
    G4Material* Material){


  if(Beam->GetParticleName()=="neutron"){
    return IncidentEnergy;
  }

  double dedx,de;
  int NbSlice = 100;
  static G4EmCalculator emCalculator;

  if(Thickness!=0){
    for (G4int i = 0; i < NbSlice; i++){
      dedx = emCalculator.ComputeTotalDEDX(IncidentEnergy, Beam, Material);
      de   = dedx * Thickness / NbSlice;
      IncidentEnergy -= de;
      if(IncidentEnergy<0){
        IncidentEnergy = 0;
        break;
      }
    }
  }

  return IncidentEnergy;

}

