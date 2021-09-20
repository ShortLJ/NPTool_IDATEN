/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This event Generator is used to simulated Radioactive Ion beam           *
 *  Emmitance, Energy spread and dispersion are taken into account           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++
#include <limits>

// G4 header
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
// G4 headers including CLHEP headers
// for generating random numbers
#include "Randomize.hh"

// NPSimulation header
#include "EventGeneratorBeam.hh"
#include "Particle.hh"
// NPL header
#include "RootOutput.h"

using namespace CLHEP;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorBeam::EventGeneratorBeam(){
  m_ParticleStack  = ParticleStack::getInstance();
  m_Target         = NULL    ;
  m_particle       = NULL    ;
  m_Beam           = new NPL::Beam();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorBeam::~EventGeneratorBeam(){
  delete m_Beam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void   EventGeneratorBeam::SetTarget(Target* Target){
  if(Target!=0)
    m_Target = Target;
  
  // Set the target parameter for the internal event generator of m_Beam
  m_Beam->SetTargetSize(m_Target->GetTargetRadius());
  m_Beam->SetTargetAngle(m_Target->GetTargetAngle());
  m_Beam->SetTargetThickness(m_Target->GetTargetThickness());
  m_Beam->SetTargetZ(m_Target->GetTargetZ());
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorBeam::ReadConfiguration(NPL::InputParser parser){
  m_Beam->ReadConfigurationFile(parser);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorBeam::GenerateEvent(G4Event* anEvent){

  if( anEvent->GetEventID()==0){
    // Define the particle to be shoot
    if(m_Beam->GetZ()==0 &&  m_Beam->GetA()==1)
      m_particle = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    else if(m_Beam->GetName()=="electron")
      m_particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
   else if(m_Beam->GetName()=="gamma")
      m_particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

   else
      m_particle = 
        G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(m_Beam->GetZ(), m_Beam->GetA() ,m_Beam->GetExcitationEnergy());

  }
  
  ///////////////////////////////////////////////////////////////////////
  ///// Calculate the incident beam direction as well as the vertex /////
  ///// of interaction in target and Energy Loss of the beam within /////
  ///// the target.                                                 /////
  ///////////////////////////////////////////////////////////////////////
  //G4ThreeVector InterCoord;
  
  double InitialBeamEnergy, x0, y0, z0, Beam_thetaX, Beam_phiY;

  m_Beam->GenerateRandomEvent(InitialBeamEnergy, x0, y0, z0, Beam_thetaX, Beam_phiY);
  //Set the direction cosines: alpha=90-Beam_thetaX, beta=90-Beam_phiY
  double Xdir = sin(Beam_thetaX); // cos(90-x) = sin(x)
  double Ydir = sin(Beam_phiY); 
  double Zdir = sqrt(1-Xdir*Xdir-Ydir*Ydir); // alpha^2 + beta^2 + gamma^2 = 1
  G4ThreeVector BeamDir(Xdir,Ydir,Zdir);
  G4ThreeVector BeamPos(x0,y0,z0);
  ///////////////////////////////////////////////////////
  ///// Add the Beam particle to the particle Stack /////
  ///////////////////////////////////////////////////////
  Particle BeamParticle( m_particle,
                        InitialBeamEnergy,
                        InitialBeamEnergy,
                        BeamDir.unit(),
                        BeamPos,
                        1);
  
  m_ParticleStack->AddBeamParticleToStack(BeamParticle);

  }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorBeam::InitializeRootOutput(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
