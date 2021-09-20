/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : February 2013                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  File old the scorer to record Hit energy,time and position               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "InteractionScorers.hh"
#include "G4UnitsTable.hh"
using namespace InteractionScorers ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<InteractionData>::iterator InteractionDataVector::find(const unsigned int& index){
  for(vector<InteractionData>::iterator it= m_Data.begin()  ; it !=m_Data.end() ; it++){
    if((*it).GetIndex()==index)
      return it;
  }
  return m_Data.end();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Interactions::PS_Interactions(G4String name,TInteractionCoordinates* Inter, int depth)  :G4VPrimitiveScorer(name, depth){
  m_Level = depth;
  m_InterractionCoordinates=Inter;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Interactions::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  static G4StepPoint* point;
  point = aStep->GetPreStepPoint();
  t_Energy = aStep->GetTotalEnergyDeposit();
  t_Time = point->GetGlobalTime();
  t_Position  = point->GetPosition();
  t_Index = aStep->GetTrack()->GetTrackID(); 
  t_ParticleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName(); 
  t_A = aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass();
  t_Z = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
  t_Mass = aStep->GetTrack()->GetDynamicParticle()->GetMass();
  t_Charge = aStep->GetTrack()->GetDynamicParticle()->GetCharge();
  double KineticEnergy = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
  t_Brho = sqrt(KineticEnergy * KineticEnergy + 2 * KineticEnergy * t_Mass)
    / (c_light * t_Charge);

  // add it to check the theta of momentum
  // MOMENT  = aStep->GetPreStepPoint()->GetMomentumDirection();

  vector<InteractionData>::iterator it;
  it = m_DataVector.find(t_Index); 
  if(it!=m_DataVector.end())
    it->Add(t_Energy);
  else
  { m_DataVector.Set(t_Index,t_Energy,t_Time,t_Position.x(),t_Position.y(),t_Position.z(), /*MOMENT*/ t_Position.theta(),t_Position.phi(), t_ParticleName, t_A, t_Z, t_Mass, t_Charge, t_Brho);

  }

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Interactions::Initialize(G4HCofThisEvent*){
  // Clear is called by EventAction
  clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Interactions::EndOfEvent(G4HCofThisEvent*){
  unsigned int size = m_DataVector.size();

  for(unsigned int i = 0 ; i < size ; i++){
    int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    int Index = t_Index + eventID; 
    m_InterractionCoordinates->SetInteraction(Index,m_DataVector[i]->GetEnergy(),m_DataVector[i]->GetTime(),
        m_DataVector[i]->GetPositionX(),m_DataVector[i]->GetPositionY(),m_DataVector[i]->GetPositionZ(),
        m_DataVector[i]->GetTheta()/deg,m_DataVector[i]->GetPhi()/deg, m_DataVector[i]->GetParticleName(),
        m_DataVector[i]->GetA(), m_DataVector[i]->GetZ(), m_DataVector[i]->GetMass(), m_DataVector[i]->GetCharge(), 
        m_DataVector[i]->GetBrho());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Interactions::clear(){
  m_DataVector.clear();
  m_InterractionCoordinates->Clear();
}
