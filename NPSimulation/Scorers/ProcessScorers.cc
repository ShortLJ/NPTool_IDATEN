/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. MORFOUACE  contact address: pierre.morfouace2@cea.fr  *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Scorer specific to the porcesses occuring in the simulation              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * The scorer hold the processes name                                        *
 *****************************************************************************/
#include "ProcessScorers.hh"
#include "G4UnitsTable.hh"
#include "G4SteppingManager.hh"
using namespace ProcessScorers ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Process::PS_Process(G4String name,G4int depth)
  :G4VPrimitiveScorer(name, depth){
    //m_NestingLevel = NestingLevel;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Process::~PS_Process(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Process::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  // Contain Process Name
  G4String processname;
  if(aStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
    processname = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    
    t_processname.push_back(processname);
    t_processtime.push_back(aStep->GetPreStepPoint()->GetGlobalTime());
  }

  /*G4int trackID;
  //G4int parentID;
  G4String volumename, particlename;
  //G4double step_length;
  G4double track_kineE;
  if(aStep->GetTrack()->GetNextVolume() != 0){
    trackID = aStep->GetTrack()->GetTrackID();
    //parentID = aStep->GetTrack()->GetParentID();
    volumename = aStep->GetTrack()->GetVolume()->GetName();
    particlename = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    //step_length = aStep->GetTrack()->GetStepLength();
    track_kineE = aStep->GetTrack()->GetKineticEnergy();
    
    if(volumename=="av_1_impr_1_gas_box_logic_pv_0"){
      if(processname=="hadElastic")t_FC_process.push_back(1);
      if(processname=="neutronInelastic")t_FC_process.push_back(2);
    }

    if(HasBeenTracked[trackID]==0){
      if(particlename=="gamma"){ 
        HasBeenTracked[trackID]=1;
        //cout << trackID  << " " << track_kineE << endl;
        t_gamma_energy.push_back(track_kineE);
      }
      if(particlename=="proton"){
        HasBeenTracked[trackID]=1;
        if(track_kineE>0.1){
          t_proton_energy.push_back(track_kineE);
          t_proton_time.push_back(aStep->GetPreStepPoint()->GetGlobalTime());
        }
      }
    }
  }*/

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::Initialize(G4HCofThisEvent*){
  clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::clear(){
  t_processname.clear();
  t_processtime.clear();
  t_gamma_energy.clear();
  t_proton_energy.clear();
  t_proton_time.clear();
  t_FC_process.clear();
  
  for(unsigned int i=0; i<100; i++){
    HasBeenTracked[i] = 0;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::DrawAll(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::PrintAll(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

