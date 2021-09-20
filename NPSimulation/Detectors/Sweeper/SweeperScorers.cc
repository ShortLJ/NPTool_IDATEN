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
 *  File old the scorer specific to the Sharc Detector                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new type of scorer is aim to become the standard for DSSD,SSSD and   *
 * PAD detector (any Silicon Detector)                                       *
 *****************************************************************************/
#include "SweeperScorers.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"

using namespace SweeperScorers;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Sweeper::PS_Sweeper(G4String name, vector<G4int> NestingLevel,G4int depth)
	:G4VPrimitiveScorer(name, depth),HCID(-1){
	m_NestingLevel = NestingLevel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Sweeper::~PS_Sweeper(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Sweeper::ProcessHits(G4Step* aStep, G4TouchableHistory*){
	G4double edep = aStep->GetTotalEnergyDeposit();
        bool fGeomDetectedNeutron=true;
	
	// if(edep <= 0) {
	//   if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && aStep->GetTrack()->GetParentID()==1 && aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge()==0)cout<<"GeomDetected"<<endl;
	//   return false; }
	
	//G4bool isFragment = (aStep->GetTrack()->GetParentID() ==1 && aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge()!=0) ? true : false;
	// if(isFragment) {
	//   aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
	//   //cout<<"det frag   "<<  aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_NestingLevel[0])<<endl;
	//   return false;
	//}
	//cout<<"IS NEUTRON : " << aStep->GetTrack()->GetParentID() <<" "<< aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge()<<endl;
	// Contain Energy, Time + as many copy number as nested volume
	unsigned int mysize = m_NestingLevel.size();
	StepInfo step_info;
	auto& Infos = step_info.Infos;
	Infos.resize(7+mysize);
	Infos[0] = edep;
	Infos[1] = aStep->GetPreStepPoint()->GetGlobalTime();
    
    
	// Interaction coordinates (used to fill the InteractionCoordinates branch)
	m_Position  = aStep->GetPreStepPoint()->GetPosition();
    
	Infos[2] = m_Position.x();
	Infos[3] = m_Position.y();
	Infos[4] = m_Position.z();
	Infos[5] = m_Position.theta();
	Infos[6] = m_Position.phi();
	
	for(unsigned int i = 0 ; i < mysize ; i++){
		Infos[i+7] = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_NestingLevel[i]);
	}
    
	m_Index = 0 ;
	G4int multiplier = 1;
	for(unsigned int i = 0 ; i < mysize ; i++){
		m_Index+= aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_NestingLevel[i])*multiplier;
		multiplier*=10;
	}
    
	// Check if the particle has interact before, if yes, add up the energies.
	//       inline std::map<G4int,T*>* GetMap() const
	//
	map<G4int, StepInfo*>::iterator it;
	it= EvtMap->GetMap()->find(m_Index);
	if(it!=EvtMap->GetMap()->end()){
		StepInfo* dummy = it->second;
		Infos[0]+=dummy->Infos[0];
	}

	// Record particle information
	{
		const G4ParticleDefinition* p = aStep->GetTrack()->GetParticleDefinition();
		step_info.Particle_Name = p->GetParticleName();
		if(p->GetParticleName() == "e-") {
			step_info.Particle_Z = -1;
			step_info.Particle_A =  0;
		} else if(p->GetParticleName() == "e+") {
			step_info.Particle_Z = +1;
			step_info.Particle_A =  0;				
		}	else {
			step_info.Particle_Z = p->GetPDGCharge();
			step_info.Particle_A = p->GetAtomicNumber();
		}
	}

        
	// // Record process information
	if(aStep->GetTrack()->GetCreatorProcess()) {
		step_info.ProcessName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	} else {
		step_info.ProcessName = "NULL";
	}
	
	EvtMap->set(m_Index, step_info);
	return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Sweeper::Initialize(G4HCofThisEvent* HCE){
	EvtMap = new NPS::HitsMap<StepInfo>(GetMultiFunctionalDetector()->GetName(), GetName());
	if (HCID < 0) {
		HCID = GetCollectionID(0);
	}
	HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Sweeper::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Sweeper::clear(){   
	EvtMap->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Sweeper::DrawAll(){
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Sweeper::PrintAll(){
	G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
	G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
	G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}

