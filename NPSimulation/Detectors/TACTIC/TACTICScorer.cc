//#define USE_Garfield //only use if compiled with Garfield

#include "TACTICScorer.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "TACTIC.hh"

#ifdef USE_Garfield
#include "GARFDRIFT.h"
#endif

double excess;

using namespace TACTICScorer;

Gas_Scorer::Gas_Scorer(G4String name,G4int Level,G4double ScorerLength,G4int NumberOfSegments, G4int depth, G4double p0, G4double p1, G4double p2, G4double p3,string Shape) //what do level and depth do?       
:G4VPrimitiveScorer(name, depth),HCID(-1){
  m_ScorerLength = ScorerLength;
  m_NumberOfSegments = NumberOfSegments;
  m_SegmentLength = m_ScorerLength / m_NumberOfSegments;
  m_Level = Level;
  m_Position = G4ThreeVector(-1000,-1000,-1000);
  m_SegmentNumber = -1;
  m_Index = -1;
  m_p0 = p0;
  m_p1 = p1;
  m_p2 = p2;
  m_p3 = p3;
  m_Shape = Shape;
}

Gas_Scorer::~Gas_Scorer(){}

G4bool Gas_Scorer::ProcessHits(G4Step* aStep, G4TouchableHistory*){

  G4double* Infos = new G4double[15];
  //bool first_step = true;
  m_Position  = aStep->GetPreStepPoint()->GetPosition();
  
  Infos[0] = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  Infos[1] = aStep->GetTrack()->GetTrackID();

  Infos[2] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();;
  
  Infos[3] = aStep->GetPreStepPoint()->GetGlobalTime();
  Infos[4] = aStep->GetPreStepPoint()->GetKineticEnergy();
  Infos[5] = aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit();

  m_SegmentNumber = (int)((m_Position.z() + m_ScorerLength / 2.) / m_SegmentLength ) + 1; //Pad number 
  Infos[6] = m_SegmentNumber;
  //prepad = Infos[6]; 
  Infos[7] = m_Position.z();
  if(m_Shape == "Cylindrical") Infos[8] = pow(pow(m_Position.x(),2) + pow(m_Position.y(),2),0.5); //R
  if(m_Shape == "Long_Chamber") Infos[8] = m_Position.y();
  Infos[9] = aStep->GetTrack()->GetVertexPosition()[2];
  Infos[10] = aStep->GetTrack()->GetVertexKineticEnergy();
  G4ThreeVector p_vec = aStep->GetTrack()->GetVertexMomentumDirection();
  Infos[11] = acos(p_vec[2]/pow(pow(p_vec[0],2)+pow(p_vec[1],2)+pow(p_vec[2],2),0.5))/deg; //angle relative to z axis (theta);   
  Infos[12] = aStep->GetTrack()->GetTrackLength();
  Infos[13] = m_p0 + m_p1*Infos[8] + m_p2*Infos[8]*Infos[8] + m_p3*Infos[8]*Infos[8]*Infos[8];
  
  //Infos[14] = excess;
  
  G4ThreeVector delta_Position = aStep->GetDeltaPosition();

  m_DetectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  m_Index = m_DetectorNumber * 1e3 + m_SegmentNumber * 1e6;
  
  if(isnan(Infos[10])) {
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    return 0;
  }

  if(aStep->IsFirstStepInVolume() == true) excess = 0.;
  
  map<G4int, G4double**>::iterator it;
  it= EvtMap->GetMap()->find(m_Index);
  if(it!=EvtMap->GetMap()->end()){
    G4double* dummy = *(it->second);
    if(Infos[1]==dummy[1]) Infos[5]+=dummy[5]; //accumulate ionisation energy deposit to get total accross pad
    delete dummy;
  }
    
#ifdef USE_Garfield

  Infos[14] = GARFDRIFT(((aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit())/eV+excess), Infos[3], m_Position/cm, delta_Position/cm, Infos[8]/cm, Infos[6], Infos[2], m_ScorerLength/cm, m_SegmentLength/cm, Infos[0], (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit())/eV, m_Shape)*eV;
  /*  
  file.open("excess_test.dat",std::ios::app);
  file << Infos[6] << "\t"  << "\t" <<  aStep->IsFirstStepInVolume() << "\t" << excess  << "\t" << Infos[14]/eV << "\t" << (int)((((aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit())/eV+excess) / 41.1)*0.01) << "\t" <<  Infos[8] <<  endl;
  file.close();
  */
  excess = Infos[14]/eV;
 
#endif

  EvtMap->set(m_Index, Infos);

  return TRUE;

}

void Gas_Scorer::Initialize(G4HCofThisEvent* HCE){
    EvtMap = new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
    if (HCID < 0) {
        HCID = GetCollectionID(0);
    }
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void Gas_Scorer::EndOfEvent(G4HCofThisEvent*){}

void Gas_Scorer::clear(){
    std::map<G4int, G4double**>::iterator    MapIterator;
    for (MapIterator = EvtMap->GetMap()->begin() ; MapIterator != EvtMap->GetMap()->end() ; MapIterator++){
        delete *(MapIterator->second);
    }

    EvtMap->clear();
}


void Gas_Scorer::DrawAll(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                              

void Gas_Scorer::PrintAll(){
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
    G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}
