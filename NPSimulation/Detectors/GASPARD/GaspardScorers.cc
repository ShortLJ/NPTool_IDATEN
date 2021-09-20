/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 11/07/09                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: This class holds all the scorers needed by the                *
 *             GaspardTracker*** objects.                                    *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// G4 headers
#include "G4UnitsTable.hh"

// NPTool headers
#include "ObsoleteGeneralScorers.hh"
#include "GaspardScorers.hh"

#include "GaspardTrackerDummyShape.hh"
#include "GaspardTrackerSquare.hh"
#include "GaspardTrackerTrapezoid.hh"
#include "GaspardTrackerAnnular.hh"
using namespace GPDSQUARE;
using namespace GPDTRAP;
using namespace GPDANNULAR;
using namespace GPDDUMMYSHAPE;

using namespace GPDSCORERS;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Added by Adrien MATTA:
// Those Scorer use TrackID as map index. This way ones can rebuild energy deposit,
// time of flight or position,... particle by particle for each event. Because standard
// scorer provide by G4 don't work this way but using a global ID for each event you should
// not use those scorer with some G4 provided ones or being very carefull doing so.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// FirstStage Energy Scorer (deal with multiple particle hit)
GPDScorerFirstStageEnergy::GPDScorerFirstStageEnergy(G4String name, G4String volumeName, G4int depth)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_VolumeName = volumeName;
}

GPDScorerFirstStageEnergy::~GPDScorerFirstStageEnergy()
{
}

G4bool GPDScorerFirstStageEnergy::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, m_VolumeName);

   // get energy
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index = aStep->GetTrack()->GetTrackID();
   EvtMap->add(DetNbr + index, edep);
   return TRUE;
}

void GPDScorerFirstStageEnergy::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerFirstStageEnergy::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerFirstStageEnergy::clear()
{
   EvtMap->clear();
}

void GPDScorerFirstStageEnergy::DrawAll()
{
}

void GPDScorerFirstStageEnergy::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
   G4cout << " PrimitiveScorer " << GetName() << G4endl;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl;
   std::map<G4int, G4double*>::iterator itr = EvtMap->GetMap()->begin();
   for (; itr != EvtMap->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first
      << "  energy deposit: " << G4BestUnit(*(itr->second), "Energy")
      << G4endl;
   }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// SecondStage Energy Scorer (deal with multiple particle hit)
GPDScorerSecondStageEnergy::GPDScorerSecondStageEnergy(G4String name, G4String volumeName, G4int depth)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_VolumeName = volumeName;
}

GPDScorerSecondStageEnergy::~GPDScorerSecondStageEnergy()
{
}

G4bool GPDScorerSecondStageEnergy::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, m_VolumeName);

   // get energy
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index = aStep->GetTrack()->GetTrackID();
   EvtMap->add(DetNbr + index, edep);
   return TRUE;
}

void GPDScorerSecondStageEnergy::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerSecondStageEnergy::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerSecondStageEnergy::clear()
{
   EvtMap->clear();
}

void GPDScorerSecondStageEnergy::DrawAll()
{
}

void GPDScorerSecondStageEnergy::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
   G4cout << " PrimitiveScorer " << GetName() << G4endl;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl;
   std::map<G4int, G4double*>::iterator itr = EvtMap->GetMap()->begin();
   for (; itr != EvtMap->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first
      << "  energy deposit: " << G4BestUnit(*(itr->second), "Energy")
      << G4endl;
   }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ThirdStage Energy Scorer (deal with multiple particle hit)
GPDScorerThirdStageEnergy::GPDScorerThirdStageEnergy(G4String name, G4String volumeName, G4int depth)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_VolumeName = volumeName;
}

GPDScorerThirdStageEnergy::~GPDScorerThirdStageEnergy()
{
}

G4bool GPDScorerThirdStageEnergy::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, m_VolumeName);

   // get energy
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index = aStep->GetTrack()->GetTrackID();
   EvtMap->add(DetNbr + index, edep);
   return TRUE;
}

void GPDScorerThirdStageEnergy::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerThirdStageEnergy::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerThirdStageEnergy::clear()
{
   EvtMap->clear();
}

void GPDScorerThirdStageEnergy::DrawAll()
{
}

void GPDScorerThirdStageEnergy::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
   G4cout << " PrimitiveScorer " << GetName() << G4endl;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl;
   std::map<G4int, G4double*>::iterator itr = EvtMap->GetMap()->begin();
   for (; itr != EvtMap->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first
      << "  energy deposit: " << G4BestUnit(*(itr->second), "Energy")
      << G4endl;
   }
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// FirstStage Front Strip position Scorer for DummyShape geometry
GPDScorerFirstStageFrontStripDummyShape::GPDScorerFirstStageFrontStripDummyShape(G4String name, G4int depth, G4int NumberOfStrip)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_NumberOfStrip    = NumberOfStrip  ;
}

GPDScorerFirstStageFrontStripDummyShape::~GPDScorerFirstStageFrontStripDummyShape()
{
}

G4bool GPDScorerFirstStageFrontStripDummyShape::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, "GPDDummyShape");

   // get front strip number
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
//   G4cout << "POS world: " << POS << G4endl;
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);
//   G4cout << "POS local: " << POS << G4endl;

   G4double StripPitch = GPDDUMMYSHAPE::FirstStageFace / m_NumberOfStrip;

   G4double temp = (POS(0) + GPDDUMMYSHAPE::FirstStageFace / 2.) / StripPitch   ;
   G4int X = int(temp) + 1 ;
//   G4cout << "strip X: " << X << G4endl;

   //Rare case where particle is close to edge of silicon plan
   if (X == m_NumberOfStrip+1) X = m_NumberOfStrip;
   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index =  aStep->GetTrack()->GetTrackID();
   EvtMap->set(DetNbr + index, X);
   return TRUE;
}

void GPDScorerFirstStageFrontStripDummyShape::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4int>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerFirstStageFrontStripDummyShape::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerFirstStageFrontStripDummyShape::clear()
{
   EvtMap->clear();
}

void GPDScorerFirstStageFrontStripDummyShape::DrawAll()
{
}

void GPDScorerFirstStageFrontStripDummyShape::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
   G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// FirstStage Back Strip position Scorer for DummyShape geometry
GPDScorerFirstStageBackStripDummyShape::GPDScorerFirstStageBackStripDummyShape(G4String name, G4int depth, G4int NumberOfStrip)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_NumberOfStrip    = NumberOfStrip  ;
}

GPDScorerFirstStageBackStripDummyShape::~GPDScorerFirstStageBackStripDummyShape()
{
}

G4bool GPDScorerFirstStageBackStripDummyShape::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, "GPDDummyShape");

   // get back strip number
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double StripPitch = GPDDUMMYSHAPE::FirstStageFace / m_NumberOfStrip;

   G4double temp = (POS(1) + GPDDUMMYSHAPE::FirstStageFace / 2.) / StripPitch   ;
   G4int X = int(temp) + 1 ;
//   G4cout << "strip Y: " << X << G4endl;

   //Rare case where particle is close to edge of silicon plan
   if (X == m_NumberOfStrip+1) X = m_NumberOfStrip;
   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index =  aStep->GetTrack()->GetTrackID();
   EvtMap->set(DetNbr + index, X);
   return TRUE;
}

void GPDScorerFirstStageBackStripDummyShape::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4int>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerFirstStageBackStripDummyShape::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerFirstStageBackStripDummyShape::clear()
{
   EvtMap->clear();
}

void GPDScorerFirstStageBackStripDummyShape::DrawAll()
{
}

void GPDScorerFirstStageBackStripDummyShape::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
   G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// FirstStage Front Strip position Scorer for Square geometry
GPDScorerFirstStageFrontStripSquare::GPDScorerFirstStageFrontStripSquare(G4String name, G4int depth, G4int NumberOfStrip)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_NumberOfStrip    = NumberOfStrip  ;
}

GPDScorerFirstStageFrontStripSquare::~GPDScorerFirstStageFrontStripSquare()
{
}

G4bool GPDScorerFirstStageFrontStripSquare::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, "GPDSquare");

   // get front strip
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double StripPitch = GPDSQUARE::FirstStageFace / m_NumberOfStrip;

   G4double temp = (POS(0) + GPDSQUARE::FirstStageFace / 2.) / StripPitch   ;
   G4int X = int(temp) + 1 ;
   //Rare case where particle is close to edge of silicon plan
   if (X == 129) X = 128;
   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index =  aStep->GetTrack()->GetTrackID();
   EvtMap->set(DetNbr + index, X);
   return TRUE;
}

void GPDScorerFirstStageFrontStripSquare::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4int>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerFirstStageFrontStripSquare::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerFirstStageFrontStripSquare::clear()
{
   EvtMap->clear();
}

void GPDScorerFirstStageFrontStripSquare::DrawAll()
{
}

void GPDScorerFirstStageFrontStripSquare::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
   G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// FirstStage Back Strip position Scorer for Square geometry
GPDScorerFirstStageBackStripSquare::GPDScorerFirstStageBackStripSquare(G4String name, G4int depth, G4int NumberOfStrip)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_NumberOfStrip    = NumberOfStrip  ;
}

GPDScorerFirstStageBackStripSquare::~GPDScorerFirstStageBackStripSquare()
{
}

G4bool GPDScorerFirstStageBackStripSquare::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, "GPDSquare");

   // get back strip
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double StripPitch = GPDSQUARE::FirstStageFace / m_NumberOfStrip;

   G4double temp = (POS(1) + GPDSQUARE::FirstStageFace / 2.) / StripPitch   ;
   G4int temp2 = temp ;
   G4int Y = temp2 + 1                    ;
   //Rare case where particle is close to edge of silicon plan
   if (Y == 129) Y = 128;

   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index =  aStep->GetTrack()->GetTrackID();
   EvtMap->set(DetNbr + index, Y);
   return TRUE;
}

void GPDScorerFirstStageBackStripSquare::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4int>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerFirstStageBackStripSquare::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerFirstStageBackStripSquare::clear()
{
   EvtMap->clear();
}

void GPDScorerFirstStageBackStripSquare::DrawAll()
{
}

void GPDScorerFirstStageBackStripSquare::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
   G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// FirstStage Front Strip position Scorer for Trapezoid geometry
GPDScorerFirstStageFrontStripTrapezoid::GPDScorerFirstStageFrontStripTrapezoid(G4String name, G4int depth, G4int NumberOfStrip)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_NumberOfStrip = NumberOfStrip;
}

GPDScorerFirstStageFrontStripTrapezoid::~GPDScorerFirstStageFrontStripTrapezoid()
{
}

G4bool GPDScorerFirstStageFrontStripTrapezoid::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, "GPDTrapezoid");

   // get front strip
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double StripPitch = GPDTRAP::FirstStageBaseLarge / m_NumberOfStrip;

   G4double temp = (POS(0) + GPDTRAP::FirstStageBaseLarge / 2.) / StripPitch   ;
   G4int X = int(temp) + 1;
   //Rare case where particle is close to edge of silicon plan
   if (X == 129) X = 128;
   
   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index =  aStep->GetTrack()->GetTrackID();
   EvtMap->set(DetNbr + index, X);
   return TRUE;
}

void GPDScorerFirstStageFrontStripTrapezoid::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4int>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerFirstStageFrontStripTrapezoid::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerFirstStageFrontStripTrapezoid::clear()
{
   EvtMap->clear();
}

void GPDScorerFirstStageFrontStripTrapezoid::DrawAll()
{
}

void GPDScorerFirstStageFrontStripTrapezoid::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
   G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// FirstStage Back Strip position Scorer for Trapezoid geometry
GPDScorerFirstStageBackStripTrapezoid::GPDScorerFirstStageBackStripTrapezoid(G4String name, G4int depth, G4int NumberOfStrip)
      : G4VPrimitiveScorer(name, depth), HCID(-1)
{
   m_NumberOfStrip = NumberOfStrip;
}

GPDScorerFirstStageBackStripTrapezoid::~GPDScorerFirstStageBackStripTrapezoid()
{
}

G4bool GPDScorerFirstStageBackStripTrapezoid::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   // get detector number
   int DetNbr = OBSOLETEGENERALSCORERS::PickUpDetectorNumber(aStep, "GPDTrapezoid");

   // get back strip
   G4ThreeVector POS  = aStep->GetPreStepPoint()->GetPosition();
   POS = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(POS);

   G4double StripPitch = GPDTRAP::FirstStageHeight / m_NumberOfStrip;

   G4double temp = (POS(1) + GPDTRAP::FirstStageHeight / 2.) / StripPitch   ;
   G4int temp2 = temp ;
   G4int Y = temp2 + 1                    ;
   //Rare case where particle is close to edge of silicon plan
   if (Y == 129) Y = 128;

   G4double edep = aStep->GetTotalEnergyDeposit();
   if (edep < TriggerThreshold) return FALSE;
   G4int  index =  aStep->GetTrack()->GetTrackID();
   EvtMap->set(DetNbr + index, Y);
   return TRUE;
}

void GPDScorerFirstStageBackStripTrapezoid::Initialize(G4HCofThisEvent* HCE)
{
   EvtMap = new G4THitsMap<G4int>(GetMultiFunctionalDetector()->GetName(), GetName());
   if (HCID < 0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void GPDScorerFirstStageBackStripTrapezoid::EndOfEvent(G4HCofThisEvent*)
{
}

void GPDScorerFirstStageBackStripTrapezoid::clear()
{
   EvtMap->clear();
}

void GPDScorerFirstStageBackStripTrapezoid::DrawAll()
{
}

void GPDScorerFirstStageBackStripTrapezoid::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
   G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
   G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}
