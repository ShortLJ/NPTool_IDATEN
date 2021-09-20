#ifndef GPDScorer_h
#define GPDScorer_h 1
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

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

namespace GPDSCORERS
{
   // This Threshold is used in all scorers
   // Any energy deposit under this threshold will not create an entry
   const double TriggerThreshold = 0.1 * keV;

class GPDScorerFirstStageEnergy : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageEnergy(G4String name, G4String volumeName, G4int depth = 0);
   virtual ~GPDScorerFirstStageEnergy();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4String m_VolumeName;
   G4int HCID;
   G4THitsMap<G4double>* EvtMap;
};



class GPDScorerSecondStageEnergy : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerSecondStageEnergy(G4String name, G4String volumeName, G4int depth = 0);
   virtual ~GPDScorerSecondStageEnergy();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4String m_VolumeName;
   G4int HCID;
   G4THitsMap<G4double>* EvtMap;
};



class GPDScorerThirdStageEnergy : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerThirdStageEnergy(G4String name, G4String volumeName, G4int depth = 0);
   virtual ~GPDScorerThirdStageEnergy();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4String m_VolumeName;
   G4int HCID;
   G4THitsMap<G4double>* EvtMap;
};




class GPDScorerFirstStageFrontStripDummyShape : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageFrontStripDummyShape(G4String name, G4int depth = 0, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageFrontStripDummyShape();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};



class GPDScorerFirstStageBackStripDummyShape : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageBackStripDummyShape(G4String name, G4int depth = 0, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageBackStripDummyShape();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};



class GPDScorerFirstStageFrontStripSquare : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageFrontStripSquare(G4String name, G4int depth = 0, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageFrontStripSquare();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};



class GPDScorerFirstStageBackStripSquare : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageBackStripSquare(G4String name, G4int depth = 0, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageBackStripSquare();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};



class GPDScorerFirstStageFrontStripTrapezoid : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageFrontStripTrapezoid(G4String name, G4int depth = 0, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageFrontStripTrapezoid();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};



class GPDScorerFirstStageBackStripTrapezoid : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageBackStripTrapezoid(G4String name, G4int depth = 0, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageBackStripTrapezoid();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};



class GPDScorerFirstStageFrontStripAnnular : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageFrontStripAnnular(G4String name, G4int depth = 0, G4double StripPlaneSize = 98, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageFrontStripAnnular();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4double  m_StripPlaneSize;
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};



class GPDScorerFirstStageBackStripAnnular : public G4VPrimitiveScorer
{
public: // with description
   GPDScorerFirstStageBackStripAnnular(G4String name, G4int depth = 0, G4double StripPlaneSize = 98, G4int NumberOfStrip = 128);
   virtual ~GPDScorerFirstStageBackStripAnnular();

protected: // with description
   virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

public:
   virtual void Initialize(G4HCofThisEvent*);
   virtual void EndOfEvent(G4HCofThisEvent*);
   virtual void clear();
   virtual void DrawAll();
   virtual void PrintAll();

private:
   G4double  m_StripPlaneSize;
   G4int     m_NumberOfStrip ;
   G4int HCID;
   G4THitsMap<G4int>* EvtMap;
};

}
#endif
