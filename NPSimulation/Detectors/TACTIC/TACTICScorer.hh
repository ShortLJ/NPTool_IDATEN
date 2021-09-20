#ifndef TACTICScorer_h
#define TACTICScorer_h 1

#include "G4VPrimitiveScorer.hh"
#include "NPSHitsMap.hh"
#include <map>

using namespace std;
using namespace CLHEP;

namespace TACTICScorer {
  
  class Gas_Scorer : public G4VPrimitiveScorer{
    
  public:
    Gas_Scorer(G4String name,
	       G4int Level,
	       G4double ScorerLength,
	       G4int NumberOfSegments,
	       G4int depth=0,
	       G4double p0 = 0,
	       G4double p1 = 0,
	       G4double p2 = 0,
	       G4double p3 = 0,
	       string Shape = "default"
	       );
    
    ~Gas_Scorer();
    
  protected: // with description                                                                                                                            
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    
  public:
    void Initialize(G4HCofThisEvent*);
    void EndOfEvent(G4HCofThisEvent*);
    void clear();
    void DrawAll();
    void PrintAll();
    
    
  private: // Geometry of the detector                                                                                                                      
    G4double m_ScorerLength;
    G4int    m_NumberOfSegments;
    G4double m_SegmentLength;
    // Level at which to find the copy number linked to the detector number                                                                                 
    G4int    m_Level;
    G4double m_p0;
    G4double m_p1;
    G4double m_p2;
    G4double m_p3;
    string m_Shape; 
    ofstream file;

    //    G4double excess;
    
  private: // inherited from G4VPrimitiveScorer                                                                                                             
    G4int HCID;
    NPS::HitsMap<G4double*>* EvtMap;
    
    
  private: // Needed for intermediate calculation (avoid multiple instantiation in Processing Hit)                                                          
    G4ThreeVector m_Position  ;
    G4int m_DetectorNumber    ;
    G4int m_SegmentNumber ;
    G4long m_Index             ;
    
  };
  
}

#endif
