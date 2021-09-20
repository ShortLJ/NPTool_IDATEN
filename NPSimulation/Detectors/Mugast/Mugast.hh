#ifndef Mugast_h
#define Mugast_h 1
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2019                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Mugast simulation                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TMugastData.h"
#include "NPInputParser.h"

class Mugast : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Mugast() ;
    virtual ~Mugast() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(int DetectorNumber,
      string Shape,
      G4ThreeVector PX1_Y1 ,
      G4ThreeVector PX1_Y128=G4ThreeVector() ,
      G4ThreeVector PX128_Y1=G4ThreeVector(),
      G4ThreeVector PX128_Y128=G4ThreeVector());

    G4LogicalVolume* BuildSquareDetector();
    G4LogicalVolume* BuildTrapezoidDetector();
    G4LogicalVolume* BuildAnnularDetector();
  
  private:
    G4LogicalVolume* m_SquareDetector;
    G4LogicalVolume* m_TrapezoidDetector;
    G4LogicalVolume* m_AnnularDetector;
    
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetextorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;

    // Construct detector and inialise sensitive part.
    // Called After DetecorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Detector branch to the EventTree.
    // Called After DetecorConstruction::AddDetector Method
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event) ;

  public:   // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_SquareScorer;
    G4MultiFunctionalDetector* m_TrapezoidScorer;
    G4MultiFunctionalDetector* m_AnnularScorer;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TMugastData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    // Used for "By Point Definition"
    vector<G4ThreeVector>   m_X1_Y1     ; // Top Left Corner Position Vector
    vector<G4ThreeVector>   m_X1_Y128   ; // Bottom Left Corner Position Vector
    vector<G4ThreeVector>   m_X128_Y1   ; // Bottom Right Corner Position Vector
    vector<G4ThreeVector>   m_X128_Y128 ; // Center Corner Position Vector
    
    //   Shape type
    vector<string> m_Shape ;
    // DetectorNumber
    vector<int>    m_DetectorNumber;
    // Visualisation Attribute
    G4VisAttributes* m_VisSquare;
    G4VisAttributes* m_VisCylinder;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
