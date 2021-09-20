#ifndef Scone_h
#define Scone_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : March 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Scone simulation                             *
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
#include "G4AssemblyVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TSconeData.h"
#include "NPInputParser.h"

class Scone : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Scone() ;
    virtual ~Scone() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    void AddDetector(G4ThreeVector POS);

    G4LogicalVolume* Build2x2Assembly(int DetNumber);
    G4LogicalVolume* Build6x6Assembly(int DetNumber, double plastic_length);
    G4LogicalVolume* BuildSquareDetector();

    void Build2x2Block(G4LogicalVolume* world);
    void BuildRing1(G4LogicalVolume* world);
    void BuildRing2(G4LogicalVolume* world);
    G4AssemblyVolume* BuildFissionChamber();

  private:
    G4LogicalVolume* m_2x2Assembly;
    G4LogicalVolume* m_6x6Assembly;
    G4LogicalVolume* m_SquareDetector;
    G4AssemblyVolume* m_FissionChamberVolume;
    
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
    //   Initialize all Scorer used 
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_SconeScorer ;
    G4MultiFunctionalDetector* m_GdScorer ;
    G4MultiFunctionalDetector* m_FCScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TSconeData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_R; 
    vector<double>  m_Theta;
    vector<double>  m_Phi; 
    
    int m_BuildRing1;
    int m_BuildRing2;
    int m_BuildFissionChamber;
    int m_NumberOfInnerDetector;
    int m_NumberOfRing1Detector;
    int m_NumberOfRing2Detector;
    int m_Assembly;

  private: // Initalise material used in detector definition
    void InitializeMaterial();
    G4Material* m_MaterialVaccuum;


    //   Shape type
    vector<string> m_Shape ;
   
    // Visualisation Attribute
    G4VisAttributes* m_VisSquare;
    G4VisAttributes* m_Vis2x2;
    G4VisAttributes* m_Vis6x6R1;
    G4VisAttributes* m_Vis6x6R2;
    G4VisAttributes* m_VisGd;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
