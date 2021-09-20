#ifndef TACTIC_h
#define TACTIC_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Warren Lynch  contact address: Warren.Lynch@york.ac.uk                        *
 *                                                                           *
 * Creation Date  : June 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  TACTIC simulation                             *
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
#include "G4UserLimits.hh"
#include "G4VFastSimulationModel.hh"
#include "G4FastSimulationManager.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TTACTICData.h"
#include "NPInputParser.h"
#include "Decay.hh"
#include "BeamReaction.hh"

extern double excess;

class TACTIC : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
public:
  TACTIC() ;
  virtual ~TACTIC() ;
  
  ////////////////////////////////////////////////////
  /////// Specific Function of this Class ///////////
  ////////////////////////////////////////////////////
public:
  // Cartesian
  void AddDetector(G4ThreeVector POS, string Shape);
  // Spherical
  void AddDetector(double R,double Theta,double Phi,string Shape);  
  
  G4LogicalVolume* BuildCylindricalDetector();
  
private:
  G4LogicalVolume* m_CylindricalDetector;
  G4LogicalVolume* gas_volume_log;
  G4LogicalVolume* window_log;
  //G4LogicalVolume* window_log_2;
  //G4LogicalVolume* vacuum_log;
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
  G4MultiFunctionalDetector* m_Scorer ;
  ////////////////////////////////////////////////////
  ///////////Event class to store Data////////////////
  ////////////////////////////////////////////////////
private:
  TTACTICData* m_Event ;
  
  ////////////////////////////////////////////////////
  ///////////////Private intern Data//////////////////
  ////////////////////////////////////////////////////
private: // Geometry
  // Detector Coordinate 
  vector<double>  m_R; 
  vector<double>  m_Theta;
  vector<double>  m_Phi; 
  
  vector<string> m_GasMaterial;
  vector<int> m_GasFraction;
  
  double m_Pressure;
  double m_Temperature;
  //   Shape type
  vector<string> m_Shape ;
  string m_Active;
  double m_p0, m_p1, m_p2, m_p3;
  string Shape;
  
  //int NumberOfStrips;
  
  // Visualisation Attribute
  G4VisAttributes* m_VisChamber;
  G4VisAttributes* m_VisWindows;
  G4VisAttributes* m_VisGas;
  G4VisAttributes* m_VisVacuum;
  
private:
  // Region were reaction can occure:
  G4Region* m_ReactionRegion;
  vector<G4VFastSimulationModel*> m_ReactionModel;
  
  // Needed for dynamic loading of the library
public:
  static NPS::VDetector* Construct();

};

#endif
