#ifndef FissionChamber_h
#define FissionChamber_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : September 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  FissionChamber simulation                             *
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
#include "G4AssemblyVolume.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TFissionChamberData.h"
#include "NPInputParser.h"

class FissionChamber : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    FissionChamber() ;
    virtual ~FissionChamber() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(G4ThreeVector POS);
    // Spherical
    void AddDetector(double R,double Theta,double Phi);  


    G4AssemblyVolume* BuildDetector();
    G4AssemblyVolume* BuildFissionChamber();
    void BuildAnode(double PosZ);
    void BuildCathode(double PosZ);

  private:
    G4LogicalVolume* m_Detector;
    G4AssemblyVolume* m_AssemblyVolume;

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
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_FissionChamberScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TFissionChamberData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_R; 
    vector<double>  m_Theta;
    vector<double>  m_Phi; 

    string m_GasMaterial;
    double m_Pressure;

    // Visualisation Attribute
    G4VisAttributes* m_VisFCWall;
    G4VisAttributes* m_VisAl;
    G4VisAttributes* m_VisCu;
    G4VisAttributes* m_VisGas;
    G4VisAttributes* m_VisTi;
    G4VisAttributes* m_VisRogers4003C;

    // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
