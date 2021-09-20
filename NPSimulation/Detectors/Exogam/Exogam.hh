#ifndef Exogam_h
#define Exogam_h 1
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Goigoux  contact address: thomas.goigoux@cea.fr                        *
 *                                                                           *
 * Creation Date  : july 2019                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Exogam simulation                             *
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
#include "TExogamData.h"
#include "NPInputParser.h"

class Exogam : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Exogam() ;
    virtual ~Exogam() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Spherical
    void AddDetector(double X,double Y, double Z, double ThetaX, double ThetaY, double ThetaZ);  

    G4int InitializeMaterials();
    void BuildClover(int i_clo, G4LogicalVolume* world);
    void BuildSideCatcher();
    void BuildBackCatcher();
    void BuildSideShield();
    void BuildCollimator();
  
 
    
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
    G4MultiFunctionalDetector* m_ExogamScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TExogamData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate (of the front of the Germanium)
    vector<double>  m_X; 
    vector<double>  m_Y; 
    vector<double>  m_Z; 
    // Detector orientation
    vector<double>  m_ThetaX; //rotation angles to X, Y, Z axis
    vector<double>  m_ThetaY;
    vector<double>  m_ThetaZ;
    
 private:
    //materials
    G4Material* m_Vacuum;
    G4Material* m_Aluminum;
    G4Material* m_Copper;
    G4Material* m_Germanium;
    G4Material* m_BGO;
    G4Material* m_CsI;

    // Some rotation matrices
    G4RotationMatrix rm90;
    G4RotationMatrix rm90m;
    G4RotationMatrix rmCut2;
    G4RotationMatrix rm180;
    G4RotationMatrix rm270;

    G4double HalfLengthCan;
    G4double TaperLengthCan;
    G4double distCollimatorToBGOSShield;
    G4double Offset;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
