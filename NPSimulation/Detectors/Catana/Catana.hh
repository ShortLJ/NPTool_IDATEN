#ifndef Catana_h
#define Catana_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Catana simulation                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * Geometry of crystal based on official catana simulation from Samurai      *
 * Collaboration package 5.2                                                 *
 * http://be.nucl.ap.titech.ac.jp/~nebula/simulator.php                      *
 * Thanks to Togano-san                                                      *
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
#include "TCatanaData.h"
#include "NPInputParser.h"

class Catana : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Catana() ;
    virtual ~Catana() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(double X, double Y, double Z, double Theta, double Phi, int ID,int Type,double Rshift=0);
    void ReadCSV(string path,double Rshift);

    G4LogicalVolume* BuildDetector(int Type);

  private:
    G4LogicalVolume* m_DetectorType1;
    G4LogicalVolume* m_DetectorType2;
    G4LogicalVolume* m_DetectorType3;
    G4LogicalVolume* m_DetectorType4;
    G4LogicalVolume* m_DetectorType5;
    

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
    G4MultiFunctionalDetector* m_CatanaScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TCatanaData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_X; 
    vector<double>  m_Y; 
    vector<double>  m_Z; 
    vector<double>  m_Theta; 
    vector<double>  m_Phi; 
    vector<int>     m_ID;
    vector<int>     m_Type;
    G4ThreeVector   m_Ref;
    // this parameter is here because some csv file have very small overlap
    // due to difference between mechanical design and reality of the detector
    // a shift is apply to the position of the crystal to slightly icrease the radius
    // and avoid shift. Typical value shoulde be < 100um
    vector<double>  m_Rshift;// additional shift to apply to csv file
    // relative shift of crystal w/r to the housing
    map<int,double>  m_Zoffset;

    // Visualisation Attribute
    G4VisAttributes* m_VisCrystal1;
    G4VisAttributes* m_VisCrystal2;
    G4VisAttributes* m_VisCrystal3;
    G4VisAttributes* m_VisCrystal4;
    G4VisAttributes* m_VisCrystal5;

    G4VisAttributes* m_VisHousing;
    G4VisAttributes* m_VisReflector;



  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
