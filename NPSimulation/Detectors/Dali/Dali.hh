#ifndef Dali_h
#define Dali_h 1
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: elidiano.tronchin@studenti.unipd.it                        *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Dali simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
#include <cmath>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TDaliData.h"
#include "NPInputParser.h"

class Dali : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Dali() ;
    virtual ~Dali() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(G4ThreeVector POS);
    // Spherical
    void AddDetector(double R,double Theta,double Phi);  
    //Cylindrical
  void AddDetector(double R,double Alpha,double Zeta, int Ring);

    G4LogicalVolume* BuildSquareDetector();
    G4LogicalVolume* BuildCylindricalDetector();
  
  private:
  
    G4LogicalVolume* m_SquareDetector;
    G4LogicalVolume* Log_Dali_box;
    G4LogicalVolume* Log_Al_Cryst_can;
    G4LogicalVolume* Log_MgO_Cryst_can;
    G4LogicalVolume*  Log_Crystal;
    G4LogicalVolume* lAlPMT;
    G4LogicalVolume* Log_Dali_1Volume;
    G4LogicalVolume* Log_Dali_3Volume;
    G4LogicalVolume* m_CylindricalDetector;
    G4LogicalVolume* lMuPMT;
    G4LogicalVolume* lTopPlatePMT;
    G4LogicalVolume* lGlassPMT;
    G4LogicalVolume* AriaExtrude;

    G4Material* NaI_Tl;

    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetextorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;

    // Definition of materials
    // Called in ConstructDetector()
    void DefinitionMaterials();
  
    // Construct detector and inialise sensitive part.
    // Called After DetecorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Detector branch to the EventTree.
    // Called After DetecorConstruction::AddDetector Method
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* ) ;

  public:   // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_DaliScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TDaliData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_Zeta;
    vector<double>  m_R; 
    vector<double>  m_Alpha; 
    vector<int>  m_Ring; 
    
    // Visualisation Attribute
    G4VisAttributes* m_VisSquare;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
