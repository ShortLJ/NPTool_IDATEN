#ifndef eAGanil_h
#define eAGanil_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  eAGanil simulation                                  *
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
#include "TeAGanilData.h"
#include "NPInputParser.h"

class eAGanil : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    eAGanil() ;
    virtual ~eAGanil() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:

    void AddDetector(double  R, double  Theta, double  Phi, double EntranceWidth,double EntranceHeigh,double MR);
    void SetTrap(double Length,double  InnerRadius, double  OuterRadius, double BladesThickness, int NumberOfBlades,double  Phi, double WindowsThickness);

    G4LogicalVolume* BuildDetector(unsigned int i);
    G4LogicalVolume* BuildTrap();

  private:
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
    G4MultiFunctionalDetector* m_eAGanilScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TeAGanilData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_SpecR; 
    vector<double>  m_SpecTheta;
    vector<double>  m_SpecPhi; 
    vector<double>  m_SpecEntranceWidth;
    vector<double>  m_SpecEntranceHeigh; 
    vector<double>  m_SpecMomentumResolution; 

    // Trap
    double m_Length;
    double m_InnerRadius; 
    double m_OuterRadius;
    double m_BladesThickness;
    double m_NumberOfBlades;
    double m_Phi;       
    double m_WindowsThickness;

    // Visualisation Attribute
    G4VisAttributes* m_VisDetector;
    G4VisAttributes* m_VisTrap;
    // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
