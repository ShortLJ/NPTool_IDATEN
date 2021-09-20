#ifndef Strasse_h
#define Strasse_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta  contact address: matta@lpccaen.in2p3.fr        *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Strasse simulation                                  *
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
#include "TStrasseData.h"
#include "NPInputParser.h"

class Strasse : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Strasse() ;
    virtual ~Strasse() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cylindrical coordinate
    void AddInnerDetector(double R,double Z,double Phi, double Shift, G4ThreeVector Ref);  
    void AddOuterDetector(double R,double Z,double Phi, double Shift, G4ThreeVector Ref);  
    void AddChamber(double Z);

    G4LogicalVolume* BuildInnerDetector();
    G4LogicalVolume* BuildOuterDetector();
    G4LogicalVolume* BuildElectronic();
    G4LogicalVolume* BuildChamber();
    G4LogicalVolume* BuildChamberFromCAD(string path);
    G4LogicalVolume* BuildStars(string path);
    G4LogicalVolume* BuildBlades(string path);
    G4LogicalVolume* BuildBase(string path);

  private:
    G4LogicalVolume* m_InnerDetector;
    G4LogicalVolume* m_OuterDetector;
    G4LogicalVolume* m_Electronic;
    G4LogicalVolume* m_Stars;
    G4LogicalVolume* m_Chamber;
    G4LogicalVolume* m_Blades;
    G4LogicalVolume* m_Base;

    string ChamberPath;
    string BasePath;
    string StarsPath;
    string BladesPath;
    bool found_chamber;
    bool found_blades;
    bool found_stars;
    bool found_base;

  private:
    //    Initialize material used in detector definition
    void InitializeMaterial();


    //   List of material
    G4Material* m_MaterialSilicon ;
    G4Material* m_MaterialAl      ;
    G4Material* m_MaterialVacuum  ;
    G4Material* m_MaterialPCB     ;
    G4Material* m_MaterialCu     ;

    // calculated dimension
    double m_Active_InnerWafer_Width;
    double m_Active_InnerWafer_Length; 
    double m_Active_OuterWafer_Width;
    double m_Active_OuterWafer_Length; 


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
    G4MultiFunctionalDetector* m_InnerScorer1 ;
    G4MultiFunctionalDetector* m_OuterScorer1 ;
    G4MultiFunctionalDetector* m_InnerScorer2 ;
    G4MultiFunctionalDetector* m_OuterScorer2 ;

    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TStrasseData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_Inner_R; 
    vector<double>  m_Inner_Z;
    vector<double>  m_Inner_Phi; 
    vector<double>  m_Inner_Shift; 
    vector<G4ThreeVector> m_Inner_Ref;

    vector<double>  m_Outer_R; 
    vector<double>  m_Outer_Z;
    vector<double>  m_Outer_Phi; 
    vector<double>  m_Outer_Shift; 
    vector<G4ThreeVector> m_Outer_Ref;

    vector<double>  m_Chamber_Z;


    // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();


  private: // Visualisation
    G4VisAttributes* SiliconVisAtt  ;
    G4VisAttributes* PCBVisAtt;
    G4VisAttributes* PADVisAtt  ;
    G4VisAttributes* StarsVisAtt ;
    G4VisAttributes* ChamberVisAtt ;
    G4VisAttributes* GuardRingVisAtt ;
    G4VisAttributes* BladeVisAtt ;


};
#endif
