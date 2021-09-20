#ifndef Minos_h
#define Minos_h 1
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin                                        *
 * Maintainer : Cyril Lenain       contact address: lenain@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date  : October 2018                                             *
 * Last update    : April 2021                                               *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Minos simulation                                    *
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
#include "G4VFastSimulationModel.hh"
#include "G4UserLimits.hh"
#include "G4FastSimulationManager.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TMinosData.h"
#include "NPInputParser.h"
#include "Decay.hh"
#include "BeamReaction.hh"
#include "NPXmlParser.h"

//Root
#include "TF1.h"

class Minos : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Minos() ;
    virtual ~Minos() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
  // With TargetLenght
    void AddDetector(G4ThreeVector POS, double TargetLength, G4String MaterialOfTarget,G4String MaterialOfCell, int TPCOnly);
    /* void AddDetector(G4ThreeVector POS, double TargetLength, int TPCOnly); */
 
  private:
  //For material definition
    void DefineMaterials();
  
  public:
    G4LogicalVolume* BuildSquareDetector();
    G4LogicalVolume* BuildCylindricalDetector();
    G4LogicalVolume* BuildTarget();
    G4LogicalVolume* BuildChamber();
    G4LogicalVolume* BuildInnerRohacell();
    G4LogicalVolume* BuildOuterRohacell();
    G4LogicalVolume* BuildKapton();
    G4LogicalVolume* BuildOuterKapton();
    G4LogicalVolume* BuildTPC();
    G4LogicalVolume* BuildWindow0();
    G4LogicalVolume* BuildWindow1();
    G4LogicalVolume* BuildWindow2();
 
 public:
     
     /* G4double    GetTargetLength()      {return TargetLength*2.;}; */
     /* G4Material* GetTargetMaterial()    {return TargetMaterial;}; */
     /* G4double    GetTargetRadius()      {return TargetRadius;}; */

  private:
     G4Material*        TargetMaterial;
     G4Material*        WindowMaterial;
     G4Material*        ChamberMaterial;
     G4Material*        RohacellMaterial;
     G4Material*        KaptonMaterial;
     G4Material*        TPCMaterial;
     G4Material*        defaultMaterial;
     G4Material*        Mylar;


    G4LogicalVolume* m_SquareDetector;
    G4LogicalVolume* m_CylindricalDetector;

     G4Tubs*             solidTarget;   
     G4LogicalVolume*   logicTarget;   
  // G4VPhysicalVolume* physiTarget;   
     
     G4Tubs*             solidChamber;  
     G4LogicalVolume*   logicChamber;  
  // G4VPhysicalVolume* physiChamber;  
         
     G4Tubs*             solidTPC; 
     G4LogicalVolume*   logicTPC; 
  // G4VPhysicalVolume* physiTPC; 
 
  ////////////////
     G4Tubs*             solidWindow0; 
     G4LogicalVolume*   logicWindow0; 
  // G4VPhysicalVolume* physiWindow0; 
     
     G4Tubs*             solidRohacell;   
     G4LogicalVolume*   logicRohacell;   
  // G4VPhysicalVolume* physiInnerRohacell;   
      
     /* G4Tubs*             solidOuterRohacell; */   
     /* G4LogicalVolume*   logicOuterRohacell; */   
  /* // G4VPhysicalVolume* physiOuterRohacell; */   
     
     G4Tubs*             solidKapton;   
     G4LogicalVolume*   logicKapton;   
  // G4VPhysicalVolume* physiKapton;   
    
     G4double TargetLength;
     G4int    TPCOnly;

     G4double start,end,time, DriftTime;
     TH1F* Raw_Signal ;      
     TH1F* Elec_Signal;
     TF1* fa1;   
     vector<int> Q, T;
     
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
    void SimulateGainAndDigitizer(vector<double>* rawT, vector<double>* rawQ,vector<int>& Q,vector<int>& T);

  private: // parameter for the digitization  
    double m_TimeBin;
    double m_ShapingTime;
    double m_Baseline;
    unsigned int m_Sampling;
    double m_ZOffset;


    //   Associated Scorer
    G4MultiFunctionalDetector* m_MinosPadScorer ;
  
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TMinosData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: 
    
  // Geometry
  // Detector Coordinate 
  vector<G4ThreeVector>  m_POS; 
  vector<double>         m_TargetLength;
  vector<G4String>       m_TargetMaterial;
  vector<G4String>       m_CellMaterial;
  vector<int>            m_TPCOnly;

  // Visualisation Attribute
  G4VisAttributes* m_VisTarget;
  G4VisAttributes* m_VissimpleBox;
  G4VisAttributes* m_VisTPC;
  G4VisAttributes* m_VisRohacell;
  G4VisAttributes* m_VisKapton;
  G4VisAttributes* m_VisTargetCell;
  G4VisAttributes* m_VisOuterKapton;

  private:
  // XML ID to G4ID
  // map of ID, coordinate
  map<unsigned int,pair<double,double> > m_XY;
  // map of G4ID to XML ID
  map<unsigned int,unsigned int> m_ID;
  void ReadXML(NPL::XmlParser& xml);
  unsigned int FindPadID(unsigned int G4ID, double X,double Y);

  private:
    
  // Region were reaction can occure:
  G4Region* m_ReactionRegion;
  vector<G4VFastSimulationModel*> m_ReactionModel;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};

#endif
