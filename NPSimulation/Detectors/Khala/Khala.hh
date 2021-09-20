#ifndef Khala_h
#define Khala_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: M. Labiche  contact address: marc.labiche@stfc.ac.uk     *
 *                                                                           *
 * Creation Date  : December 2009                                            *
 * Last update    : December 2014                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe the Khala scintillator array                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>
#include <vector>

// G4 header defining G4 types
#include "globals.hh"

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4MultiFunctionalDetector.hh"
#include "NPInputParser.h"
#include "G4Box.hh"
// NPSimulation header
#include "NPSVDetector.hh"

// NPLib
#include "TKhalaData.h"
using namespace std;
using namespace CLHEP;

class Khala : public NPS::VDetector{
	////////////////////////////////////////////////////
	/////// Default Constructor and Destructor /////////
	////////////////////////////////////////////////////
	public:
		Khala() ;
		~Khala() ;

		////////////////////////////////////////////////////
		//////// Specific Function of this Class ///////////
		////////////////////////////////////////////////////
	public:
		// To add a Detector 
		// By corner position
		void AddDetector(G4ThreeVector Pos1,G4ThreeVector Pos2,G4ThreeVector Pos3,G4ThreeVector Pos4);
		// By Center position
		void AddDetector(G4ThreeVector Pos, double beta_u=0, double beta_v=0, double beta_w=0);

		// Return a logical volume of the detector
		G4LogicalVolume* ConstructDetector();

	private: // Guarranty that each volume is created only once
		G4LogicalVolume* m_LogicalDetector;

		////////////////////////////////////////////////////
		/////////  Inherite from NPS::VDetector class ///////////
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

		////////////////////////////////////////////////////
		///////////Event class to store Data////////////////
		////////////////////////////////////////////////////
	private:
		TKhalaData* m_Event ;

		////////////////////////////////////////////////////
		///////////////// Scorer Related ///////////////////
		////////////////////////////////////////////////////

	private:
		//   Initialize all Scorer
		void InitializeScorers() ;

		//   Scorer Associate to the Silicon
		G4MultiFunctionalDetector* m_LaBr3Scorer ;

	private:
		////////////////////////////////////////////////////
		///////////////Private intern Data//////////////////
		////////////////////////////////////////////////////
	private:
		// Detector position
		vector<G4ThreeVector> m_Pos;
		vector<G4RotationMatrix*> m_Rot;

	private:/// Visualisation Attribute:
		G4VisAttributes* m_LaBr3VisAtt;
		G4VisAttributes* m_DetectorCasingVisAtt  ;
		G4VisAttributes* m_PMTShieldVisAtt;
		G4VisAttributes* m_PMTVisAtt;
		G4VisAttributes* m_ShieldCanVisAtt;
		G4VisAttributes* m_HolderCanVisAtt;
		G4VisAttributes* m_HolderVisAtt;
		G4VisAttributes* m_Holder1VisAtt;
		G4VisAttributes* m_Holder2VisAtt;
		G4VisAttributes* m_Holder3VisAtt;
		G4VisAttributes* m_Holder4VisAtt;
		G4VisAttributes* m_Holder5VisAtt;
		G4VisAttributes* m_Holder6VisAtt;
		G4VisAttributes* m_Holder7VisAtt;



		G4double    m_ChamberHmin;
		G4double    m_ChamberHmax;
		G4double    m_ChamberWmin;
		G4double    m_ChamberWmax;
		G4double    m_ChamberDmin; 
		G4double    m_ChamberDmax;

	public:
		static NPS::VDetector* Construct();
};

namespace KHALA{
	// Resolution
	const G4double EnergyResolution = 0.0099;	// = 3.5% at .662MeV of Resolution   //   Unit is MeV/2.35
	const G4double EnergyThreshold = .01*MeV;   

	// Geometry for the mother volume 
	const G4double FaceFront = 4.5*cm;
	const G4double Length = 22.0*cm; 

	// LaBr3
	const G4double LaBr3Face = 3.81*cm;
	const G4double LaBr3Thickness  = 3.81*cm;

	// Al Can
	const G4double CanOuterDiameter = 4.3*cm;
	const G4double CanInnerDiameter = 4.2*cm;
	const G4double CanLength = 4.0*cm; 

	// Al front Window
	const G4double WinOuterDiameter = 4.2*cm;
	const G4double WinInnerDiameter = 0*cm;
	const G4double WinLength  = 0.05*cm;

	// PMT 
	const G4double PMTThickness = 11.0*cm;
	const G4double PMTIn = 4.0*cm;
	const G4double PMTOut = 4.5*cm;

	// PMT back Window
	const G4double PMTCoverOuterDiameter = 4.0*cm;
	const G4double PMTCoverInnerDiameter = 0*cm;
	const G4double PMTCoverLength = 0.1*cm;

	// PM Base
	const G4double PMThickness = 7.0*cm;
	const G4double PMIn = 4.5*cm;
	const G4double PMOut = 4.8*cm;

	// PM Base back Window
	const G4double PMCoverOuterDiameter = 4.5*cm;
	const G4double PMCoverInnerDiameter = 0*cm;
	const G4double PMCoverLength = 0.1*cm;

	// Position 
	const G4double LaBr3_PosZ  = -Length*0.5 + 0.5*LaBr3Thickness + 0.1*cm;
	const G4double LaBr3Can_PosZ    = -Length*0.5 + 0.5*CanLength;
	const G4double LaBr3Win_PosZ    = -Length*0.5 + 0.5*WinLength;
	const G4double PMT_PosZ    = -Length*0.5 + 0.5*PMTThickness + CanLength;
	const G4double PMTCover_PosZ = -Length*0.5 + PMTThickness + CanLength + 0.5*PMTCoverLength;
	const G4double PM_PosZ    = -Length*0.5 + (Length-PMThickness) + 0.5*PMThickness;
	const G4double PMCover_PosZ = -Length*0.5 + PMTThickness + PMThickness + CanLength + 0.5*PMCoverLength;

}

#endif
