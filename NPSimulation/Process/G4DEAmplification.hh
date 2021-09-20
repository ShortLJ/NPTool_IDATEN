//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4DEAmplification.hh 69576 2013-05-08 13:48:13Z gcosmo $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Amplification Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4DEAmplification.hh
// Description: Discrete Process -- Bulk absorption of Drift Electron
// Version:     1.0
// Created:     2016-12-19 
// Author:      Adrien Matta
// Updated:     
//             
//             
//              
//              
// mail:        matta@lpccaen.in2p3.fr 
//
////////////////////////////////////////////////////////////////////////

#ifndef G4DEAmplification_h
#define G4DEAmplification_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4DriftElectron.hh"

// Class Description:
// Discrete Process -- Bulk amplification of Drift Electron.
// Class inherits publicly from G4VDiscreteProcess
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4DEAmplification : public G4VRestDiscreteProcess 
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4DEAmplification(const G4String& processName = "DEAmplification",
                                G4ProcessType type = fTransportation);
	~G4DEAmplification();

private:

        G4DEAmplification(const G4DEAmplification &right);

        //////////////
        // Operators
        //////////////

        G4DEAmplification& operator=(const G4DEAmplification &right);

public:

	////////////
	// Methods
  ////////////

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
  // Returns true -> 'is applicable' only for an drift electron.

	G4double GetMeanFreePath(const G4Track& aTrack,
				 G4double ,
				 G4ForceCondition* );
  // Returns true and strongly enforced

	G4double GetMeanLifeTime(const G4Track& aTrack,
				 G4ForceCondition* );
  // Returns true and strongly enforced

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
 				        const G4Step&  aStep);
   // This is the method implementing bulk amplification of drift
   // electron.

};

////////////////////
// Inline methods
////////////////////

inline
G4bool G4DEAmplification::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return ( &aParticleType == G4DriftElectron::DriftElectron() );
}

#endif /* G4DEAmplification_h */
