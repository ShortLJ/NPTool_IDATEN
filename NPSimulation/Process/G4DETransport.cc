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
// $Id: G4DETransport.hh 69576 2016-12-19 13:48:13Z gcosmo $
//
////////////////////////////////////////////////////////////////////////
// Optical DriftElectron Transport Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4DETransport.cc
// Description: Discrete Process -- Bulk transport of Drift Electron
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

#include "G4ios.hh"
#include "G4DEProcessSubType.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DETransport.hh"
#include "G4FieldManager.hh"
#include "G4ElectroMagneticField.hh"

/////////////////////////
// Class Implementation
/////////////////////////

//////////////
// Operators
//////////////

// G4DETransport::operator=(const G4DETransport &right)
// {
// }

/////////////////
// Constructors
/////////////////

G4DETransport::G4DETransport(const G4String& processName, G4ProcessType type)
: G4VDiscreteProcess(processName, type)
{
    if (verboseLevel>0) {
        G4cout << GetProcessName() << " is created " << G4endl;
    }
    
    SetProcessSubType(fDETransport);
    G4TransportationManager* TransportationManager = G4TransportationManager::GetTransportationManager();
    m_SafetyHelper = TransportationManager->GetSafetyHelper();
}

// G4DETransport::G4DETransport(const G4OpAbsorpton &right)
// {
// }

////////////////
// Destructors
////////////////

G4DETransport::~G4DETransport(){}

////////////
// Methods

// PostStepDoIt
// -------------
//
G4VParticleChange*

G4DETransport::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep){
  aParticleChange.Initialize(aTrack);
  
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  //G4StepPoint* pPostStepPoint  = aStep.GetPostStepPoint();
  
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double      t0 = pPreStepPoint->GetGlobalTime();
  
  // The time scale is imposed by the distance travelled
  G4double step_length = aStep.GetDeltaPosition().mag();
  G4double delta_time = aStep.GetDeltaTime();

  // allow internal relocation of the track by the kernel taking a 0 length intermediate step
  // suppress also parasite infinetismal step that slow down the tracking
  if(step_length<100*micrometer){
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
  // Get the material table
  const G4Material* aMaterial = aTrack.GetMaterial();
  
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
  aMaterial->GetMaterialPropertiesTable();

  if (!aMaterialPropertiesTable)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  if(!aMaterialPropertiesTable->ConstPropertyExists("DE_PAIRENERGY") ||
     !aMaterialPropertiesTable->ConstPropertyExists("DE_TRANSVERSALSPREAD") ||
     !aMaterialPropertiesTable->ConstPropertyExists("DE_LONGITUDINALSPREAD") ||
     !aMaterialPropertiesTable->ConstPropertyExists("DE_DRIFTSPEED") )
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
 
  // Electron follow the field direction
  // The field direction is taken from the field manager
  static G4double* fieldArr = new G4double[6];
  G4double  Point[4]={x0.x(),x0.y(),x0.z(),t0};
  G4FieldManager*  fMng = pPreStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetFieldManager();
    
  G4ElectroMagneticField* field = (G4ElectroMagneticField*)fMng->GetDetectorField();
  field->GetFieldValue(Point,fieldArr) ;
  
  // Electron move opposite to the field direction, hance the minus sign
  G4ThreeVector driftDir(-fieldArr[3],-fieldArr[4],-fieldArr[5]);
  // Normalised the drift direction
  driftDir = driftDir.unit();
  
  G4double v_drift = aMaterialPropertiesTable->GetConstProperty("DE_DRIFTSPEED");
  // Speed sigma in transversal and longitudinal direction based on the 
  // step length, drift speed and long/trans spread parameters
  G4double sigmaTrans  = sqrt(2*step_length*aMaterialPropertiesTable->GetConstProperty("DE_TRANSVERSALSPREAD")/v_drift)/delta_time;
  G4double sigmaLong   = sqrt(2*step_length*aMaterialPropertiesTable->GetConstProperty("DE_LONGITUDINALSPREAD")/v_drift)/delta_time;
   
  // Building the modified drift vector
  // i.e. the speed vector of the particle at the end of the step
  G4double d_trans = G4RandGauss::shoot(0,sigmaTrans);
  G4double d_long  = G4RandGauss::shoot(v_drift,sigmaLong);
  
  // The transverse component is perpendicular to the drift direction,
  // to build an arbritary trans vector, we have to make a cross product with
  // an arbritary, non parallel vector, and  driftDir
  G4ThreeVector arbritary(1,0,0);
  G4ThreeVector trans;
  if(driftDir.dot(arbritary)!=1) // not orthogonal (both are unitary)
    trans=d_trans*(driftDir.cross(arbritary)).unit();
  else{
    arbritary=G4ThreeVector(0,1,0);
    trans=d_trans*(driftDir.cross(arbritary)).unit();
    }
 
  // Rotate randomly around driftDir for the Phi component
  trans.rotate(driftDir,twopi*G4UniformRand()); 

  // new position is drift vector
  G4ThreeVector drift_vector=d_long*driftDir.unit()+trans.unit()*d_trans;
  double speed = drift_vector.mag(); 
//  std::cout << " " << step_length << " " << speed << " " << d_trans << " " << d_long-v_drift << std::endl;

  G4double Energy= (0.5*((electron_mass_c2/joule)/(c_squared/((m*m)/(s*s))))*(speed/(m/s))*(speed/(m/s)))*joule;
  aParticleChange.ProposeMomentumDirection(drift_vector.unit());
  aParticleChange.ProposeEnergy(Energy);
  
  return &aParticleChange;
}


// GetMeanFreePath
// ---------------
//
G4double G4DETransport::GetMeanFreePath(const G4Track& track,
                                        G4double ,
                                        G4ForceCondition* )
{
  // Typical distance after which the electron position should be reevaluated
  // to take into account the diffusivity of the charge cloud inside the 
  // medium
  double d = m_SafetyHelper->ComputeSafety(track.GetPosition());
  if(d>1*mm)
    return 1*mm;
  else if (d>1*micrometer)
    return d;
  else
    return 1*pc;
}
