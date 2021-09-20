/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: B. Monteagudo  contact address: monteagu@frib.msu.edu    *                    
 *                                                                           *
 * Creation Date  : May 2020                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Sweeper simulation                                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: Sweeper setup (NSCL): 2 drift chambers (CRDC1 & CRDC2) for       *
 * tracking, an ionization chamber (dE) and a scintillator (ToF)             *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4AutoDelete.hh"

// NPTool header
#include "Sweeper.hh"
#include "MagField.hh"

#include "SweeperScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Sweeper_NS{
  
  // Mother Volume
  const double Phi = 0*deg;
  const double VolumeWidth = 1.000*m;
  const double VolumeLength = 0.800*m;
  const double VolumeThickness = 3.5*m;
  //Sweeper
  const double Rmin = 0.53*m ; //1.03-0.5=0.53
  const double Rmax = 1.53*m ; //1.03+0.5=1.53
  const double Rmin_MagField = 0.83*m ; //1.03+0.2=1.23
  const double Rmax_MagField = 1.23*m ; //1.03-0.2=0.83
  const double Theta = 43.3*deg ;
  const double Width = 100*mm ;
  const double Thickness = 25*cm ;
  const string Material = "Fe";

  const double Thickness_MagField = 14*cm ;
  // const double Thickness_MagField = 50*cm ;
  
}
namespace CRDC_NS{
  const double ResoDriftTime = 0.0001*ns ; //? 
  const double ResoPosX = 0.5*mm ;
  const double ResoPosY = 0.5*mm ;
  const double Width = 30*cm ;
  const double Thickness = 10*cm ;
  const G4double DriftSpeed = 6 * cm / microsecond;//?
  const G4ThreeVector DriftDir   = G4ThreeVector(0, 1, 0);
  const string Material[2] = {"Isobutane", "CF4"};
  const int GasProportion[2] = {20, 80};
  
}
namespace IonChamber_NS{
  const double EnergyThreshold = 0.1*MeV;//?
  const double ResoEnergy = 3*perCent;//%
  const double Width = 45*cm ;
  const double Thickness = 65*cm ;
  const string Material = "P10_1atm";

}
namespace ThinScint_NS{
  const double EnergyThreshold = 0.001*MeV;//?
  const double ResoTime = 0.300*ns ;
  const double ResoEnergy = 0.001*MeV ;//?
  const double Width = 55*cm ;
  const double Thickness = 5*cm ;
  const string Material = "BC400";//?
}
namespace OldHodo_NS{
  const double EnergyThreshold = 0.001*MeV;//?
  //const double ResoTime = 0.94*ns ;
  const double ResoEnergy = 0.001*MeV ;//?
  const double Width = 412.5*cm ;
  const double Thickness = 55*cm ;
  const string Material = "CsI_Scintillator";//?
}
namespace NewHodo_NS{
  const double EnergyThreshold = 0.001*MeV;//?
  //const double ResoTime = 0.94*ns ;
  const double ResoEnergy = 0.03; //3% resolution
  const double Width = 250*cm ; //50*5
  const double Thickness = 30*cm ;
  const string Material = "CsI_Scintillator";
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Sweeper Specific Method
Sweeper::Sweeper(){
  m_Event = new TSweeperData() ;
  m_SweeperScorer = 0;
  // m_DriftChamberScorer=0;
  m_SweeperLog = 0;
  m_MotherLog = 0;
  m_CRDCLog=0;
  m_IonChamberLog=0;
  m_ThinScintLog=0;
  m_OldHodoLog=0;
  m_NewHodoLog=0;
  m_SweeperPhys = 0;
  m_SweeperMagFieldLog = 0;
  
  
  // RGB Color + Transparency
 
  m_VisSweeper = new G4VisAttributes(G4Colour(0, 0, 1, 0.5)); //Blue   
  m_VisCRDC = new G4VisAttributes(G4Colour(0, 1, 0, 0.5)); //Green 
  m_VisIonChamber = new G4VisAttributes(G4Colour(1, 0, 1, 0.5)); //Magenta 
  m_VisThinScint = new G4VisAttributes(G4Colour(1, 1, 0, 0.5)); //Yellow
  m_VisHodo = new G4VisAttributes(G4Colour(1, 0, 0, 0.5)); //Red
  
}

Sweeper::~Sweeper(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Sweeper::AddDetector(G4ThreeVector POS, double  Theta, double  Brho, double *Dist){
  m_Pos.push_back(POS);
  m_R.push_back(POS.mag());
  m_Theta.push_back(Theta);
  m_Brho.push_back(Brho);
  
  //Dist
  m_DistToExit.push_back(Dist[0]);
  m_DistToDC1.push_back(Dist[1]);
  m_DistToDC2.push_back(Dist[2]);
  m_DistToIC.push_back(Dist[3]);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Sweeper::BuildMotherVolume(){
  if(!m_MotherLog){
    G4Box* box = new G4Box("MotherVolume",Sweeper_NS::VolumeWidth*0.5,Sweeper_NS::VolumeLength*0.5,Sweeper_NS::VolumeThickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    m_MotherLog = new G4LogicalVolume(box,DetectorMaterial,"logic_Mother",0,0,0);
    m_MotherLog->SetVisAttributes(G4VisAttributes::GetInvisible);
    //m_MotherLog->SetSensitiveDetector(m_SweeperScorer);

  }
  return m_MotherLog;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Sweeper::BuildCRDC(){
  if(!m_CRDCLog){
    G4Box* box = new G4Box("CRDC",CRDC_NS::Width*0.5,CRDC_NS::Width*0.5,CRDC_NS::Thickness*0.5);

    G4Material* DetectorMaterial = new G4Material("CRDC_Material",1.6347*mg/cm3,2);
    int NumMaterial = sizeof(CRDC_NS::GasProportion)/sizeof(CRDC_NS::GasProportion[0]);
    
    for(int i=0;i<NumMaterial;i++)
      {  G4Material *material = MaterialManager::getInstance()->GetGasFromLibrary(CRDC_NS::Material[i],1.01325*bar,288.15*kelvin);
        DetectorMaterial->AddMaterial(material,CRDC_NS::GasProportion[i]*perCent);
      }
    
    m_CRDCLog = new G4LogicalVolume(box,DetectorMaterial,"logic_CRDC",0,0,0);
    m_CRDCLog->SetVisAttributes(m_VisCRDC);

    m_CRDCLog->SetSensitiveDetector(m_SweeperScorer);

  }
  return m_CRDCLog;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Sweeper::BuildNewHodo(){
  if(!m_NewHodoLog){
    G4Box* box = new G4Box("NewHodo",NewHodo_NS::Width*0.5,NewHodo_NS::Width*0.5,NewHodo_NS::Thickness*0.5);
    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(NewHodo_NS::Material);
 
    m_NewHodoLog = new G4LogicalVolume(box,DetectorMaterial,"logic_NewHodo",0,0,0);
    m_NewHodoLog->SetVisAttributes(m_VisHodo);
    m_NewHodoLog->SetSensitiveDetector(m_SweeperScorer);
  }
  return m_NewHodoLog;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Sweeper::BuildOldHodo(){
  if(!m_OldHodoLog){
    G4Box* box = new G4Box("OldHodo",OldHodo_NS::Width*0.5,OldHodo_NS::Width*0.5,OldHodo_NS::Thickness*0.5);
   G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(NewHodo_NS::Material);
   
    m_OldHodoLog = new G4LogicalVolume(box,DetectorMaterial,"logic_OldHodo",0,0,0);
    m_OldHodoLog->SetVisAttributes(m_VisHodo);
    m_OldHodoLog->SetSensitiveDetector(m_SweeperScorer);

  }
  return m_OldHodoLog;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Sweeper::BuildSweeper(double theta){
  if(!m_SweeperLog){
    G4Tubs* tub = new G4Tubs("Sweeper_Cyl",Sweeper_NS::Rmin,Sweeper_NS::Rmax,Sweeper_NS::Thickness*0.5,0,/*Sweeper_NS::Theta*/theta);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    m_SweeperLog = new G4LogicalVolume(tub,DetectorMaterial,"logic_SweeperMagField",0,0,0);
    m_SweeperLog->SetVisAttributes(m_VisSweeper);
    //m_SweeperLog->SetSensitiveDetector(m_SweeperScorer);

  }
  return m_SweeperLog;
}
G4LogicalVolume* Sweeper::BuildSweeperMagField(double theta){
  if(!m_SweeperMagFieldLog){
    G4Tubs* tub = new G4Tubs("SweeperMagFieldCyl",Sweeper_NS::Rmin_MagField,Sweeper_NS::Rmax_MagField,Sweeper_NS::Thickness_MagField*0.5,0,/*Sweeper_NS::Theta*/ theta);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    m_SweeperMagFieldLog = new G4LogicalVolume(tub,DetectorMaterial,"logic_SweeperMagField",0,0,0);
    m_SweeperMagFieldLog->SetVisAttributes(/*G4VisAttributes::Invisible*/ m_VisHodo);
    //m_SweeperMagFieldLog->SetSensitiveDetector(m_SweeperMagFieldScorer);
    
  }
  return m_SweeperMagFieldLog;
}

G4LogicalVolume* Sweeper::BuildIonChamber(){

  if(!m_IonChamberLog){
    G4Box* box = new G4Box("IonChamber_Box",IonChamber_NS::Width*0.5,IonChamber_NS::Width*0.5,IonChamber_NS::Thickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(IonChamber_NS::Material);
    
    m_IonChamberLog = new G4LogicalVolume(box,DetectorMaterial,"logic_IonChamber",0,0,0);
    m_IonChamberLog->SetVisAttributes(m_VisIonChamber);
    m_IonChamberLog->SetSensitiveDetector(m_SweeperScorer);

  }
  return m_IonChamberLog;
											  
  
}
G4LogicalVolume* Sweeper::BuildThinScint(){

  if(!m_ThinScintLog){
    G4Box* box = new G4Box("ThinScint_Box",ThinScint_NS::Width*0.5,ThinScint_NS::Width*0.5,ThinScint_NS::Thickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ThinScint_NS::Material);
    
    m_ThinScintLog = new G4LogicalVolume(box,DetectorMaterial,"logic_ThinScint",0,0,0);
    m_ThinScintLog->SetVisAttributes(m_VisThinScint);
    m_ThinScintLog->SetSensitiveDetector(m_SweeperScorer);

  }
  return m_ThinScintLog;
											  
  
}

void Sweeper::SetSweeperField(bool kMap, double B_Field=0){
  
  //G4FieldManager* SweeperMagFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager(); //?? not working

  G4FieldManager* SweeperMagFieldMgr;

  if(kMap){ //3D map
    MagField *FieldMap = new MagField();
    FieldMap->LoadMagneticField(/*filename*/"MFmap/mapfile.csv");
    FieldMap->SetMagAngle(/*fAngle*/43.3*deg);
    FieldMap->SetFieldFactor(1.);

    SweeperMagFieldMgr = new G4FieldManager(FieldMap);
    SweeperMagFieldMgr->SetDetectorField(FieldMap);
    SweeperMagFieldMgr->CreateChordFinder(FieldMap);

    G4AutoDelete::Register(FieldMap);
    G4AutoDelete::Register(SweeperMagFieldMgr);
  }else { //Uniform map
   
    cout<<"Uniform magnetic field set for: "<<B_Field/tesla<<" "<<B_Field*tesla<<endl;
    
    fSweeperMagField = new G4UniformMagField(G4ThreeVector(0.,B_Field*tesla,0.));
    SweeperMagFieldMgr = new G4FieldManager(fSweeperMagField);
    SweeperMagFieldMgr->SetDetectorField(fSweeperMagField);
    SweeperMagFieldMgr->CreateChordFinder(fSweeperMagField);
    //G4AutoDelete::Register(fSweeperMagField);
  }


  //step for magnetic field integration ???
  //SweeperMagFieldMgr->GetChordFinder()->SetDeltaChord(0.001*mm); // 10 micron m order accuracy

  m_SweeperMagFieldLog->SetFieldManager(SweeperMagFieldMgr, true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Sweeper::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Sweeper");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> basic = {"Pos","Theta","Brho"};
  vector<string> dist_par = {"DistToExit", "DistToDC1", "DistToDC2", "DistToIC"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
     if(blocks[i]->HasTokenList(basic)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Sweeper setup " << i+1 <<  endl;
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("Pos","mm"));
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Brho = blocks[i]->GetDouble("Brho","tesla*m");
      //double Phi = blocks[i]->GetDouble("Phi");

      //Get distance parameters
      double Dist[]={-1,-1,-1,-1};
      for(unsigned int j=0; j<dist_par.size(); j++){
	if(blocks[i]->HasToken(dist_par[j])){
	  Dist[j]=blocks[i]->GetDouble(dist_par[j], "mm");
	}
      }

      AddDetector(Pos,Theta,Brho,Dist);      
      
     }
   
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Sweeper::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {
    
    ///// Sweeper
    G4double Sweeper_R_min = Sweeper_NS::Rmin; //0.53*m; // 1.03 - 0.2
    G4double Sweeper_R_max = Sweeper_NS::Rmax; //1.53*m; // 1.03 + 0.2
    G4double Sweeper_x=-(Sweeper_R_min+Sweeper_R_max)/2;
    G4double Sweeper_y = Sweeper_NS::Thickness_MagField*0.5 + Sweeper_NS::Thickness*0.5; //7cm ?
    G4double Sweeper_z =  m_R[i];
    G4ThreeVector Sweeper_pos = G4ThreeVector(Sweeper_x,Sweeper_y,Sweeper_z);
    
    G4RotationMatrix *RotSweeper = new G4RotationMatrix();
    RotSweeper->rotateX(90*deg);
    //RotSweeper->rotateY(-90*deg);
    
    // Upper part of sweeper in (x,y,z)
    m_SweeperPhys = new G4PVPlacement(G4Transform3D(*RotSweeper,Sweeper_pos),
      BuildSweeper(m_Theta[i]),
      "Sweeper",world,false,0);
    //Lower part of sweeper in (x,-y,z)
    m_SweeperPhys = new G4PVPlacement(G4Transform3D(*RotSweeper,G4ThreeVector(Sweeper_x,-Sweeper_y,Sweeper_z)),
    				      BuildSweeper(m_Theta[i]),
    				      "sweeper_phys",world,false,0);
    ///// Magnetic Field Sweeper

    G4double SweeperMagField_x = Sweeper_x;
    G4double SweeperMagField_y = 0.0;
    G4double SweeperMagField_z = m_R[i];
    G4ThreeVector SweeperMagField_pos = G4ThreeVector(SweeperMagField_x, SweeperMagField_y,SweeperMagField_z);

    new G4PVPlacement(G4Transform3D(*RotSweeper,SweeperMagField_pos),
    				      BuildSweeperMagField(m_Theta[i]),"sweeperMagField_phys",world,false,0);

    /////////////// MAGNETIC FIELD

    
    G4double Sweeper_ArcLength = ((Sweeper_NS::Rmin_MagField+Sweeper_NS::Rmax_MagField)/2.)*CLHEP::pi*m_Theta[i]/(180*deg);
    //G4double B_Field = (m_Brho[i]*tesla*m)/(Sweeper_ArcLength/m);
    G4double B_Field = (m_Brho[i]*tesla*m)/(((Sweeper_NS::Rmin_MagField/m+Sweeper_NS::Rmax_MagField/m)/2.));

    cout<<m_Brho[i]<<" "<<Sweeper_ArcLength/m<<" "<<B_Field/tesla<<" "<< (Sweeper_NS::Rmin_MagField+Sweeper_NS::Rmax_MagField)/2 << " "<< CLHEP::pi*m_Theta[i]/(180*deg)<<endl;
      
    bool k3DMap = false;
    SetSweeperField(k3DMap, B_Field); //3D case or uniform
    
    //// Volume containing detectors of sweeper ///////////
    
    G4double wX = Sweeper_x-sin(m_Theta[i])*(Sweeper_NS::VolumeThickness*0.5)-cos(m_Theta[i])*Sweeper_x/*Sweeper_NS::Rmin+Sweeper_NS::VolumeWidth*0.5)*/;
    G4double wY = 0;
    G4double wZ =Sweeper_z+cos(m_Theta[i])*(Sweeper_NS::VolumeThickness*0.5)-sin(m_Theta[i])*Sweeper_x/*+Sweeper_NS::Rmin-Sweeper_NS::VolumeWidth*0.5)*/; 
  
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    
    G4RotationMatrix *RotMother = new G4RotationMatrix();
    RotMother->rotateY(-m_Theta[i]);

    new G4PVPlacement(G4Transform3D(*RotMother,Det_pos),BuildMotherVolume(),"MotherVolume",world, false,0);
    
    //// CRDCs

    G4double CRDC1_x = 0.0;
    G4double CRDC1_y = 0.0;
    G4double CRDC1_z = -Sweeper_NS::VolumeThickness*0.5+47*cm;
    cout<<CRDC1_z<<endl;
    G4ThreeVector CRDC1_pos = G4ThreeVector(CRDC1_x,CRDC1_y,CRDC1_z);
    m_CRDCPhys = new G4PVPlacement(0,CRDC1_pos,
    				      BuildCRDC(),"crdc1_phys",m_MotherLog,false,0);

    G4double Dist_CRDCs;
    if(m_DistToExit[i]!=-1)Dist_CRDCs=m_DistToExit[i];
    else Dist_CRDCs=1.54*m;
    
    G4double CRDC2_x = 0.0;
    G4double CRDC2_y = 0.0;
    G4double CRDC2_z = CRDC1_z + Dist_CRDCs;

    G4ThreeVector CRDC2_pos = G4ThreeVector(CRDC2_x,CRDC2_y,CRDC2_z);
    m_CRDCPhys = new G4PVPlacement(0,CRDC2_pos,
     				      BuildCRDC(),"crdc2_phys",m_MotherLog,false,1);

    ///// Ion Chamber
    G4double Dist_CRDC2ToIonChamber;
    if(m_DistToDC1[i]!=-1) Dist_CRDC2ToIonChamber=m_DistToDC1[i];
    else  Dist_CRDC2ToIonChamber=0.09*m;
    
    G4double IonChamber_x = 0.0;
    G4double IonChamber_y = 0.0;
    G4double IonChamber_z = CRDC2_z + (CRDC_NS::Thickness*0.5 + IonChamber_NS::Thickness*0.5 + Dist_CRDC2ToIonChamber);
    G4ThreeVector IonChamber_pos = G4ThreeVector(IonChamber_x,IonChamber_y,IonChamber_z);
    new G4PVPlacement(0,IonChamber_pos,
		      BuildIonChamber(),"IonChamber_phys",m_MotherLog,false,2);

    ////// Thin Scint
    G4double Dist_IonChamberToScint;
    if(m_DistToDC2[i]!=-1) Dist_IonChamberToScint=m_DistToDC2[i];
    else Dist_IonChamberToScint=9.525*cm;
    
    G4double Thin_x = 0;
    G4double Thin_y = 0;
    G4double Thin_z = IonChamber_z+IonChamber_NS::Thickness*0.5+ThinScint_NS::Thickness*0.5+Dist_IonChamberToScint;

    
    G4ThreeVector ThinScint_pos = G4ThreeVector(Thin_x,Thin_y,Thin_z);
    new G4PVPlacement(0,ThinScint_pos,
    		      BuildThinScint(),"ThinScint_phys",m_MotherLog,false,3);

    //////// New Hodoscope
    // G4double Dist_ScintToNewHodo=100*cm;
    // if(m_DistToIC[i]!=-1)Dist_ScintToNewHodo=m_DistToIC[i];
    // else Dist_ScintToNewHodo=100*cm;
    
    // G4double NewHodo_x = 0;
    // G4double NewHodo_y = 0;
    // G4double NewHodo_z = Dist_ScintToNewHodo + ThinScint_NS::Thickness*0.5 + NewHodo_NS::Thickness*0.5;
    
    // G4ThreeVector NewHodoScint_pos = G4ThreeVector(NewHodo_x,NewHodo_y,NewHodo_z);
    // new G4PVPlacement(0,NewHodoScint_pos,
    // 		      BuildNewHodo(),"NewHodoScint_phys",m_MotherLog,false,5);

    // //////// Old Hodoscope
    // G4double Dist_ScintToOldHodo=50*cm;
    // if(m_DistToIC[i]!=-1)Dist_ScintToNewHodo=m_DistToIC[i];
    // else Dist_ScintToNewHodo=50*cm;
    
    // 
    // G4double OldHodo_x = 0;
    // G4double OldHodo_y = 0;
    // G4double OldHodo_z = Dist_ScintToOldHodo + NewHodo_NS::Thickness*0.5 + OldHodo_NS::Thickness*0.5;
    
    // G4ThreeVector OldHodoScint_pos = G4ThreeVector(OldHodo_x,OldHodo_y,OldHodo_z);
    // new G4PVPlacement(0,OldHodoScint_pos,
    // 		      BuildOldHodo(),"OldHodoScint_phys",m_MotherLog,false,6);

    
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Sweeper::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Sweeper")){
    pTree->Branch("Sweeper", "TSweeperData", &m_Event) ;
  }
  pTree->SetBranchAddress("Sweeper", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Sweeper::ReadSensitive(const G4Event* event){
  m_Event->Clear();

  ////// Sweeper Collection Hits
  NPS::HitsMap<G4double*>*     SweeperHitMap;
  std::map<G4int, G4double**>::iterator    Sweeper_itr;
  G4int SweeperCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("SweeperScorer/SweeperHits");

  if(SweeperCollectionID == -1){
    G4cerr<< " ERROR "<<G4endl;
    return;
  }
  
  SweeperHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(SweeperCollectionID));
  
  //Loop in map
  for(Sweeper_itr = SweeperHitMap->GetMap()->begin(); Sweeper_itr != SweeperHitMap->GetMap()->end();Sweeper_itr++){
    
    G4double *Info = *(Sweeper_itr->second); 
    cout<<"Energy: "<< Info[0]<<endl;
    cout<<"time: "<< Info[1]<<endl;
    cout<<"Detector Nbr: "<< Info[7]<<" "<<endl;

    double energy = Info[0];
    double time = Info[1];
    double xpos = Info[2];
    double ypos = Info[3];
    unsigned short detnum = Info[7];

    if(detnum<2){//DCs
      double x = RandGauss::shoot(xpos, CRDC_NS::ResoPosX);
      double y = RandGauss::shoot(ypos, CRDC_NS::ResoPosY);
  
      m_Event->SetPosition(detnum,x,y);

    }else if(detnum==2){//IC
      double e = RandGauss::shoot(energy, IonChamber_NS::ResoEnergy*energy);
      m_Event->SetEnergy(detnum,e);
    }else if(detnum==3){//Thin
      double t = RandGauss::shoot(time, ThinScint_NS::ResoTime);
      m_Event->SetTime(detnum,t);
    }
    //else if(detnum==5){ //Hodos

  }
  cout<<SweeperHitMap->GetSize()<<endl;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Sweeper::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 

  m_SweeperScorer = CheckScorer("SweeperScorer",already_exist);
  //m_DriftChamberScorer = CheckScorer("DriftChamberScorer",already_exist);
  
  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* SweeperHits = new SweeperScorers::PS_Sweeper("SweeperHits", level, 0);
  
  //and register it to the multifunctionnal detector
  m_SweeperScorer->RegisterPrimitive(SweeperHits);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SweeperScorer) ;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Sweeper::Construct(){
  return  (NPS::VDetector*) new Sweeper();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Sweeper{
    public:
      proxy_nps_Sweeper(){
        NPS::DetectorFactory::getInstance()->AddToken("Sweeper","Sweeper");
        NPS::DetectorFactory::getInstance()->AddDetector("Sweeper",Sweeper::Construct);
      }
  };

  proxy_nps_Sweeper p_nps_Sweeper;
}
