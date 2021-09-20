/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 15/07/09                                                 *
 * Last update    : 12/10/09                                                 *
 *---------------------------------------------------------------------------*
 * Decription: Define a module of trapezoidal shape for the Gaspard tracker  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + 12/10/09: Change scorer scheme (N. de Sereville)                     *
 *    + 01/10/10: Fix bug with TInteractionCoordinate map size in Read       *
 *                Sensitive (N. de Sereville)                                *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <string>
#include <cmath>

// G4 Geometry headers
#include "G4Box.hh"
#include "G4Trap.hh"

// G4 various headers
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4PVDivision.hh"

// G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool headers
#include "GaspardTrackerTrapezoid.hh"
#include "ObsoleteGeneralScorers.hh"
#include "GaspardScorers.hh"
#include "RootOutput.h"
#include "NPSVDetector.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"

// CLHEP
#include "CLHEP/Random/RandGauss.h"
using namespace std;
using namespace CLHEP;

using namespace GPDTRAP;
using namespace GPDSCORERS;




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GaspardTrackerTrapezoid::GaspardTrackerTrapezoid()
{
   ms_InterCoord = 0;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GaspardTrackerTrapezoid::~GaspardTrackerTrapezoid()
{
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GaspardTrackerTrapezoid::AddModule(G4ThreeVector X1_Y1     ,
      G4ThreeVector X128_Y1   ,
      G4ThreeVector X1_Y128   ,
      G4ThreeVector X128_Y128 ,
      bool wFirstStage        ,
      bool wSecondStage       ,
      bool wThirdStage)
{
   m_DefinitionType.push_back(true) ;

   m_X1_Y1.push_back(X1_Y1)               ;
   m_X128_Y1.push_back(X128_Y1)           ;
   m_X1_Y128.push_back(X1_Y128)           ;
   m_X128_Y128.push_back(X128_Y128)       ;
   m_wFirstStage.push_back(wFirstStage)   ;
   m_wSecondStage.push_back(wSecondStage) ;
   m_wThirdStage.push_back(wThirdStage)   ;

   m_R.push_back(0)      ;
   m_Theta.push_back(0)  ;
   m_Phi.push_back(0)    ;
   m_beta_u.push_back(0) ;
   m_beta_v.push_back(0) ;
   m_beta_w.push_back(0) ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GaspardTrackerTrapezoid::AddModule(G4double R        ,
      G4double Theta    ,
      G4double Phi      ,
      G4double beta_u   ,
      G4double beta_v   ,
      G4double beta_w   ,
      bool wFirstStage  ,
      bool wSecondStage ,
      bool wThirdStage)
{
   G4ThreeVector empty = G4ThreeVector(0, 0, 0);

   m_DefinitionType.push_back(false);

   m_R.push_back(R)                       ;
   m_Theta.push_back(Theta)               ;
   m_Phi.push_back(Phi)                   ;
   m_beta_u.push_back(beta_u)             ;
   m_beta_v.push_back(beta_v)             ;
   m_beta_w.push_back(beta_w)             ;
   m_wFirstStage.push_back(wFirstStage)   ;
   m_wSecondStage.push_back(wSecondStage) ;
   m_wThirdStage.push_back(wThirdStage)   ;

   m_X1_Y1.push_back(empty)     ;
   m_X128_Y1.push_back(empty)   ;
   m_X1_Y128.push_back(empty)   ;
   m_X128_Y128.push_back(empty) ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GaspardTrackerTrapezoid::VolumeMaker(G4int DetectorNumber,
                                          G4ThreeVector MMpos,
                                          G4RotationMatrix* MMrot,
                                          bool wFirstStage,
                                          bool wSecondStage,
                                          bool wThirdStage,
                                          G4LogicalVolume* world)
{
   G4double NbrTelescopes = DetectorNumber  ;
   G4String DetNumber                   ;
   ostringstream Number                      ;
   Number << NbrTelescopes                   ;
   DetNumber = Number.str()             ;

   ////////////////////////////////////////////////////////////////
   ////////////// Starting Volume Definition //////////////////////
   ////////////////////////////////////////////////////////////////
   G4String Name = "GPDTrapezoid" + DetNumber ;

   // Definition of the volume containing the sensitive detector
   G4Trap* solidGPDTrapezoid = new G4Trap(Name, 
                                          Length/2, 0*deg, 0*deg, 
                                          Height/2, BaseLarge/2, BaseSmall/2, 0*deg, 
                                          Height/2, BaseLarge/2, BaseSmall/2, 0*deg);
   G4LogicalVolume* logicGPDTrapezoid = new G4LogicalVolume(solidGPDTrapezoid, m_MaterialVacuum, Name, 0, 0, 0);

   new G4PVPlacement(G4Transform3D(*MMrot, MMpos), logicGPDTrapezoid, Name, world, false, DetectorNumber);

   G4VisAttributes* TrapezoideVisAtt = new G4VisAttributes(G4Colour(0.90, 0.90, 0.90));
   TrapezoideVisAtt->SetForceWireframe(true); 
    logicGPDTrapezoid->SetVisAttributes(TrapezoideVisAtt);

   //Place two marker to identify the u and v axis on silicon face:
   //marker are placed a bit before the silicon itself so they don't perturbate simulation
   //Uncomment to help debugging or if you want to understand the way the code work.
   //I should recommand to Comment it during simulation to avoid perturbation of simulation
   //Remember G4 is limitationg step on geometry constraints.
  /* 
         G4ThreeVector positionMarkerU = CT*0.98 + MMu*SiliconFace/4;
         G4Box*          solidMarkerU = new G4Box( "solidMarkerU" , SiliconFace/4 , 1*mm , 1*mm )              ;
         G4LogicalVolume* logicMarkerU = new G4LogicalVolume( solidMarkerU , m_MaterialVacuum , "logicMarkerU",0,0,0)       ;
         PVPBuffer = new G4PVPlacement(G4Transform3D(*MMrot,positionMarkerU),logicMarkerU,"MarkerU",world,false,0) ;

         G4VisAttributes* MarkerUVisAtt= new G4VisAttributes(G4Colour(0.,0.,0.5));//blue
         logicMarkerU->SetVisAttributes(MarkerUVisAtt);

         G4ThreeVector positionMarkerV = CT*0.98 + MMv*SiliconFace/4;
         G4Box*          solidMarkerV = new G4Box( "solidMarkerU" , 1*mm , SiliconFace/4 , 1*mm )              ;
         G4LogicalVolume* logicMarkerV = new G4LogicalVolume( solidMarkerV , m_MaterialVacuum , "logicMarkerV",0,0,0)       ;
         PVPBuffer = new G4PVPlacement(G4Transform3D(*MMrot,positionMarkerV),logicMarkerV,"MarkerV",world,false,0) ;

         G4VisAttributes* MarkerVVisAtt= new G4VisAttributes(G4Colour(0.,0.5,0.5));//green
         logicMarkerV->SetVisAttributes(MarkerVVisAtt);
   */

   ////////////////////////////////////////////////////////////////
   /////////////////// First Stage Construction////////////////////
   ////////////////////////////////////////////////////////////////
   if (wFirstStage) {
      // Silicon detector itself
      G4ThreeVector  positionFirstStage = G4ThreeVector(0, 0, FirstStage_PosZ);

      G4Trap* solidFirstStage = new G4Trap("solidFirstStage", 
                                           FirstStageThickness/2, 0*deg, 0*deg, 
                                           FirstStageHeight/2, FirstStageBaseLarge/2, FirstStageBaseSmall/2, 0*deg, 
                                           FirstStageHeight/2, FirstStageBaseLarge/2, FirstStageBaseSmall/2, 0*deg);
      G4LogicalVolume* logicFirstStage = new G4LogicalVolume(solidFirstStage, m_MaterialSilicon, "logicFirstStage", 0, 0, 0);

      new G4PVPlacement(0,
                        positionFirstStage,
                        logicFirstStage,
                        Name + "_FirstStage",
                        logicGPDTrapezoid,
                        false,
                        DetectorNumber);

      // Set First Stage sensible
      logicFirstStage->SetSensitiveDetector(m_FirstStageScorer);

      ///Visualisation of FirstStage Strip
      G4VisAttributes* FirstStageVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));   // blue
      logicFirstStage->SetVisAttributes(FirstStageVisAtt);
   }

   ////////////////////////////////////////////////////////////////
   //////////////// Second Stage  Construction ////////////////////
   ////////////////////////////////////////////////////////////////
   if (wSecondStage) {
      // Second stage silicon detector
      G4ThreeVector  positionSecondStage = G4ThreeVector(0, 0, SecondStage_PosZ);

      G4Trap* solidSecondStage = new G4Trap("solidSecondStage", 
                                            SecondStageThickness/2, 0*deg, 0*deg, 
                                           FirstStageHeight/2, FirstStageBaseLarge/2, FirstStageBaseSmall/2, 0*deg, 
                                           FirstStageHeight/2, FirstStageBaseLarge/2, FirstStageBaseSmall/2, 0*deg);
      G4LogicalVolume* logicSecondStage = new G4LogicalVolume(solidSecondStage, m_MaterialSilicon, "logicSecondStage", 0, 0, 0);

      new G4PVPlacement(0,
                        positionSecondStage,
                        logicSecondStage,
                        Name + "_SecondStage",
                        logicGPDTrapezoid,
                        false,
                        DetectorNumber);

      // Set Second Stage sensible
      logicSecondStage->SetSensitiveDetector(m_SecondStageScorer);

      ///Visualisation of SecondStage Strip
      G4VisAttributes* SecondStageVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
      logicSecondStage->SetVisAttributes(SecondStageVisAtt);
   }

   ////////////////////////////////////////////////////////////////
   ///////////////// Third Stage Construction /////////////////////
   ////////////////////////////////////////////////////////////////
   if (wThirdStage) {
      // Third stage silicon detector
      G4ThreeVector  positionThirdStage = G4ThreeVector(0, 0, ThirdStage_PosZ);

      G4Trap* solidThirdStage = new G4Trap("solidThirdStage", 
                                           ThirdStageThickness/2, 0*deg, 0*deg, 
                                           FirstStageHeight/2, FirstStageBaseLarge/2, FirstStageBaseSmall/2, 0*deg, 
                                           FirstStageHeight/2, FirstStageBaseLarge/2, FirstStageBaseSmall/2, 0*deg);
      G4LogicalVolume* logicThirdStage = new G4LogicalVolume(solidThirdStage, m_MaterialSilicon, "logicThirdStage", 0, 0, 0);

      new G4PVPlacement(0,
                        positionThirdStage,
                        logicThirdStage,
                        Name + "_ThirdStage",
                        logicGPDTrapezoid,
                        false,
                        DetectorNumber);

      // Set Third Stage sensible
      logicThirdStage->SetSensitiveDetector(m_ThirdStageScorer);

      ///Visualisation of Third Stage
      G4VisAttributes* ThirdStageVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));   // red
      logicThirdStage->SetVisAttributes(ThirdStageVisAtt);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void GaspardTrackerTrapezoid::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("GaspardTracker");
  vector<string> token_cart= {"X1_Y1","X128_Y1","X1_Y128","X128_Y128"};
  vector<string> token_sphe= {"R","THETA","PHI","BETA"};

  vector<string> token={"FIRSTSTAGE","SECONDSTAGE","THIRDSTAGE"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->GetMainValue() == "Trapezoid" && blocks[i]->HasTokenList(token) ){
      cout << "Gaspard Trapezoid " << i+1 << ":"  << endl; 

       
      bool first = blocks[i]->GetInt("FIRSTSTAGE");
      bool second = blocks[i]->GetInt("SECONDSTAGE");
      bool third = blocks[i]->GetInt("THIRDSTAGE");
      if(blocks[i]->HasToken("VIS"))
        m_non_sensitive_part_visiualisation =  blocks[i]->GetInt("VIS");

      if(blocks[i]->HasTokenList(token_cart)){
        // Add module
        G4ThreeVector A = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1","mm"));
        G4ThreeVector B = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y1","mm"));
        G4ThreeVector C = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y128","mm"));
        G4ThreeVector D = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y128","mm"));
      
        AddModule(A,B,C,D,first,second,third);
      }
     else if(blocks[i]->HasTokenList(token_sphe)){
        // Add module
        double R = blocks[i]->GetDouble("R","mm");
        double Theta = blocks[i]->GetDouble("THETA","deg");
        double Phi = blocks[i]->GetDouble("PHI","deg");
        vector<double> beta = blocks[i]->GetVectorDouble("BETA","deg");
      
        AddModule(R,Theta,Phi,beta[0],beta[1],beta[2],first,second,third);
      }
    }
  } 
}

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void GaspardTrackerTrapezoid::ConstructDetector(G4LogicalVolume* world)
{
   G4RotationMatrix* MMrot    = NULL;
   G4ThreeVector     MMpos    = G4ThreeVector(0, 0, 0);
   G4ThreeVector     MMu      = G4ThreeVector(0, 0, 0);
   G4ThreeVector     MMv      = G4ThreeVector(0, 0, 0);
   G4ThreeVector     MMw      = G4ThreeVector(0, 0, 0);
   G4ThreeVector     MMCenter = G4ThreeVector(0, 0, 0);

   bool FirstStage  = true ;
   bool SecondStage = true ;
   bool ThirdStage  = true ;

   G4int NumberOfModule = m_DefinitionType.size() ;

   for (G4int i = 0; i < NumberOfModule; i++) {
      // By Point
      if (m_DefinitionType[i]) {
         // (u,v,w) unitary vector associated to trapezoidal referencial
         // (u,v) // to silicon plan
         //      -------
         //     /       \              ^
         //    /         \             |  v
         //   /           \            |
         //  ---------------     <------
         //                         u
         // w perpendicular to (u,v) plan and pointing ThirdStage
         MMu = m_X128_Y1[i] - m_X1_Y1[i];
         MMu = MMu.unit();

         MMv = 0.5 * (m_X1_Y128[i] + m_X128_Y128[i] - m_X1_Y1[i] - m_X128_Y1[i]);
         MMv = MMv.unit();

         MMw = MMu.cross(MMv);
         MMw = MMw.unit();

         // Center of the module
         MMCenter = (m_X1_Y1[i] + m_X1_Y128[i] + m_X128_Y1[i] + m_X128_Y128[i]) / 4;

         // Passage Matrix from Lab Referential to Module Referential
         MMrot = new G4RotationMatrix(MMu, MMv, MMw);
         // translation to place Module
         MMpos = MMw * Length * 0.5 + MMCenter;
      }

      // By Angle
      else {
         G4double Theta = m_Theta[i];
         G4double Phi   = m_Phi[i];

         // (u,v,w) unitary vector associated to telescope referencial
         // (u,v) // to silicon plan
         //      -------
         //     /       \              ^
         //    /         \             |  v
         //   /           \            |
         //  ---------------     <------
         //                         u
         // w perpendicular to (u,v) plan and pointing ThirdStage
         // Phi is angle between X axis and projection in (X,Y) plan
         // Theta is angle between  position vector and z axis
         G4double wX = m_R[i] * sin(Theta) * cos(Phi);
         G4double wY = m_R[i] * sin(Theta) * sin(Phi);
         G4double wZ = m_R[i] * cos(Theta);
         MMw = G4ThreeVector(wX, wY, wZ);

         // vector corresponding to the center of the module
         MMCenter = MMw;

         // vector parallel to one axis of silicon plane
         // in fact, this is vector u
         G4double ii = cos(Theta) * cos(Phi);
         G4double jj = cos(Theta) * sin(Phi);
         G4double kk = -sin(Theta);
         G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

         MMw = MMw.unit();
         MMv = MMw.cross(Y);
         MMu = MMv.cross(MMw);
         MMv = MMv.unit();
         MMu = MMu.unit();

         // Passage Matrix from Lab Referential to Telescope Referential
         MMrot = new G4RotationMatrix(MMu, MMv, MMw);
         // Telescope is rotate of Beta angle around MMv axis.
         MMrot->rotate(m_beta_u[i], MMu);
         MMrot->rotate(m_beta_v[i], MMv);
         MMrot->rotate(m_beta_w[i], MMw);
         // translation to place Telescope
         MMpos = MMw * Length + MMCenter;

      }

      FirstStage  = m_wFirstStage[i];
      SecondStage = m_wSecondStage[i];
      ThirdStage  = m_wThirdStage[i];

      VolumeMaker(i + 1, MMpos, MMrot, FirstStage, SecondStage, ThirdStage, world);
   }

   delete MMrot;
}



// Connect the GaspardTrackingData class to the output TTree
// of the simulation
void GaspardTrackerTrapezoid::InitializeRootOutput()
{
}



// Set the TinteractionCoordinates object from NPS::VDetector to the present class
void GaspardTrackerTrapezoid::SetInterCoordPointer(TInteractionCoordinates* interCoord)
{
   ms_InterCoord = interCoord;
}



// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void GaspardTrackerTrapezoid::ReadSensitive(const G4Event* event)
{
   //////////////////////////////////////////////////////////////////////////////////////
   //////////////////////// Used to Read Event Map of detector //////////////////////////
   //////////////////////////////////////////////////////////////////////////////////////
   // First Stage
   std::map<G4int, G4int*>::iterator    DetectorNumber_itr;
   std::map<G4int, G4double*>::iterator Energy_itr;
   std::map<G4int, G4double*>::iterator Time_itr;
   std::map<G4int, G4int*>::iterator    X_itr;
   std::map<G4int, G4int*>::iterator    Y_itr;
   std::map<G4int, G4double*>::iterator Pos_X_itr;
   std::map<G4int, G4double*>::iterator Pos_Y_itr;
   std::map<G4int, G4double*>::iterator Pos_Z_itr;
   std::map<G4int, G4double*>::iterator Ang_Theta_itr;
   std::map<G4int, G4double*>::iterator Ang_Phi_itr;

   G4THitsMap<G4int>*    DetectorNumberHitMap;
   G4THitsMap<G4double>* EnergyHitMap;
   G4THitsMap<G4double>* TimeHitMap;
   G4THitsMap<G4int>*    XHitMap;
   G4THitsMap<G4int>*    YHitMap;
   G4THitsMap<G4double>* PosXHitMap;
   G4THitsMap<G4double>* PosYHitMap;
   G4THitsMap<G4double>* PosZHitMap;
   G4THitsMap<G4double>* AngThetaHitMap;
   G4THitsMap<G4double>* AngPhiHitMap;

   // NULL pointer are given to avoid warning at compilation
   // Second Stage
   std::map<G4int, G4double*>::iterator SecondStageEnergy_itr;
   G4THitsMap<G4double>* SecondStageEnergyHitMap = NULL;
   // Third Stage
   std::map<G4int, G4double*>::iterator ThirdStageEnergy_itr;
   G4THitsMap<G4double>* ThirdStageEnergyHitMap = NULL;

   // Read the Scorer associated to the first Stage
   //Detector Number
   G4int StripDetCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/DetectorNumber")    ;
   DetectorNumberHitMap = (G4THitsMap<G4int>*)(event->GetHCofThisEvent()->GetHC(StripDetCollectionID))         ;
   DetectorNumber_itr =  DetectorNumberHitMap->GetMap()->begin()                                               ;

   //Energy
   G4int StripEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/StripEnergy")   ;
   EnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(StripEnergyCollectionID))                    ;
   Energy_itr = EnergyHitMap->GetMap()->begin()                                                          ;

   //Time of Flight
   G4int StripTimeCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/StripTime")    ;
   TimeHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(StripTimeCollectionID))                        ;
   Time_itr = TimeHitMap->GetMap()->begin()                                                              ;

   //Strip Number X
   G4int StripXCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/StripNumberX")    ;
   XHitMap = (G4THitsMap<G4int>*)(event->GetHCofThisEvent()->GetHC(StripXCollectionID))                              ;
   X_itr = XHitMap->GetMap()->begin()                                                                    ;

   //Strip Number Y
   G4int StripYCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/StripNumberY")    ;
   YHitMap = (G4THitsMap<G4int>*)(event->GetHCofThisEvent()->GetHC(StripYCollectionID))                              ;
   Y_itr = YHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate X
   G4int InterCoordXCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/InterCoordX")    ;
   PosXHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordXCollectionID))                              ;
   Pos_X_itr = PosXHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Y
   G4int InterCoordYCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/InterCoordY")    ;
   PosYHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordYCollectionID))                              ;
   Pos_Y_itr = PosYHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Z
   G4int InterCoordZCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/InterCoordZ")    ;
   PosZHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordZCollectionID))                              ;
   Pos_Z_itr = PosXHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Angle Theta
   G4int InterCoordAngThetaCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/InterCoordAngTheta")    ;
   AngThetaHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordAngThetaCollectionID))                              ;
   Ang_Theta_itr = AngThetaHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Angle Phi
   G4int InterCoordAngPhiCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerGPDTrapezoid/InterCoordAngPhi")    ;
   AngPhiHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordAngPhiCollectionID))                              ;
   Ang_Phi_itr = AngPhiHitMap->GetMap()->begin()                                                                    ;

   // Read the Scorer associated to the Second and Third Stage 
   // Energy second stage
   G4int SecondStageEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("SecondStageScorerGPDTrapezoid/SecondStageEnergy")      ;
   SecondStageEnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(SecondStageEnergyCollectionID))                      ;
   SecondStageEnergy_itr = SecondStageEnergyHitMap->GetMap()->begin()                                                       ;
   // Energy third stage
   G4int ThirdStageEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("ThirdStageScorerGPDTrapezoid/ThirdStageEnergy")      ;
   ThirdStageEnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(ThirdStageEnergyCollectionID))                      ;
   ThirdStageEnergy_itr = ThirdStageEnergyHitMap->GetMap()->begin()                                                       ;

   // Check the size of different map
   G4int sizeN = DetectorNumberHitMap->entries();
   G4int sizeE = EnergyHitMap->entries();
   G4int sizeT = TimeHitMap->entries();
   G4int sizeX = XHitMap->entries();
   G4int sizeY = YHitMap->entries();

   if (sizeE != sizeT || sizeT != sizeX || sizeX != sizeY) {
      G4cout << "No match size Si Event Map: sE:"
         << sizeE << " sT:" << sizeT << " sX:" << sizeX << " sY:" << sizeY << G4endl ;
      return;
   }

   // Loop on FirstStage number
   for (G4int l = 0; l < sizeN; l++) {
      G4double N     = *(DetectorNumber_itr->second);
      G4int NTrackID =   DetectorNumber_itr->first - N;

      if (N > 0) {
         // Fill detector number
         ms_Event->SetGPDTrkFirstStageFrontEDetectorNbr(m_index["Trapezoid"] + N);
         ms_Event->SetGPDTrkFirstStageFrontTDetectorNbr(m_index["Trapezoid"] + N);
         ms_Event->SetGPDTrkFirstStageBackEDetectorNbr(m_index["Trapezoid"] + N);
         ms_Event->SetGPDTrkFirstStageBackTDetectorNbr(m_index["Trapezoid"] + N);

         // Energy
         Energy_itr = EnergyHitMap->GetMap()->begin();
         for (G4int ll = 0 ; ll < sizeE ; ll++) {
            G4int ETrackID  =   Energy_itr->first - N;
            G4double E     = *(Energy_itr->second);
            if (ETrackID == NTrackID) {
               ms_Event->SetGPDTrkFirstStageFrontEEnergy(RandGauss::shoot(E, ResoFirstStage));
               ms_Event->SetGPDTrkFirstStageBackEEnergy(RandGauss::shoot(E, ResoFirstStage));
            }
            Energy_itr++;
         }

         //  Time
         Time_itr = TimeHitMap->GetMap()->begin();
         for (G4int h = 0 ; h < sizeT ; h++) {
            G4int TTrackID  =   Time_itr->first - N;
            G4double T     = *(Time_itr->second);

            if (TTrackID == NTrackID) {
               T = RandGauss::shoot(T, ResoTimePPAC)   ;
               ms_Event->SetGPDTrkFirstStageFrontTTime(RandGauss::shoot(T, ResoTimeGpd)) ;
               ms_Event->SetGPDTrkFirstStageBackTTime(RandGauss::shoot(T, ResoTimeGpd)) ;
            }
            Time_itr++;
         }

         // X
         X_itr = XHitMap->GetMap()->begin();
         for (G4int h = 0 ; h < sizeX ; h++) {
            G4int XTrackID  =   X_itr->first - N;
            G4int X         = *(X_itr->second);
            if (XTrackID == NTrackID) {
               ms_Event->SetGPDTrkFirstStageFrontEStripNbr(X);
               ms_Event->SetGPDTrkFirstStageFrontTStripNbr(X);
            }
            X_itr++;
         }

         // Y
         Y_itr = YHitMap->GetMap()->begin()  ;
         for (G4int h = 0 ; h < sizeY ; h++) {
            G4int YTrackID  =   Y_itr->first - N;
            G4int     Y     = *(Y_itr->second);
            if (YTrackID == NTrackID) {
               ms_Event->SetGPDTrkFirstStageBackEStripNbr(Y);
               ms_Event->SetGPDTrkFirstStageBackTStripNbr(Y);
            }
            Y_itr++;
         }

         // Pos X
         Pos_X_itr = PosXHitMap->GetMap()->begin();
         for (unsigned int h = 0; h < PosXHitMap->entries(); h++) {
            G4int PosXTrackID =   Pos_X_itr->first - N    ;
            G4double PosX     = *(Pos_X_itr->second)      ;
            if (PosXTrackID == NTrackID) {
               ms_InterCoord->SetDetectedPositionX(PosX) ;
            }
            Pos_X_itr++;
         }

         // Pos Y
         Pos_Y_itr = PosYHitMap->GetMap()->begin();
         for (unsigned int h = 0; h < PosYHitMap->entries(); h++) {
            G4int PosYTrackID =   Pos_Y_itr->first  - N   ;
            G4double PosY     = *(Pos_Y_itr->second)      ;
            if (PosYTrackID == NTrackID) {
               ms_InterCoord->SetDetectedPositionY(PosY) ;
            }
            Pos_Y_itr++;
         }

         // Pos Z
         Pos_Z_itr = PosZHitMap->GetMap()->begin();
         for (unsigned int h = 0; h < PosZHitMap->entries(); h++) {
            G4int PosZTrackID =   Pos_Z_itr->first - N    ;
            G4double PosZ     = *(Pos_Z_itr->second)      ;
            if (PosZTrackID == NTrackID) {
               ms_InterCoord->SetDetectedPositionZ(PosZ) ;
            }
            Pos_Z_itr++;
         }

         // Angle Theta
         Ang_Theta_itr = AngThetaHitMap->GetMap()->begin();
         for (unsigned int h = 0; h < AngThetaHitMap->entries(); h++) {
            G4int AngThetaTrackID =   Ang_Theta_itr->first - N    ;
            G4double AngTheta     = *(Ang_Theta_itr->second)      ;
            if (AngThetaTrackID == NTrackID) {
               ms_InterCoord->SetDetectedAngleTheta(AngTheta) ;
            }
            Ang_Theta_itr++;
         }

         // Angle Phi
         Ang_Phi_itr = AngPhiHitMap->GetMap()->begin();
         for (unsigned int h = 0; h < AngPhiHitMap->entries(); h++) {
            G4int AngPhiTrackID =   Ang_Phi_itr->first - N    ;
            G4double AngPhi     = *(Ang_Phi_itr->second)      ;
            if (AngPhiTrackID == NTrackID) {
               ms_InterCoord->SetDetectedAnglePhi(AngPhi) ;
            }
            Ang_Phi_itr++;
         }

         // Second Stage
         SecondStageEnergy_itr = SecondStageEnergyHitMap->GetMap()->begin()  ;
         for (unsigned int h = 0 ; h < SecondStageEnergyHitMap->entries() ; h++) {
            G4int SecondStageEnergyTrackID  =   SecondStageEnergy_itr->first - N;
            G4double SecondStageEnergy      = *(SecondStageEnergy_itr->second);

            if (SecondStageEnergyTrackID == NTrackID) {
               ms_Event->SetGPDTrkSecondStageEEnergy(RandGauss::shoot(SecondStageEnergy, ResoSecondStage));
               ms_Event->SetGPDTrkSecondStageEPadNbr(1);
               ms_Event->SetGPDTrkSecondStageTPadNbr(1);
               ms_Event->SetGPDTrkSecondStageTTime(1);
               ms_Event->SetGPDTrkSecondStageTDetectorNbr(m_index["Trapezoid"] + N);
               ms_Event->SetGPDTrkSecondStageEDetectorNbr(m_index["Trapezoid"] + N);
            }
            SecondStageEnergy_itr++;
         }

         // Third Stage
         ThirdStageEnergy_itr = ThirdStageEnergyHitMap->GetMap()->begin()  ;
         for (unsigned int h = 0 ; h < ThirdStageEnergyHitMap->entries() ; h++) {
            G4int ThirdStageEnergyTrackID  =   ThirdStageEnergy_itr->first - N;
            G4double ThirdStageEnergy      = *(ThirdStageEnergy_itr->second);

            if (ThirdStageEnergyTrackID == NTrackID) {
               ms_Event->SetGPDTrkThirdStageEEnergy(RandGauss::shoot(ThirdStageEnergy, ResoThirdStage));
               ms_Event->SetGPDTrkThirdStageEPadNbr(1);
               ms_Event->SetGPDTrkThirdStageTPadNbr(1);
               ms_Event->SetGPDTrkThirdStageTTime(1);
               ms_Event->SetGPDTrkThirdStageTDetectorNbr(m_index["Trapezoid"] + N);
               ms_Event->SetGPDTrkThirdStageEDetectorNbr(m_index["Trapezoid"] + N);
            }
            ThirdStageEnergy_itr++;
         }

      }
      DetectorNumber_itr++;
   }

   // clear map for next event
   DetectorNumberHitMap ->clear();
   EnergyHitMap   ->clear();
   TimeHitMap     ->clear();
   XHitMap        ->clear();
   YHitMap        ->clear();
   PosXHitMap     ->clear();
   PosYHitMap     ->clear();
   PosZHitMap     ->clear();
   AngThetaHitMap ->clear();
   AngPhiHitMap   ->clear();
   SecondStageEnergyHitMap ->clear();
   ThirdStageEnergyHitMap ->clear();
}



void GaspardTrackerTrapezoid::InitializeScorers()
{
   bool already_exist = false;
   m_FirstStageScorer = NPS::VDetector::CheckScorer("FirstStageScorerGPDTrapezoid", already_exist);
   m_SecondStageScorer = NPS::VDetector::CheckScorer("SecondStageScorerGPDTrapezoid",already_exist);
   m_ThirdStageScorer = NPS::VDetector::CheckScorer("ThirdStageScorerGPDTrapezoid",already_exist);
   if(already_exist) return;



   // First stage Associate Scorer
   G4VPrimitiveScorer* DetNbr                           = new OBSOLETEGENERALSCORERS::PSDetectorNumber("DetectorNumber", "GPDTrapezoid", 0);
   G4VPrimitiveScorer* TOF                              = new OBSOLETEGENERALSCORERS::PSTOF("StripTime","GPDTrapezoid", 0);
   G4VPrimitiveScorer* InteractionCoordinatesX          = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesX("InterCoordX","GPDTrapezoid", 0);
   G4VPrimitiveScorer* InteractionCoordinatesY          = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesY("InterCoordY","GPDTrapezoid", 0);
   G4VPrimitiveScorer* InteractionCoordinatesZ          = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesZ("InterCoordZ","GPDTrapezoid", 0);
   G4VPrimitiveScorer* InteractionCoordinatesAngleTheta = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesAngleTheta("InterCoordAngTheta","GPDTrapezoid", 0);
   G4VPrimitiveScorer* InteractionCoordinatesAnglePhi   = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesAnglePhi("InterCoordAngPhi","GPDTrapezoid", 0);
   G4VPrimitiveScorer* Energy                           = new GPDScorerFirstStageEnergy("StripEnergy", "GPDTrapezoid", 0);
   G4VPrimitiveScorer* StripPositionX                   = new GPDScorerFirstStageFrontStripTrapezoid("StripNumberX", 0, NumberOfStripsX);
   G4VPrimitiveScorer* StripPositionY                   = new GPDScorerFirstStageBackStripTrapezoid("StripNumberY",  0, NumberOfStripsY);

   //and register it to the multifunctionnal detector
   m_FirstStageScorer->RegisterPrimitive(DetNbr);
   m_FirstStageScorer->RegisterPrimitive(Energy);
   m_FirstStageScorer->RegisterPrimitive(TOF);
   m_FirstStageScorer->RegisterPrimitive(StripPositionX);
   m_FirstStageScorer->RegisterPrimitive(StripPositionY);
   m_FirstStageScorer->RegisterPrimitive(InteractionCoordinatesX);
   m_FirstStageScorer->RegisterPrimitive(InteractionCoordinatesY);
   m_FirstStageScorer->RegisterPrimitive(InteractionCoordinatesZ);
   m_FirstStageScorer->RegisterPrimitive(InteractionCoordinatesAngleTheta);
   m_FirstStageScorer->RegisterPrimitive(InteractionCoordinatesAnglePhi);

   // Second stage Associate Scorer
   G4VPrimitiveScorer* SecondStageEnergy = new GPDScorerSecondStageEnergy("SecondStageEnergy", "GPDTrapezoid", 0);
   m_SecondStageScorer->RegisterPrimitive(SecondStageEnergy);

   //  Third stage Associate Scorer 
   G4VPrimitiveScorer* ThirdStageEnergy = new GPDScorerThirdStageEnergy("ThirdStageEnergy", "GPDTrapezoid", 0);
   m_ThirdStageScorer->RegisterPrimitive(ThirdStageEnergy);

   //  Add All Scorer to the Global Scorer Manager
   G4SDManager::GetSDMpointer()->AddNewDetector(m_FirstStageScorer);
   G4SDManager::GetSDMpointer()->AddNewDetector(m_SecondStageScorer);
   G4SDManager::GetSDMpointer()->AddNewDetector(m_ThirdStageScorer);
}
