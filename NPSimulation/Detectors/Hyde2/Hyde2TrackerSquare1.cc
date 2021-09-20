/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 10/06/09                                                 *
 * Last update    : 12/10/09                                                 *
 *---------------------------------------------------------------------------*
 * Decription: Define a module of square shape for the Hyde2 tracker       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + 07/09/09: Fix bug for placing module with (r,theta,phi) method.      *
 *                (N. de Sereville)                                          *
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
#include "G4Trd.hh"
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
#include "Hyde2TrackerSquare1.hh"
#include "ObsoleteGeneralScorers.hh"
#include "Hyde2Scorers.hh"
#include "RootOutput.h"
#include "NPSVDetector.hh"
// CLHEP
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace HYD2SQUARE1;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Hyde2TrackerSquare1::Hyde2TrackerSquare1()
{
   ms_InterCoord = 0;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Hyde2TrackerSquare1::~Hyde2TrackerSquare1()
{
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hyde2TrackerSquare1::AddModule(G4ThreeVector X1_Y1     ,
                                     G4ThreeVector X128_Y1   ,
                                     G4ThreeVector X1_Y128   ,
                                     G4ThreeVector X128_Y128 ,
                                     bool wFirstStage        ,
                                     bool wSecondStage       ,
                                     bool wThirdStage        ,
                                     bool wFourthStage       ,
                                     bool wFifthStage        ,
                                     bool wSixthStage)
{
   m_DefinitionType.push_back(true) ;

   m_X1_Y1.push_back(X1_Y1)               ;
   m_X128_Y1.push_back(X128_Y1)           ;
   m_X1_Y128.push_back(X1_Y128)           ;
   m_X128_Y128.push_back(X128_Y128)       ;
   m_wFirstStage.push_back(wFirstStage)   ;
   m_wSecondStage.push_back(wSecondStage) ;
   m_wThirdStage.push_back(wThirdStage)   ;
   m_wFourthStage.push_back(wFourthStage) ;
   m_wFifthStage.push_back(wFifthStage)   ;
   m_wSixthStage.push_back(wSixthStage)   ;

   m_R.push_back(0)      ;
   m_Theta.push_back(0)  ;
   m_Phi.push_back(0)    ;
   m_beta_u.push_back(0) ;
   m_beta_v.push_back(0) ;
   m_beta_w.push_back(0) ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hyde2TrackerSquare1::AddModule(G4double R        ,
                                     G4double Theta    ,
                                     G4double Phi      ,
                                     G4double beta_u   ,
                                     G4double beta_v   ,
                                     G4double beta_w   ,
                                     bool wFirstStage  ,
                                     bool wSecondStage ,
                                     bool wThirdStage  ,
                                     bool wFourthStage       ,
                                     bool wFifthStage        ,
                                     bool wSixthStage)
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
   m_wFourthStage.push_back(wFourthStage) ;
   m_wFifthStage.push_back(wFifthStage)   ;
   m_wSixthStage.push_back(wSixthStage)   ;

   m_X1_Y1.push_back(empty)     ;
   m_X128_Y1.push_back(empty)   ;
   m_X1_Y128.push_back(empty)   ;
   m_X128_Y128.push_back(empty) ;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hyde2TrackerSquare1::VolumeMaker(G4int TelescopeNumber,
                                       G4ThreeVector MMpos,
                                       G4RotationMatrix* MMrot,
                                       bool wFirstStage,
                                       bool wSecondStage,
                                       bool wThirdStage,
                                       bool wFourthStage,
                                       bool wFifthStage        ,
                                       bool wSixthStage        ,
                                       G4LogicalVolume* world)
{
   G4double NbrTelescopes = TelescopeNumber  ;
   G4String DetectorNumber                   ;
   ostringstream Number                      ;
   Number << NbrTelescopes                   ;
   DetectorNumber = Number.str()             ;

   ////////////////////////////////////////////////////////////////
   ////////////// Starting Volume Definition //////////////////////
   ////////////////////////////////////////////////////////////////
   G4String Name = "HYD2Square1" + DetectorNumber;

   G4Box*           solidHYD2Square1 = new G4Box(Name, 0.5*FaceFront, 0.5*FaceFront, 0.5*Length);
   G4LogicalVolume* logicHYD2Square1 = new G4LogicalVolume(solidHYD2Square1, m_MaterialVacuum, Name, 0, 0, 0);

   new G4PVPlacement(G4Transform3D(*MMrot, MMpos), logicHYD2Square1, Name, world, false, 0);

   logicHYD2Square1->SetVisAttributes(G4VisAttributes::Invisible);
   if (m_non_sensitive_part_visiualisation) logicHYD2Square1->SetVisAttributes(G4VisAttributes(G4Colour(0.90, 0.90, 0.90)));

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
   //////////////// First Stage Construction //////////////////////
   ////////////////////////////////////////////////////////////////
   if (wFirstStage) {
      // Silicon detector itself
      G4ThreeVector  positionFirstStage = G4ThreeVector(0, 0, FirstStage_PosZ);

      G4Box*           solidFirstStage = new G4Box("solidFirstStage", 0.5*FirstStageFace, 0.5*FirstStageFace, 0.5*FirstStageThickness);
      G4LogicalVolume* logicFirstStage = new G4LogicalVolume(solidFirstStage, m_MaterialSilicon, "logicFirstStage", 0, 0, 0);

      new G4PVPlacement(0,
                                    positionFirstStage,
                                    logicFirstStage,
                                    Name + "_FirstStage",
                                    logicHYD2Square1,
                                    false,
                                    0);

      // Set First Stage sensible
      logicFirstStage->SetSensitiveDetector(m_FirstStageScorer);

      ///Visualisation of FirstStage Strip
      G4VisAttributes* FirstStageVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));   // blue
      logicFirstStage->SetVisAttributes(FirstStageVisAtt);
   }

   ////////////////////////////////////////////////////////////////
   //////////////////// Second Stage  Construction ////////////////
   ////////////////////////////////////////////////////////////////
   if (wSecondStage) {
      // Second stage silicon detector
      G4ThreeVector  positionSecondStage = G4ThreeVector(0, 0, SecondStage_PosZ);

      G4Box*           solidSecondStage = new G4Box("solidSecondStage", 0.5*SecondStageFace, 0.5*SecondStageFace, 0.5*SecondStageThickness);
      G4LogicalVolume* logicSecondStage = new G4LogicalVolume(solidSecondStage, m_MaterialSilicon, "logicSecondStage", 0, 0, 0);

      new G4PVPlacement(0,
                                    positionSecondStage,
                                    logicSecondStage,
                                    Name + "_SecondStage",
                                    logicHYD2Square1,
                                    false,
                                    0);

      // Set Second Stage sensible
      logicSecondStage->SetSensitiveDetector(m_SecondStageScorer);

      ///Visualisation of SecondStage Strip
      G4VisAttributes* SecondStageVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
      logicSecondStage->SetVisAttributes(SecondStageVisAtt)                        ;
   }

   ////////////////////////////////////////////////////////////////
   ///////////////// Third Stage Construction /////////////////////
   ////////////////////////////////////////////////////////////////
   if (wThirdStage) {
      // Third stage silicon detector
      G4ThreeVector  positionThirdStage = G4ThreeVector(0, 0, ThirdStage_PosZ);

      G4Box*           solidThirdStage = new G4Box("solidThirdStage", 0.5*ThirdStageFace, 0.5*ThirdStageFace, 0.5*ThirdStageThickness);
      G4LogicalVolume* logicThirdStage = new G4LogicalVolume(solidThirdStage, m_MaterialSilicon, "logicThirdStage", 0, 0, 0);

      new G4PVPlacement(0,
                                    positionThirdStage,
                                    logicThirdStage,
                                    Name + "_ThirdStage",
                                    logicHYD2Square1,
                                    false,
                                    0);

      // Set Third Stage sensible
      logicThirdStage->SetSensitiveDetector(m_ThirdStageScorer);

      ///Visualisation of Third Stage
      G4VisAttributes* ThirdStageVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.0));   // green
      logicThirdStage->SetVisAttributes(ThirdStageVisAtt);
   }

   ////////////////////////////////////////////////////////////////
   ///////////////// Fourth Stage Construction/////////////////////
   ////////////////////////////////////////////////////////////////
   if (wFourthStage) {
      // Fourth stage silicon detector
      G4ThreeVector  positionFourthStage = G4ThreeVector(0, 0, FourthStage_PosZ);

      G4Box*           solidFourthStage = new G4Box("solidFourthStage", 0.5*FourthStageFace, 0.5*FourthStageFace, 0.5*FourthStageThickness);
      G4LogicalVolume* logicFourthStage = new G4LogicalVolume(solidFourthStage, m_MaterialSilicon, "logicFourthStage", 0, 0, 0);

      new G4PVPlacement(0,
                                    positionFourthStage,
                                    logicFourthStage,
                                    Name + "_FourthStage",
                                    logicHYD2Square1,
                                    false,
                                    0);

      // Set Fourth Stage sensible
      logicFourthStage->SetSensitiveDetector(m_FourthStageScorer);

      ///Visualisation of Fourth Stage
      G4VisAttributes* FourthStageVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.0));   // green
      logicFourthStage->SetVisAttributes(FourthStageVisAtt);
   }

   ////////////////////////////////////////////////////////////////
   ///////////////// Fifth Stage Construction/////////////////////
   ////////////////////////////////////////////////////////////////
   if (wFifthStage) {
      // Fifth stage silicon detector
      G4ThreeVector  positionFifthStage = G4ThreeVector(0, 0, FifthStage_PosZ);

      G4Box*           solidFifthStage = new G4Box("solidFifthStage", 0.5*FifthStageFace, 0.5*FifthStageFace, 0.5*FifthStageThickness);
      G4LogicalVolume* logicFifthStage = new G4LogicalVolume(solidFifthStage, m_MaterialSilicon, "logicFifthStage", 0, 0, 0);

      new G4PVPlacement(0,
                                   positionFifthStage,
                                    logicFifthStage,
                                    Name + "_FifthStage",
                                    logicHYD2Square1,
                                    false,
                                    0);

      // Set Fifth Stage sensible
      logicFifthStage->SetSensitiveDetector(m_FifthStageScorer);

      ///Visualisation of Fifth Stage
      G4VisAttributes* FifthStageVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.0));   // green
      logicFifthStage->SetVisAttributes(FifthStageVisAtt);
   }

   ////////////////////////////////////////////////////////////////
   ///////////////// Sixth Stage Construction/////////////////////
   ////////////////////////////////////////////////////////////////
   if (wSixthStage) {
      // Sixth stage silicon detector
      G4ThreeVector  positionSixthStage = G4ThreeVector(0, 0, SixthStage_PosZ);

      G4Box*           solidSixthStage = new G4Box("solidSixthStage", 0.5*SixthStageFace, 0.5*SixthStageFace, 0.5*SixthStageThickness);
      G4LogicalVolume* logicSixthStage = new G4LogicalVolume(solidSixthStage, m_MaterialSilicon, "logicSixthStage", 0, 0, 0);

      new G4PVPlacement(0,
                                    positionSixthStage,
                                    logicSixthStage,
                                    Name + "_SixthStage",
                                    logicHYD2Square1,
                                    false,
                                    0);

      // Set Sixth Stage sensible
      logicSixthStage->SetSensitiveDetector(m_SixthStageScorer);

      ///Visualisation of Sixth Stage
      G4VisAttributes* SixthStageVisAtt = new G4VisAttributes(G4Colour(0.0, 0.9, 0.0));   // green
      logicSixthStage->SetVisAttributes(SixthStageVisAtt);
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Hyde2TrackerSquare1::ReadConfiguration(string Path)
{
   ifstream ConfigFile           ;
   ConfigFile.open(Path.c_str()) ;
   string LineBuffer          ;
   string DataBuffer          ;

   // A:X1_Y1     --> X:1    Y:1
   // B:X128_Y1   --> X:128  Y:1
   // C:X1_Y128   --> X:1    Y:128
   // D:X128_Y128    --> X:128  Y:128

   G4double Ax , Bx , Cx , Dx , Ay , By , Cy , Dy , Az , Bz , Cz , Dz          ;
   G4ThreeVector A , B , C , D                                                 ;
   G4double Theta = 0 , Phi = 0 , R = 0 , beta_u = 0 , beta_v = 0 , beta_w = 0 ;
   int FIRSTSTAGE = 0 , SECONDSTAGE = 0 , THIRDSTAGE = 0 , FOURTHSTAGE = 0 , FIFTHSTAGE = 0, SIXTHSTAGE = 0    ;

   bool ReadingStatus = false ;

   bool check_A = false ;
   bool check_C = false ;
   bool check_B = false ;
   bool check_D = false ;

   bool check_Theta = false ;
   bool check_Phi   = false ;
   bool check_R     = false ;
//   bool check_beta  = false ;
   
   bool check_FirstStage = false ;
   bool check_SecondStage = false ;
   bool check_ThirdStage = false ;
   bool check_FourthStage = false ;
   bool check_FifthStage = false ;
   bool check_SixthStage = false ;
   bool checkVis = false ;

   while (!ConfigFile.eof()) {
      getline(ConfigFile, LineBuffer);
      if (LineBuffer.compare(0, 11, "HYD2Square1") == 0) {
         G4cout << "///" << G4endl           ;
         G4cout << "Square1 element found: " << G4endl   ;
         ReadingStatus = true ;
         }
         
   while(ReadingStatus){      

         ConfigFile >> DataBuffer;
         //   Comment Line 
      if (DataBuffer.compare(0, 1, "%") == 0) {/*do nothing */;}
      
         // Position method
         else if (DataBuffer.compare(0, 6, "X1_Y1=") == 0) {
            check_A = true;
            ConfigFile >> DataBuffer ;
            Ax = atof(DataBuffer.c_str()) ;
            Ax = Ax * mm ;
            ConfigFile >> DataBuffer ;
            Ay = atof(DataBuffer.c_str()) ;
            Ay = Ay * mm ;
            ConfigFile >> DataBuffer ;
            Az = atof(DataBuffer.c_str()) ;
            Az = Az * mm ;

            A = G4ThreeVector(Ax, Ay, Az);
            G4cout << "X1 Y1 corner position : " << A << G4endl;
         }
        
         else if (DataBuffer.compare(0, 8, "X128_Y1=") == 0) {
            check_B = true;
            ConfigFile >> DataBuffer ;
            Bx = atof(DataBuffer.c_str()) ;
            Bx = Bx * mm ;
            ConfigFile >> DataBuffer ;
            By = atof(DataBuffer.c_str()) ;
            By = By * mm ;
            ConfigFile >> DataBuffer ;
            Bz = atof(DataBuffer.c_str()) ;
            Bz = Bz * mm ;

            B = G4ThreeVector(Bx, By, Bz);
            G4cout << "X128 Y1 corner position : " << B << G4endl;
         }
         
         else if (DataBuffer.compare(0, 8, "X1_Y128=") == 0) {
            check_C = true;
            ConfigFile >> DataBuffer ;
            Cx = atof(DataBuffer.c_str()) ;
            Cx = Cx * mm ;
            ConfigFile >> DataBuffer ;
            Cy = atof(DataBuffer.c_str()) ;
            Cy = Cy * mm ;
            ConfigFile >> DataBuffer ;
            Cz = atof(DataBuffer.c_str()) ;
            Cz = Cz * mm ;

            C = G4ThreeVector(Cx, Cy, Cz);
            G4cout << "X1 Y128 corner position : " << C << G4endl;
         }
        
         else if (DataBuffer.compare(0, 10, "X128_Y128=") == 0) {
            check_D = true;
            ConfigFile >> DataBuffer ;
            Dx = atof(DataBuffer.c_str()) ;
            Dx = Dx * mm ;
            ConfigFile >> DataBuffer ;
            Dy = atof(DataBuffer.c_str()) ;
            Dy = Dy * mm ;
            ConfigFile >> DataBuffer ;
            Dz = atof(DataBuffer.c_str()) ;
            Dz = Dz * mm ;

            D = G4ThreeVector(Dx, Dy, Dz);
            G4cout << "X128 Y128 corner position : " << D << G4endl;
         }
         

       // Angle method
         else if (DataBuffer.compare(0, 6, "THETA=") == 0) {
            check_Theta = true;
            ConfigFile >> DataBuffer ;
            Theta = atof(DataBuffer.c_str()) ;
            Theta = Theta * deg;
            G4cout << "Theta:  " << Theta / deg << G4endl;
         }

         else if (DataBuffer.compare(0, 4, "PHI=") == 0) {
            check_Phi = true;
            ConfigFile >> DataBuffer ;
            Phi = atof(DataBuffer.c_str()) ;
            Phi = Phi * deg;
            G4cout << "Phi:  " << Phi / deg << G4endl;
         }

         else if (DataBuffer.compare(0, 2, "R=") == 0) {
            check_R = true;
            ConfigFile >> DataBuffer ;
            R = atof(DataBuffer.c_str()) ;
            R = R * mm;
            G4cout << "R:  " << R / mm << G4endl;
         }

         else if (DataBuffer.compare(0, 5, "BETA=") == 0) {
//            check_beta = true;
            ConfigFile >> DataBuffer ;
            beta_u = atof(DataBuffer.c_str()) ;
            beta_u = beta_u * deg   ;
            ConfigFile >> DataBuffer ;
            beta_v = atof(DataBuffer.c_str()) ;
            beta_v = beta_v * deg   ;
            ConfigFile >> DataBuffer ;
            beta_w = atof(DataBuffer.c_str()) ;
            beta_w = beta_w * deg   ;
            G4cout << "Beta:  " << beta_u / deg << " " << beta_v / deg << " " << beta_w / deg << G4endl  ;
         }

         else if (DataBuffer.compare(0, 11, "FIRSTSTAGE=") == 0) {
            check_FirstStage = true ;
            ConfigFile >> DataBuffer;
            FIRSTSTAGE = atof(DataBuffer.c_str()) ;
         }

         else if (DataBuffer.compare(0, 12, "SECONDSTAGE=") == 0) {
            check_SecondStage = true ;
            ConfigFile >> DataBuffer;
            SECONDSTAGE = atof(DataBuffer.c_str()) ;
         }

         else if (DataBuffer.compare(0, 11, "THIRDSTAGE=") == 0) {
            check_ThirdStage = true ;
            ConfigFile >> DataBuffer;
            THIRDSTAGE = atof(DataBuffer.c_str()) ;
         }

         else if (DataBuffer.compare(0, 12, "FOURTHSTAGE=") == 0) {
            check_FourthStage = true ;
            ConfigFile >> DataBuffer;
            FOURTHSTAGE = atof(DataBuffer.c_str()) ;
         }

         else if (DataBuffer.compare(0, 11, "FIFTHSTAGE=") == 0) {
            check_FifthStage = true ;
            ConfigFile >> DataBuffer;
            FIFTHSTAGE = atof(DataBuffer.c_str()) ;
         }

         else if (DataBuffer.compare(0, 11, "SIXTHSTAGE=") == 0) {
            check_SixthStage = true ;
            ConfigFile >> DataBuffer;
            SIXTHSTAGE = atof(DataBuffer.c_str()) ;
         }

         else if (DataBuffer.compare(0, 4, "VIS=") == 0) {
            checkVis = true ;
            ConfigFile >> DataBuffer;
            if (DataBuffer.compare(0, 3, "all") == 0) m_non_sensitive_part_visiualisation = true;
         }
         
         else G4cout << "WARNING: Wrong Token, Hyde2TrackerSquare1: Square1 Element not added" << G4endl;

         //Add The previously define telescope
         //With position method
         if ((check_A && check_B && check_C && check_D && check_FirstStage && check_SecondStage && check_ThirdStage && check_FourthStage && check_FifthStage && check_SixthStage && checkVis) && !(check_Theta && check_Phi && check_R)) {
         
            ReadingStatus = false ;
          check_A = false ;
          check_C = false ;
          check_B = false ;
          check_D = false ;
          check_FirstStage = false ;
          check_SecondStage = false ;
          check_ThirdStage = false ;
          check_FourthStage = false ;
          check_FifthStage = false ;
          check_SixthStage = false ;
          checkVis = false ;
         
            AddModule(A                ,
                      B                ,
                      C                ,
                      D                ,
                      FIRSTSTAGE  == 1 ,
                      SECONDSTAGE == 1 ,
                      THIRDSTAGE  == 1 ,
                      FOURTHSTAGE == 1 ,
                      FIFTHSTAGE  == 1 ,
                      SIXTHSTAGE  == 1);
         }

         //with angle method
        if ((check_Theta && check_Phi && check_R && check_FirstStage && check_SecondStage && check_ThirdStage && check_FourthStage && check_FifthStage && check_SixthStage && checkVis) && !(check_A && check_B && check_C && check_D)) {
            ReadingStatus = false ;
             check_Theta = false ;
             check_Phi   = false ;
             check_R     = false ;
//             check_beta  = false ;
           check_FirstStage = false ;
          check_SecondStage = false ;
           check_ThirdStage = false ;
           check_FourthStage = false ;
           check_FifthStage = false ;
           check_SixthStage = false ;
           checkVis = false ;
           
            AddModule(R                ,
                      Theta            ,
                      Phi              ,
                      beta_u           ,
                      beta_v           ,
                      beta_w           ,
                      FIRSTSTAGE  == 1 ,
                      SECONDSTAGE == 1 ,
                      THIRDSTAGE  == 1 ,
                      FOURTHSTAGE == 1 ,
                      FIFTHSTAGE  == 1 ,
                      SIXTHSTAGE  == 1);
         }

         
      }
   }
}

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Hyde2TrackerSquare1::ConstructDetector(G4LogicalVolume* world)
{
   G4RotationMatrix* MMrot    = NULL                   ;
   G4ThreeVector     MMpos    = G4ThreeVector(0, 0, 0) ;
   G4ThreeVector     MMu      = G4ThreeVector(0, 0, 0) ;
   G4ThreeVector     MMv      = G4ThreeVector(0, 0, 0) ;
   G4ThreeVector     MMw      = G4ThreeVector(0, 0, 0) ;
   G4ThreeVector     MMCenter = G4ThreeVector(0, 0, 0) ;
   bool FirstStage  = true;
   bool SecondStage = true;
   bool ThirdStage  = true;
   bool FourthStage  = true;
   bool FifthStage  = true;
   bool SixthStage  = true;

   G4int NumberOfTelescope = m_DefinitionType.size() ;

   for (G4int i = 0; i < NumberOfTelescope; i++) {
      // By Point
      if (m_DefinitionType[i]) {
         // (u,v,w) unitary vector associated to telescope referencial
         // (u,v) // to silicon plan
         // w perpendicular to (u,v) plan and pointing ThirdStage
         MMu = m_X128_Y1[i] - m_X1_Y1[i];
         MMu = MMu.unit();

         MMv = m_X1_Y128[i] - m_X1_Y1[i];
         MMv = MMv.unit();

         MMw = MMu.cross(MMv);
         MMw = MMw.unit();

         MMCenter = (m_X1_Y1[i] + m_X1_Y128[i] + m_X128_Y1[i] + m_X128_Y128[i]) / 4;

         // Passage Matrix from Lab Referential to Telescope Referential
         MMrot = new G4RotationMatrix(MMu, MMv, MMw);
         // translation to place Telescope
         MMpos = MMw * Length * 0.5 + MMCenter;
      }

      // By Angle
      else {
         G4double Theta = m_Theta[i] ;
         G4double Phi   = m_Phi[i]   ;
         
         // (u,v,w) unitary vector associated to telescope referencial
         // (u,v) // to silicon plan
         // w perpendicular to (u,v) plan and pointing ThirdStage
         // Phi is angle between X axis and projection in (X,Y) plan
         // Theta is angle between  position vector and z axis
         G4double wX = m_R[i] * sin(Theta / rad) * cos(Phi / rad);
         G4double wY = m_R[i] * sin(Theta / rad) * sin(Phi / rad);
         G4double wZ = m_R[i] * cos(Theta / rad);
         MMw = G4ThreeVector(wX, wY, wZ);

         // vector corresponding to the center of the module
         MMCenter = MMw;

         // vector parallel to one axis of silicon plane
         G4double ii = cos(Theta / rad) * cos(Phi / rad);
         G4double jj = cos(Theta / rad) * sin(Phi / rad);
         G4double kk = -sin(Theta / rad);
         G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

         MMw = MMw.unit();
         MMu = MMw.cross(Y);
         MMv = MMw.cross(MMu);
         MMv = MMv.unit();
         MMu = MMu.unit();

         // Passage Matrix from Lab Referential to Telescope Referential
         // MUST2
         MMrot = new G4RotationMatrix(MMu, MMv, MMw);
         // Telescope is rotate of Beta angle around MMv axis.
         MMrot->rotate(m_beta_u[i], MMu);
         MMrot->rotate(m_beta_v[i], MMv);
         MMrot->rotate(m_beta_w[i], MMw);
         // translation to place Telescope
         MMpos = MMw * Length * 0.5 + MMCenter;
      }

      FirstStage  = m_wFirstStage[i]  ;
      SecondStage = m_wSecondStage[i] ;
      ThirdStage  = m_wThirdStage[i]  ;
      FourthStage  = m_wFourthStage[i]  ;
      FifthStage  = m_wFifthStage[i]  ;
      SixthStage  = m_wSixthStage[i]  ;

      VolumeMaker(i + 1, MMpos, MMrot, FirstStage, SecondStage, ThirdStage , FourthStage, FifthStage, SixthStage, world);
   }

   delete MMrot ;
}



// Connect the Hyde2TrackingData class to the output TTree
// of the simulation
void Hyde2TrackerSquare1::InitializeRootOutput()
{
}



// Set the TinteractionCoordinates object from NPS::VDetector to the present class
void Hyde2TrackerSquare1::SetInterCoordPointer(TInteractionCoordinates* interCoord)
{
   ms_InterCoord = interCoord;
}



// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Hyde2TrackerSquare1::ReadSensitive(const G4Event* event)
{
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Used to Read Event Map of detector //////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
   // First Stage
   std::map<G4int, G4int*>::iterator    DetectorNumber_itr;
   std::map<G4int, G4double*>::iterator Energy_itr;
   std::map<G4int, G4double*>::iterator Time_itr;
   std::map<G4int, G4double*>::iterator X_itr;
   std::map<G4int, G4double*>::iterator Y_itr;
   std::map<G4int, G4double*>::iterator Pos_X_itr;
   std::map<G4int, G4double*>::iterator Pos_Y_itr;
   std::map<G4int, G4double*>::iterator Pos_Z_itr;
   std::map<G4int, G4double*>::iterator Ang_Theta_itr;
   std::map<G4int, G4double*>::iterator Ang_Phi_itr;

   G4THitsMap<G4int>*    DetectorNumberHitMap;
   G4THitsMap<G4double>* EnergyHitMap;
   G4THitsMap<G4double>* TimeHitMap;
   G4THitsMap<G4double>* XHitMap;
   G4THitsMap<G4double>* YHitMap;
   G4THitsMap<G4double>* PosXHitMap;
   G4THitsMap<G4double>* PosYHitMap;
   G4THitsMap<G4double>* PosZHitMap;
   G4THitsMap<G4double>* AngThetaHitMap;
   G4THitsMap<G4double>* AngPhiHitMap;

   // NULL pointer are given to avoid warning at compilation
   // Second Stage
   std::map<G4int, G4double*>::iterator SecondStageEnergy_itr ;
   G4THitsMap<G4double>* SecondStageEnergyHitMap = NULL      ;
   // Third Stage
   std::map<G4int, G4double*>::iterator ThirdStageEnergy_itr  ;
   G4THitsMap<G4double>* ThirdStageEnergyHitMap = NULL    ;
   // Fourth Stage
   std::map<G4int, G4double*>::iterator FourthStageEnergy_itr  ;
   G4THitsMap<G4double>* FourthStageEnergyHitMap = NULL    ;
   // Fifth Stage
   std::map<G4int, G4double*>::iterator FifthStageEnergy_itr  ;
   G4THitsMap<G4double>* FifthStageEnergyHitMap = NULL    ;
   // Sixth Stage
   std::map<G4int, G4double*>::iterator SixthStageEnergy_itr  ;
   G4THitsMap<G4double>* SixthStageEnergyHitMap = NULL    ;



   // Read the Scorer associate to the Silicon Strip
   //Detector Number
   G4int StripDetCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/DetectorNumber")    ;
   DetectorNumberHitMap = (G4THitsMap<G4int>*)(event->GetHCofThisEvent()->GetHC(StripDetCollectionID))         ;
   DetectorNumber_itr =  DetectorNumberHitMap->GetMap()->begin()                                               ;

   //Energy
   G4int StripEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/StripEnergy")   ;
   EnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(StripEnergyCollectionID))                    ;
   Energy_itr = EnergyHitMap->GetMap()->begin()                                                          ;

   //Time of Flight
   G4int StripTimeCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/StripTime")    ;
   TimeHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(StripTimeCollectionID))                        ;
   Time_itr = TimeHitMap->GetMap()->begin()                                                              ;

   //Strip Number X
   G4int StripXCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/StripNumberX")    ;
   XHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(StripXCollectionID))                              ;
   X_itr = XHitMap->GetMap()->begin()                                                                    ;

   //Strip Number Y
   G4int StripYCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/StripNumberY")    ;
   YHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(StripYCollectionID))                              ;
   Y_itr = YHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate X
   G4int InterCoordXCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/InterCoordX")    ;
   PosXHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordXCollectionID))                              ;
   Pos_X_itr = PosXHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Y
   G4int InterCoordYCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/InterCoordY")    ;
   PosYHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordYCollectionID))                              ;
   Pos_Y_itr = PosYHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Z
   G4int InterCoordZCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/InterCoordZ")    ;
   PosZHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordZCollectionID))                              ;
   Pos_Z_itr = PosXHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Angle Theta
   G4int InterCoordAngThetaCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/InterCoordAngTheta")    ;
   AngThetaHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordAngThetaCollectionID))                              ;
   Ang_Theta_itr = AngThetaHitMap->GetMap()->begin()                                                                    ;

   //Interaction Coordinate Angle Phi
   G4int InterCoordAngPhiCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FirstStageScorerHYD2Square1/InterCoordAngPhi")    ;
   AngPhiHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(InterCoordAngPhiCollectionID))                              ;
   Ang_Phi_itr = AngPhiHitMap->GetMap()->begin()                                                                    ;


   // Read the Scorer associate to the SecondStage
   //Energy
   G4int SecondStageEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("SecondStageScorerHYD2Square1/SecondStageEnergy")   ;
   SecondStageEnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(SecondStageEnergyCollectionID))                 ;
   SecondStageEnergy_itr = SecondStageEnergyHitMap->GetMap()->begin()                                                     ;


   // Read the Scorer associate to the ThirdStage
   //Energy
   G4int ThirdStageEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("ThirdStageScorerHYD2Square1/ThirdStageEnergy");
   ThirdStageEnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(ThirdStageEnergyCollectionID));
   ThirdStageEnergy_itr = ThirdStageEnergyHitMap->GetMap()->begin();

   // Read the Scorer associate to the FourthStage
   //Energy
   G4int FourthStageEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FourthStageScorerHYD2Square1/FourthStageEnergy");
   FourthStageEnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(FourthStageEnergyCollectionID));
   FourthStageEnergy_itr = FourthStageEnergyHitMap->GetMap()->begin();

   // Read the Scorer associate to the FifthStage
   //Energy
   G4int FifthStageEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("FifthStageScorerHYD2Square1/FifthStageEnergy");
   FifthStageEnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(FifthStageEnergyCollectionID));
   FifthStageEnergy_itr = FifthStageEnergyHitMap->GetMap()->begin();

   // Read the Scorer associate to the SixthStage
   //Energy
   G4int SixthStageEnergyCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("SixthStageScorerHYD2Square1/SixthStageEnergy");
   SixthStageEnergyHitMap = (G4THitsMap<G4double>*)(event->GetHCofThisEvent()->GetHC(SixthStageEnergyCollectionID));
   SixthStageEnergy_itr = SixthStageEnergyHitMap->GetMap()->begin();

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
         ms_Event->SetHYD2TrkFirstStageFrontEDetectorNbr(m_index["Square1"] + N);
         ms_Event->SetHYD2TrkFirstStageFrontTDetectorNbr(m_index["Square1"] + N);
         ms_Event->SetHYD2TrkFirstStageBackEDetectorNbr(m_index["Square1"] + N);
         ms_Event->SetHYD2TrkFirstStageBackTDetectorNbr(m_index["Square1"] + N);

         // Energy
         for (G4int ll = 0 ; ll < sizeE ; ll++) {
            G4int ETrackID  =   Energy_itr->first - N;
            G4double E     = *(Energy_itr->second);
            if (ETrackID == NTrackID) {
               ms_Event->SetHYD2TrkFirstStageFrontEEnergy(RandGauss::shoot(E, ResoFirstStage));
               ms_Event->SetHYD2TrkFirstStageBackEEnergy(RandGauss::shoot(E, ResoFirstStage));
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
               ms_Event->SetHYD2TrkFirstStageFrontTTime(RandGauss::shoot(T, ResoTimeHyd2)) ;
               ms_Event->SetHYD2TrkFirstStageBackTTime(RandGauss::shoot(T, ResoTimeHyd2)) ;
            }
            Time_itr++;
         }

            // X
            X_itr = XHitMap->GetMap()->begin();
            for (G4int h = 0 ; h < sizeX ; h++) {
               G4int XTrackID  =   X_itr->first - N;
               G4double X     = *(X_itr->second);
               if (XTrackID == NTrackID) {
                  ms_Event->SetHYD2TrkFirstStageFrontEStripNbr(X);
                  ms_Event->SetHYD2TrkFirstStageFrontTStripNbr(X);
               }

               X_itr++;
            }

            // Y
            Y_itr = YHitMap->GetMap()->begin()  ;
            for (G4int h = 0 ; h < sizeY ; h++) {
               G4int YTrackID  =   Y_itr->first - N;
               G4double Y     = *(Y_itr->second);
               if (YTrackID == NTrackID) {
                  ms_Event->SetHYD2TrkFirstStageBackEStripNbr(Y);
                  ms_Event->SetHYD2TrkFirstStageBackTStripNbr(Y);
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
               G4int PosYTrackID =   Pos_Y_itr->first - N    ;
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
               SecondStageEnergy_itr = SecondStageEnergyHitMap->GetMap()->begin() ;
               for (unsigned int h = 0 ; h < SecondStageEnergyHitMap->entries() ; h++) {
                  G4int SecondStageEnergyTrackID =   SecondStageEnergy_itr->first - N;
                  G4double SecondStageEnergy     = *(SecondStageEnergy_itr->second);

                  if (SecondStageEnergyTrackID == NTrackID) {
                     ms_Event->SetHYD2TrkSecondStageEEnergy(RandGauss::shoot(SecondStageEnergy, ResoSecondStage)) ;
                     ms_Event->SetHYD2TrkSecondStageEPadNbr(1);
                     ms_Event->SetHYD2TrkSecondStageTPadNbr(1);
                     ms_Event->SetHYD2TrkSecondStageTTime(1);
                     ms_Event->SetHYD2TrkSecondStageTDetectorNbr(m_index["Square1"] + N);
                     ms_Event->SetHYD2TrkSecondStageEDetectorNbr(m_index["Square1"] + N);
                  }

                  SecondStageEnergy_itr++;
               }

            // Third Stage
               ThirdStageEnergy_itr = ThirdStageEnergyHitMap->GetMap()->begin()  ;
               for (unsigned int h = 0 ; h < ThirdStageEnergyHitMap->entries() ; h++) {
                  G4int ThirdStageEnergyTrackID  =   ThirdStageEnergy_itr->first - N;
                  G4double ThirdStageEnergy      = *(ThirdStageEnergy_itr->second)    ;

                  if (ThirdStageEnergyTrackID == NTrackID) {
                     ms_Event->SetHYD2TrkThirdStageEEnergy(RandGauss::shoot(ThirdStageEnergy, ResoThirdStage));
                     ms_Event->SetHYD2TrkThirdStageEPadNbr(1);
                     ms_Event->SetHYD2TrkThirdStageTPadNbr(1);
                     ms_Event->SetHYD2TrkThirdStageTTime(1);
                     ms_Event->SetHYD2TrkThirdStageTDetectorNbr(m_index["Square1"] + N);
                     ms_Event->SetHYD2TrkThirdStageEDetectorNbr(m_index["Square1"] + N);
                  }

                  ThirdStageEnergy_itr++;
               }

            // Fourth Stage
               FourthStageEnergy_itr = FourthStageEnergyHitMap->GetMap()->begin()  ;
               for (unsigned int h = 0 ; h < FourthStageEnergyHitMap->entries() ; h++) {
                  G4int FourthStageEnergyTrackID  =   FourthStageEnergy_itr->first - N;
                  G4double FourthStageEnergy      = *(FourthStageEnergy_itr->second)    ;

                  if (FourthStageEnergyTrackID == NTrackID) {
                     ms_Event->SetHYD2TrkFourthStageEEnergy(RandGauss::shoot(FourthStageEnergy, ResoFourthStage));
                     ms_Event->SetHYD2TrkFourthStageEPadNbr(1);
                     ms_Event->SetHYD2TrkFourthStageTPadNbr(1);
                     ms_Event->SetHYD2TrkFourthStageTTime(1);
                     ms_Event->SetHYD2TrkFourthStageTDetectorNbr(m_index["Square1"] + N);
                     ms_Event->SetHYD2TrkFourthStageEDetectorNbr(m_index["Square1"] + N);
                  }

                  FourthStageEnergy_itr++;
               }

            // Fifth Stage
               FifthStageEnergy_itr = FifthStageEnergyHitMap->GetMap()->begin()  ;
               for (unsigned int h = 0 ; h < FifthStageEnergyHitMap->entries() ; h++) {
                  G4int FifthStageEnergyTrackID  =   FifthStageEnergy_itr->first - N;
                  G4double FifthStageEnergy      = *(FifthStageEnergy_itr->second)    ;

                  if (FifthStageEnergyTrackID == NTrackID) {
                     ms_Event->SetHYD2TrkFifthStageEEnergy(RandGauss::shoot(FifthStageEnergy, ResoFifthStage));
                     ms_Event->SetHYD2TrkFifthStageEPadNbr(1);
                     ms_Event->SetHYD2TrkFifthStageTPadNbr(1);
                     ms_Event->SetHYD2TrkFifthStageTTime(1);
                     ms_Event->SetHYD2TrkFifthStageTDetectorNbr(m_index["Square1"] + N);
                     ms_Event->SetHYD2TrkFifthStageEDetectorNbr(m_index["Square1"] + N);
                  }

                  FifthStageEnergy_itr++;
               }

            // Sixth Stage
               SixthStageEnergy_itr = SixthStageEnergyHitMap->GetMap()->begin()  ;
               for (unsigned int h = 0 ; h < SixthStageEnergyHitMap->entries() ; h++) {
                  G4int SixthStageEnergyTrackID  =   SixthStageEnergy_itr->first - N;
                  G4double SixthStageEnergy      = *(SixthStageEnergy_itr->second)    ;

                  if (SixthStageEnergyTrackID == NTrackID) {
                     ms_Event->SetHYD2TrkSixthStageEEnergy(RandGauss::shoot(SixthStageEnergy, ResoSixthStage));
                     ms_Event->SetHYD2TrkSixthStageEPadNbr(1);
                     ms_Event->SetHYD2TrkSixthStageTPadNbr(1);
                     ms_Event->SetHYD2TrkSixthStageTTime(1);
                     ms_Event->SetHYD2TrkSixthStageTDetectorNbr(m_index["Square1"] + N);
                     ms_Event->SetHYD2TrkSixthStageEDetectorNbr(m_index["Square1"] + N);
                  }

                  SixthStageEnergy_itr++;
               }

         DetectorNumber_itr++;
      }

      // clear map for next event
      DetectorNumberHitMap ->clear();
      EnergyHitMap   ->clear()   ;
      TimeHitMap     ->clear()   ;
      XHitMap        ->clear()   ;
      YHitMap        ->clear()   ;
      PosXHitMap     ->clear();
      PosYHitMap     ->clear();
      PosZHitMap     ->clear();
      AngThetaHitMap ->clear();
      AngPhiHitMap   ->clear();
      SecondStageEnergyHitMap ->clear()  ;
      ThirdStageEnergyHitMap ->clear() ;
      FourthStageEnergyHitMap ->clear() ;
      FifthStageEnergyHitMap ->clear() ;
      SixthStageEnergyHitMap ->clear() ;
   }
}



void Hyde2TrackerSquare1::InitializeScorers()
{
  bool already_exist = false; 
  m_FirstStageScorer = NPS::VDetector::CheckScorer("FirstStageScorerHYD2Square1",already_exist);
  m_SecondStageScorer = NPS::VDetector::CheckScorer("SecondStageScorerHYD2Square1",already_exist);
  m_ThirdStageScorer = NPS::VDetector::CheckScorer("ThirdStageScorerHYD2Square1",already_exist);
  m_FourthStageScorer = NPS::VDetector::CheckScorer("FourthStageScorerHYD2Square1",already_exist);
  m_FifthStageScorer = NPS::VDetector::CheckScorer("FifthStageScorerHYD2Square1",already_exist);
  m_SixthStageScorer = NPS::VDetector::CheckScorer("SixthStageScorerHYD2Square1",already_exist);
  if(already_exist) return;

   // First stage Associate Scorer
   G4VPrimitiveScorer* DetNbr                           = new OBSOLETEGENERALSCORERS::PSDetectorNumber("DetectorNumber", "HYD2Square1", 0);
   G4VPrimitiveScorer* TOF                              = new OBSOLETEGENERALSCORERS::PSTOF("StripTime","HYD2Square1", 0);
   G4VPrimitiveScorer* InteractionCoordinatesX          = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesX("InterCoordX","HYD2Square1", 0);
   G4VPrimitiveScorer* InteractionCoordinatesY          = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesY("InterCoordY","HYD2Square1", 0);
   G4VPrimitiveScorer* InteractionCoordinatesZ          = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesZ("InterCoordZ","HYD2Square1", 0);
   G4VPrimitiveScorer* InteractionCoordinatesAngleTheta = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesAngleTheta("InterCoordAngTheta","HYD2Square1", 0);
   G4VPrimitiveScorer* InteractionCoordinatesAnglePhi   = new OBSOLETEGENERALSCORERS::PSInteractionCoordinatesAnglePhi("InterCoordAngPhi","HYD2Square1", 0);
   G4VPrimitiveScorer* Energy                           = new HYD2ScorerFirstStageEnergy("StripEnergy", "HYD2Square1", 0);
   G4VPrimitiveScorer* StripPositionX                   = new HYD2ScorerFirstStageFrontStripSquare1("StripNumberX", 0, NumberOfStrips);
   G4VPrimitiveScorer* StripPositionY                   = new HYD2ScorerFirstStageBackStripSquare1("StripNumberY", 0, NumberOfStrips);

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
   G4VPrimitiveScorer* SecondStageEnergy = new HYD2ScorerSecondStageEnergy("SecondStageEnergy", "HYD2Square1", 0);
   m_SecondStageScorer->RegisterPrimitive(SecondStageEnergy);

   //  Third stage Associate Scorer 
   G4VPrimitiveScorer* ThirdStageEnergy = new HYD2ScorerThirdStageEnergy("ThirdStageEnergy", "HYD2Square1", 0);
   m_ThirdStageScorer->RegisterPrimitive(ThirdStageEnergy);

   //  Fourth stage Associate Scorer 
   G4VPrimitiveScorer* FourthStageEnergy = new HYD2ScorerFourthStageEnergy("FourthStageEnergy", "HYD2Square1", 0);
   m_FourthStageScorer->RegisterPrimitive(FourthStageEnergy);

   //  Fifth stage Associate Scorer 
   G4VPrimitiveScorer* FifthStageEnergy = new HYD2ScorerFifthStageEnergy("FifthStageEnergy", "HYD2Square1", 0);
   m_FifthStageScorer->RegisterPrimitive(FifthStageEnergy);

   //  Sixth stage Associate Scorer 
   G4VPrimitiveScorer* SixthStageEnergy = new HYD2ScorerSixthStageEnergy("SixthStageEnergy", "HYD2Square1", 0);
   m_SixthStageScorer->RegisterPrimitive(SixthStageEnergy);

   //  Add All Scorer to the Global Scorer Manager
   G4SDManager::GetSDMpointer()->AddNewDetector(m_FirstStageScorer);
   G4SDManager::GetSDMpointer()->AddNewDetector(m_SecondStageScorer);
   G4SDManager::GetSDMpointer()->AddNewDetector(m_ThirdStageScorer);
   G4SDManager::GetSDMpointer()->AddNewDetector(m_FourthStageScorer);
   G4SDManager::GetSDMpointer()->AddNewDetector(m_FifthStageScorer);
   G4SDManager::GetSDMpointer()->AddNewDetector(m_SixthStageScorer);
}
