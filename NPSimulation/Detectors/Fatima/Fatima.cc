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
 *  This class describe the Fatima scintillator array                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
using namespace std;

//Geant4
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"

// NPS
#include "Fatima.hh"
using namespace FATIMA;

#include "CalorimeterScorers.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
// NPL
#include "NPOptionManager.h"
#include "RootOutput.h"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Fatima Specific Method
Fatima::Fatima(){
  m_Event = new TFatimaData();

  // Blue
  m_LaBr3VisAtt = new G4VisAttributes(G4Colour(0, 0.5, 1));

  // Dark Grey
  m_PMTVisAtt = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1));

  // Grey wireframe
  m_DetectorCasingVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5,0.2));

  m_LogicalDetector = 0;
  m_LaBr3Scorer = 0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Fatima::~Fatima(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Fatima::AddDetector(G4ThreeVector Pos1, G4ThreeVector Pos2, G4ThreeVector Pos3, G4ThreeVector Pos4){
  G4ThreeVector Pos=(Pos1+Pos2+Pos3+Pos4)/4.;
  G4ThreeVector u = Pos1-Pos2;
  G4ThreeVector v = Pos1-Pos4;
  u = u.unit(); v = v.unit();
  G4ThreeVector w = Pos.unit();
  Pos = Pos + w*Length*0.5;

  m_Pos.push_back(Pos);
  m_Rot.push_back(new G4RotationMatrix(u,v,w));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Fatima::AddDetector(G4ThreeVector Pos, double beta_u, double beta_v, double beta_w){
  double Theta = Pos.theta();
  double Phi = Pos.phi();

  // vector parallel to one axis of silicon plane
  G4double ii = cos(Theta / rad) * cos(Phi / rad);
  G4double jj = cos(Theta / rad) * sin(Phi / rad);
  G4double kk = -sin(Theta / rad);
  G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

  G4ThreeVector w = Pos.unit();
  G4ThreeVector u = w.cross(Y);
  G4ThreeVector v = w.cross(u);
  v = v.unit();
  u = u.unit();

  G4RotationMatrix* r = new G4RotationMatrix(u,v,w);
  r->rotate(beta_u,u);
  r->rotate(beta_v,v);
  r->rotate(beta_w,w);

  m_Pos.push_back(Pos);
  m_Rot.push_back(r);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Fatima::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Fatima");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 
  for(unsigned int i  = 0 ; i < blocks.size() ; i++){
    // Cartesian Case
    vector<string> cart = {"A","B","C","D"};

    // Spherical Case
    vector<string> sphe= {"R","THETA","PHI","BETA"};

    if(blocks[i]->HasTokenList(cart)){
      cout << endl << "////  Fatima " << i+1 << endl;
      G4ThreeVector A = NPS::ConvertVector(blocks[i]->GetTVector3("A","mm"));
      G4ThreeVector B = NPS::ConvertVector(blocks[i]->GetTVector3("B","mm"));
      G4ThreeVector C = NPS::ConvertVector(blocks[i]->GetTVector3("C","mm"));
      G4ThreeVector D = NPS::ConvertVector(blocks[i]->GetTVector3("D","mm"));
      AddDetector(A,B,C,D) ;
    }

    else if(blocks[i]->HasTokenList(sphe)){
      cout << endl << "////  Fatima " << i+1 << endl;
      double Theta = blocks[i]->GetDouble("THETA","deg");
      double Phi= blocks[i]->GetDouble("PHI","deg");
      double R = blocks[i]->GetDouble("R","mm");
      vector<double> beta = blocks[i]->GetVectorDouble("BETA","deg");
      R = R +  0.5*Length;
      G4ThreeVector Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
      AddDetector(Pos,beta[0],beta[1],beta[2]);
    }

    else{
      cout << "ERROR: Missing token for Fatima blocks, check your input file" << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Fatima::ConstructDetector(G4LogicalVolume* world){
	/*
       G4Material* Lead = MaterialManager::getInstance()->GetMaterialFromLibrary("Pb");
       G4Material* Chamber = MaterialManager::getInstance()->GetMaterialFromLibrary("Al" );
      //G4Material* Chamber = MaterialManager::getInstance()->GetMaterialFromLibrary("steel" );
// chamber
    // G4ThreeVector  positionchamber = G4ThreeVector(0, 0, chamber_PosZ);
    // G4Box* solidExtChamber
    // = new G4Box("solidExtChamber", m_ChamberWmax/2, m_ChamberHmax/2, m_ChamberDmax/2 );
   //G4Box* solidIntChamber
    // = new G4Box("solidIntChamber", m_ChamberWmin/2, m_ChamberHmin/2, m_ChamberDmin/2 );


     G4Box* solidExtChamber
     = new G4Box("solidExtChamber", 150*mm/2, 150*mm/2, 550*mm/2 );
   G4Box* solidIntChamber
     = new G4Box("solidIntChamber", (150*mm-1.0*mm)/2, (150*mm-1.0*mm)/2, (550*mm-1.0*mm)/2 );



   G4SubtractionSolid* solidChamber=new G4SubtractionSolid("solidChamber",solidExtChamber, solidIntChamber, 0, G4ThreeVector(0.,0.,-0.5*cm));


   G4LogicalVolume* logicChamber = new G4LogicalVolume(solidChamber, Chamber, "logicChamber" , 0,0,0);

         // rotation of target
         G4RotationMatrix *r = new G4RotationMatrix();
         //rotation->rotateY(m_ChamberAngle);

         
       new G4PVPlacement(r, G4ThreeVector(0., -1.92325*cm/5.5, 0.), logicChamber, "solidChamber",world , false, 0); 

        G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(.66, .66, .66));
         logicChamber->SetVisAttributes(ChamberVisAtt); 


//holding support left  (ring shaped)


             G4ThreeVector  positionHolder2 = G4ThreeVector(0, 0, 11*cm);

    G4Tubs* solidHolder2 = new G4Tubs("solidHolder2", 29*cm, 31.5*cm, 46*mm, 0.*deg, 360.*deg); 
    G4LogicalVolume* logicHolder2 = new G4LogicalVolume(solidHolder2, Lead, "logicHolder2", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder2, 
        logicHolder2, 
        "Fatima_Holder2", 
        world, 
        false, 
        0);

     G4VisAttributes* Holder2VisAtt = new G4VisAttributes(G4Colour(.65, .65, .65));
    logicHolder2->SetVisAttributes(Holder2VisAtt);
  

  // holding support right (ring shaped)

        G4ThreeVector  positionHolder3 = G4ThreeVector(0, 0, -14*cm);

    G4Tubs* solidHolder3 = new G4Tubs("solidHolder3", 29*cm, 31.5*cm, 42*mm, 0.*deg, 360.*deg); 
    G4LogicalVolume* logicHolder3 = new G4LogicalVolume(solidHolder3, Lead, "logicHolder3", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder3, 
        logicHolder3, 
        "Fatima_Holder3", 
        world, 
        false, 
        0);

     G4VisAttributes* Holder3VisAtt = new G4VisAttributes(G4Colour(.65, .65, .65));
    logicHolder3->SetVisAttributes(Holder3VisAtt);

    //support from the bottom rod front

              G4ThreeVector  positionHolder4 = G4ThreeVector(25*cm, -30*cm,12*cm);

    G4Box* solidHolder4 = new G4Box("solidHolder4", 40*mm, 140*mm, 20*mm); 
    G4LogicalVolume* logicHolder4 = new G4LogicalVolume(solidHolder4, Lead, "logicHolder4", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder4, 
        logicHolder4, 
        "Fatima_Holder4", 
        world, 
        false, 
        0);


     G4VisAttributes* Holder4VisAtt = new G4VisAttributes(G4Colour(.65, .65, .65));
    logicHolder4->SetVisAttributes(Holder4VisAtt);
      
        //support from the bottom rod front

             G4ThreeVector  positionHolder5 = G4ThreeVector(25*cm,-30*cm,-15*cm);

    G4Box* solidHolder5 = new G4Box("solidHolder5", 40*mm, 140*mm, 20*mm); 
    G4LogicalVolume* logicHolder5 = new G4LogicalVolume(solidHolder5, Lead, "logicHolder5", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder5, 
        logicHolder5, 
        "Fatima_Holder5", 
        world, 
        false, 
        0);


     G4VisAttributes* Holder5VisAtt = new G4VisAttributes(G4Colour(.65, .65, .65));
    logicHolder5->SetVisAttributes(Holder5VisAtt);

     // support from the bottom back1

    G4ThreeVector  positionHolder6 = G4ThreeVector(-25*cm,-30*cm,12*cm);

    G4Box* solidHolder6 = new G4Box("solidHolder6", 40*mm, 140*mm, 20*mm); 
    G4LogicalVolume* logicHolder6 = new G4LogicalVolume(solidHolder6, Lead, "logicHolder6", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder6, 
        logicHolder6, 
        "Fatima_Holder6", 
        world, 
        false, 
        0);


     G4VisAttributes* Holder6VisAtt = new G4VisAttributes(G4Colour(.65, .65, .65));
    logicHolder6->SetVisAttributes(Holder6VisAtt);

    // support from the bottom back2

    G4ThreeVector  positionHolder7 = G4ThreeVector(-25*cm,-30*cm,-12*cm);

    G4Box* solidHolder7 = new G4Box("solidHolder7", 40*mm, 140*mm, 20*mm); 
    G4LogicalVolume* logicHolder7 = new G4LogicalVolume(solidHolder7, Lead, "logicHolder7", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder7, 
        logicHolder7, 
        "Fatima_Holder7", 
        world, 
        false, 
        0);


     G4VisAttributes* Holder7VisAtt = new G4VisAttributes(G4Colour(.65, .65, .65));
    logicHolder7->SetVisAttributes(Holder7VisAtt);
*/






  unsigned int mysize = m_Pos.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    new G4PVPlacement(G4Transform3D(*m_Rot[i], m_Pos[i]), ConstructDetector(),  "FatimaDetector", world, false, i); 
    //new G4PVPPlacement(G4Transform3D(*m_Rot[i], m_Pos[i], constructdetector(), "chamber", world, false, i+1);

//new G4PVPlacement(0, G4ThreeVector(0., -1.92325*cm/5.5, 0.), logicChamber, "solidChamber",world , false, i+1);





  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Fatima::ConstructDetector(){
  if(!m_LogicalDetector){

    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    G4Material* Alu = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* Lead = MaterialManager::getInstance()->GetMaterialFromLibrary("Pb");
    G4Material* LaBr3 = MaterialManager::getInstance()->GetMaterialFromLibrary("LaBr3_Ce");
  // G4Material* chamber = MaterialManager::getInstance()->GetMaterialFromLibrary("Steel");

    // Mother Volume
    G4Tubs* solidFatimaDetector = 
      new G4Tubs("Fatima",0, 0.5*FaceFront, 0.5*Length, 0.*deg, 360.*deg);
    m_LogicalDetector = 
      new G4LogicalVolume(solidFatimaDetector, Vacuum, "Fatima", 0, 0, 0);

    m_LogicalDetector->SetVisAttributes(G4VisAttributes::Invisible);

    // Detector construction
    // LaBr3
    G4ThreeVector  positionLaBr3 = G4ThreeVector(0, 0, LaBr3_PosZ);

    G4Tubs* solidLaBr3 = new G4Tubs("solidLaBr3", 0., 0.5*LaBr3Face, 0.5*LaBr3Thickness, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLaBr3 = new G4LogicalVolume(solidLaBr3, LaBr3, "logicLaBr3", 0, 0, 0);

    new G4PVPlacement(0, 
        positionLaBr3, 
        logicLaBr3, 
        "Fatima_LaBr3", 
        m_LogicalDetector, 
        false, 
        0);



  /*  // chamber construction
    // chamber
     G4ThreeVector  positionchamber = G4ThreeVector(0, 0, chamber_PosZ);
    // G4Box* solidExtChamber
    // = new G4Box("solidExtChamber", m_ChamberWmax/2, m_ChamberHmax/2, m_ChamberDmax/2 );
   //G4Box* solidIntChamber
    // = new G4Box("solidIntChamber", m_ChamberWmin/2, m_ChamberHmin/2, m_ChamberDmin/2 );


     G4Box* solidExtChamber
     = new G4Box("solidExtChamber", 150*mm, 150*mm, 150*mm );
   G4Box* solidIntChamber
     = new G4Box("solidIntChamber", 152*mm, 152*mm, 152*mm );



   G4SubtractionSolid* solidChamber=new G4SubtractionSolid("SolidChamber",solidExtChamber, solidIntChamber, 0, G4ThreeVector(0.,0.,-0.5*cm));


   G4LogicalVolume* logicChamber = new G4LogicalVolume(solidChamber, chamber, "logicChamber");

         // rotation of target
         //G4RotationMatrix *rotation = new G4RotationMatrix();
         //rotation->rotateY(m_ChamberAngle);

         
            new G4PVPlacement(0, G4ThreeVector(0., -1.92325*cm/5.5, 0.), logicChamber, "chamber",m_LogicalDetector , false, 0);

        G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(.66, .66, .66));
         logicChamber->SetVisAttributes(ChamberVisAtt);  */





    // Set LaBr3 sensible
    logicLaBr3->SetSensitiveDetector(m_LaBr3Scorer);

    // Visualisation of LaBr3 Strip
    logicLaBr3->SetVisAttributes(m_LaBr3VisAtt);

    // Aluminium can around LaBr3
    // LaBr3 Can
    G4ThreeVector  positionLaBr3Can = G4ThreeVector(0, 0, LaBr3Can_PosZ);

    G4Tubs* solidLaBr3Can = new G4Tubs("solidLaBr3Can", 0.5*CanInnerDiameter, 0.5*CanOuterDiameter, 0.5*CanLength, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLaBr3Can = new G4LogicalVolume(solidLaBr3Can, Alu, "logicLaBr3Can", 0, 0, 0);

    new G4PVPlacement(0, 
        positionLaBr3Can, 
        logicLaBr3Can, 
        "Fatima_LaBr3Can", 
        m_LogicalDetector, 
        false, 
        0);

    // Visualisation of LaBr3Can
    logicLaBr3Can->SetVisAttributes(m_DetectorCasingVisAtt);

    // Aluminium window in front of LaBr3
    // LaBr3 Window
    G4ThreeVector  positionLaBr3Win = G4ThreeVector(0, 0, LaBr3Win_PosZ);

    G4Tubs* solidLaBr3Win = new G4Tubs("solidLaBr3Win", 0.5*WinInnerDiameter, 0.5*WinOuterDiameter, 0.5*WinLength, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLaBr3Win = new G4LogicalVolume(solidLaBr3Win, Alu, "logicLaBr3Win", 0, 0, 0);

    new G4PVPlacement(0, 
        positionLaBr3Win, 
        logicLaBr3Win, 
        "Fatima_LaBr3Win", 
        m_LogicalDetector, 
        false, 
        0);

    // Visualisation of LaBr3Win
    logicLaBr3Win->SetVisAttributes(m_DetectorCasingVisAtt);

    // PMT
    G4ThreeVector  positionPMT = G4ThreeVector(0, 0, PMT_PosZ);

   G4Tubs* solidPMout = new G4Tubs("solidPMOut", 0.5*LaBr3Face, 0.5*PMTFace, 0.5*PMTThickness, 0.*deg, 360.*deg);
    G4Tubs* solidPMin = new G4Tubs("solidPMIn", 0.5*LaBr3Face-0.1*cm, 0.5*PMTFace-0.5*cm, 0.5*(PMTThickness-2.*cm)-0.1*cm, 0.*deg, 360.*deg);
    G4RotationMatrix* RotMat=NULL;
    const G4ThreeVector &Trans= G4ThreeVector(0.,0.,1.*cm); 
    G4SubtractionSolid*           solidPMT = new G4SubtractionSolid("solidPMT", solidPMout,solidPMin, RotMat, Trans);

    G4LogicalVolume* logicPMT = new G4LogicalVolume(solidPMT, Alu, "logicPMT", 0, 0, 0); 

 


    new G4PVPlacement(0, 
        positionPMT, 
        logicPMT, 
        "Fatima_PMT", 
        m_LogicalDetector, 
        false, 
        0);

     

    // Visualisation of PMT Strip
    logicPMT->SetVisAttributes(m_PMTVisAtt);

    // Lead shielding
    // A
    G4ThreeVector  positionLeadAShield = G4ThreeVector(0, 0, LeadAShield_PosZ);
    G4Tubs* solidLeadA = new G4Tubs("solidLead", 0.5*LeadAMinR, 0.5*LeadAMaxR, 0.5*LeadALength, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLeadAShield = new G4LogicalVolume(solidLeadA, Lead, "logicLeadAShield", 0, 0, 0);

    new G4PVPlacement(0, 
        positionLeadAShield, 
        logicLeadAShield, 
        "Fatima_LeadAShield", 
        m_LogicalDetector, 
        false, 
        0);
    // B
    G4ThreeVector  positionLeadBShield = G4ThreeVector(0, 0, LeadBShield_PosZ);
    G4Tubs*           solidLeadB = new G4Tubs("solidLead", 0.5*LeadBMinR, 0.5*LeadBMaxR, 0.5*LeadBLength, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLeadBShield = new G4LogicalVolume(solidLeadB, Lead, "logicLeadBShield", 0, 0, 0);

    new G4PVPlacement(0, 
        positionLeadBShield, 
        logicLeadBShield, 
        "Fatima_LeadBShield", 
        m_LogicalDetector, 
        false, 
        0);

    // Visualisation of PMT Strip
    G4VisAttributes* LeadVisAtt = new G4VisAttributes(G4Colour(.66, .66, .66));                                  //previous in Fatima: (G4Colour(1., 1., 0.));
    logicLeadAShield->SetVisAttributes(LeadVisAtt);
    logicLeadBShield->SetVisAttributes(LeadVisAtt);
    
    
    
    
    //LeadShield on PMT


G4ThreeVector  positionPMTShield = G4ThreeVector(0, 0, PMTShield_PosZ);
    G4Tubs* solidShield = new G4Tubs("solidShield", PMTIn*0.5, PMTOut*0.5,PMTThickness *0.5, 0.*deg, 360.*deg);
    G4LogicalVolume* logicPMTShield = new G4LogicalVolume(solidShield, Alu, "logicPMTShield", 0, 0, 0);

 /*   new G4PVPlacement(0, 
        positionPMTShield, 
        logicPMTShield, 
        "Fatima_PMTShield", 
        m_LogicalDetector, 
        false, 
        0); */

    G4VisAttributes* PMTShieldVisAtt = new G4VisAttributes(G4Colour(.66, .66, .66));
    logicPMTShield->SetVisAttributes(PMTShieldVisAtt);




    
    
    
    
    
    
    
    
    
    
    
    
    
    //   Backsideof the shield
     G4ThreeVector  positionShieldCan = G4ThreeVector(0, 0, ShieldCan_PosZ);

    G4Tubs* solidShieldCan = new G4Tubs("solidShieldCan",  0.5*WinInnerDiameter, 1*WinOuterDiameter, 0.5*WinLength, 0.*deg, 360.*deg);
    G4LogicalVolume* logicShieldCan = new G4LogicalVolume(solidShieldCan, Alu, "logicShieldCan", 0, 0, 0);

    new G4PVPlacement(0, 
        positionShieldCan, 
        logicShieldCan, 
        "Fatima_ShieldCan", 
        m_LogicalDetector, 
        false, 
        0);

    // Visualisation of LaBr3Can
     G4VisAttributes* ShieldCanVisAtt = new G4VisAttributes(G4Colour(.55, .6, .6));
    logicShieldCan->SetVisAttributes(ShieldCanVisAtt);
  
/*
    // HolderCan
  G4ThreeVector  positionHolderCan = G4ThreeVector(0, 0, HolderCan_PosZ);

    G4Tubs* solidHolderCan = new G4Tubs("solidHolderCan", 0.8*LeadBMinR, 0.8*LeadBMaxR, 0.8*LeadBLength, 0.*deg, 360.*deg); 
    G4LogicalVolume* logicHolderCan = new G4LogicalVolume(solidHolderCan, Alu, "logicHolderCan", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolderCan, 
        logicHolderCan, 
        "Fatima_HolderCan", 
        m_LogicalDetector, 
        false, 
        0);

    // Visualisation of LaBr3Can
     G4VisAttributes* HolderCanVisAtt = new G4VisAttributes(G4Colour(.4, .4, .4));
    logicHolderCan->SetVisAttributes(HolderCanVisAtt);
  
    
    
    //holding1
    
   G4ThreeVector  positionHolder = G4ThreeVector(0, -3*cm, 0);

    G4Box* solidHolder = new G4Box("solidHolder", 20*mm, 50*mm, 10*mm); 
    G4LogicalVolume* logicHolder = new G4LogicalVolume(solidHolder, Lead, "logicHolder", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder, 
        logicHolder, 
        "Fatima_Holder", 
        m_LogicalDetector, 
        false, 
        0);


     G4VisAttributes* HolderVisAtt = new G4VisAttributes(G4Colour(.4, .4, .4));
    logicHolder->SetVisAttributes(HolderVisAtt);
  
   //holding2
                       G4ThreeVector  positionHolder1 = G4ThreeVector(2*cm, 5*cm, 2*cm);

    G4Box* solidHolder1 = new G4Box("solidHolder1", 16*mm, 22*mm, 8*mm); 
    G4LogicalVolume* logicHolder1 = new G4LogicalVolume(solidHolder1, Lead, "logicHolder", 0, 0, 0);

    new G4PVPlacement(0, 
        positionHolder1, 
        logicHolder1, 
        "Fatima_Holder1", 
        m_LogicalDetector, 
        false, 
        0);


     G4VisAttributes* Holder1VisAtt = new G4VisAttributes(G4Colour(.65, .65, .65));
    logicHolder1->SetVisAttributes(HolderVisAtt);
*/
         
    
  }

  return m_LogicalDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Fatima::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Fatima")){
    pTree->Branch("Fatima", "TFatimaData", &m_Event) ;
  } 
  pTree->SetBranchAddress("Fatima", &m_Event) ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Fatima::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // LaBr3
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_LaBr3Scorer->GetPrimitive(0);
  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 

    
 
     double E = Scorer->GetEnergy(i); 
      double sigma = (0.9207*pow(E*1000,0.4055));
      // double sigma = (2.154*pow(E*1000,0.6244))/2.35;
      double Energy = RandGauss::shoot(E*1000,sigma)/1000;

                               //introduced this line from nana
   // double Energy = RandGauss::shoot(E,(E*0.0325637)/(2.35*pow(E-0.00975335,0.475759)));        // [previous notation after shoot](Scorer->GetEnergy(i), EnergyResolution);

    if(Energy>EnergyThreshold){
      double Time = Scorer->GetTime(i);
      int DetectorNbr =Scorer-> GetLevel(i)[0];

      m_Event->SetFatimaLaBr3E(DetectorNbr,Energy);
      m_Event->SetFatimaLaBr3T(DetectorNbr,Time);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Fatima::InitializeScorers(){
  vector<G4int> NestingLevel;
  NestingLevel.push_back(1);

  //   LaBr3 Associate Scorer
  bool already_exist = false;
  m_LaBr3Scorer = CheckScorer("Fatima_LaBr3Scorer",already_exist);

  // if the scorer were created previously nothing else need to be made
  if(already_exist) return;

  G4VPrimitiveScorer* LaBr3Scorer =
    new  CalorimeterScorers::PS_Calorimeter("FatimaLaBr3",NestingLevel);
  //and register it to the multifunctionnal detector
  m_LaBr3Scorer->RegisterPrimitive(LaBr3Scorer);

  //   Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_LaBr3Scorer) ;
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Fatima::Construct(){
  return  (NPS::VDetector*) new Fatima();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_fatima{
    public:
      proxy_nps_fatima(){
        NPS::DetectorFactory::getInstance()->AddToken("Fatima","Fatima");
        NPS::DetectorFactory::getInstance()->AddDetector("Fatima",Fatima::Construct);
      }
  };

  proxy_nps_fatima p_nps_fatima;
}
