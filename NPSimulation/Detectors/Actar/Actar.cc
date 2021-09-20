/*****************************************************************************
 * Copyright (C) 2009-2017   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: morfouac@nscl.msu.edu *
 *                                                                           *
 * Creation Date  : September 2017                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Actar simulation                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
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

// G4 Field
#include "G4FieldManager.hh"
#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4MaterialPropertiesTable.hh"

// NPTool header
#include "Actar.hh"
#include "DSSDScorers.hh"
#include "TPCScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "FastDriftElectron.hh"
#include "NPSHitsMap.hh"
#include "CalorimeterScorers.hh"


// CLHEP header
#include "CLHEP/Random/RandGauss.h"

// ROOT
#include "TH1D.h"
#include "TF1.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Actar_NS{
    // Energy and time Resolution
    const double ChargeThreshold = 0;
    //const double ResoTime = 0.1*ns ;
    const double ResoCharge = 5./100 ;
    const double ChamberThickness = 376*mm ;
    const double ChamberWidth = 376*mm ;
    const double ChamberHeight = 40*cm ;
    //const int NumberOfPads = 16384;
    //const int PadX = 128;
    //const int PadZ = 128;

    const double Nose_Rmin = 2.5*cm;
    const double Nose_Rmax = 3.5*cm;
    const double Nose_Length = 12*cm;

    const double Mylar_Rmax = 3.5*cm;
    const double Mylar_Thickness = 7*micrometer;

    const double XGazVolume = 256.*mm;
    const double YGazVolume = 256.*mm;
    const double ZGazVolume = 256.*mm;

    const double SiliconHeight = 53.*mm;
    const double SiliconWidth = 53.*mm;
    const double SiliconThickness = 0.7*mm;
    const double DistInterSi = 1.*mm;
    const double Si_PosZ=175.*mm;
    const double ResoSilicon = 0.60/2.35;
    const double EnergyThreshold = 0.1;

    const double VamosSiliconHeight = 70*mm;
    const double VamosSiliconWidth = 50*mm;
    const double VamosSiliconThickness = 0.5*mm;
    const double VamosSiliconDistanInterSi = 1*mm;
    const double VamosSilicon_PosZ = -160*mm;

    const double CsIThickness = 1.*cm;
    const double CsIHeight = 2.5*cm;
    const double CsIWidth = 2.5*cm;
    const double DistInterCsI = 1.*mm;
    const double CsI_PosZ = 16.*cm;
    //const double ResoCsI = 0.200/2.35;

    const double BeamDumpRadius = 30*mm;
    const double BeamDumpThickness = 5*mm;
    const double BeamDump_PosZ = 160*mm;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Actar Specific Method
Actar::Actar(){
    m_Event = new TActarData() ;
    m_EventReduced = new MEventReduced();
    m_ActarScorer = 0;
    m_SquareDetector = 0;

    // RGB Color + Transparency
    m_VisChamber        = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.3));
    m_VisWindows        = new G4VisAttributes(G4Colour(1, 0, 0, 0.25));
    m_VisGas            = new G4VisAttributes(G4Colour(0, 0.5, 0.5, 0.3));
    m_VisPads           = new G4VisAttributes(G4Colour(255, 223, 50, 0.8));
    m_VisMicromegas     = new G4VisAttributes(G4Colour(100, 100, 100, 0.4));
    m_SiliconVisAtt     = new G4VisAttributes(G4Colour(0.529412, 0.807843, 0.980392, 0.95)) ;
    m_CsIVisAtt         = new G4VisAttributes(G4Colour(0.429412, 0.607843, 0.780392, 0.95));
    m_BeamDumpVisAtt    = new G4VisAttributes(G4Colour(0.9, 0.5, 0.5));
    m_VisPads->SetForceWireframe(true);

    m_build_BeamDump= 0;
    m_build_Silicon=1;
    m_build_Vamos_Silicon=0;
    m_build_CsI=1;
    m_ReactionRegion=NULL;
    // Lookup table //
   // bool ReadingLookupTable = false;
    // Opening the LookUp Table LT.dat
    G4String GlobalPath = getenv("NPTOOL");
    G4String LT_FileName = GlobalPath + "/NPSimulation/Detectors/Actar/LT.dat";
    //string LT_FileName = "./Detectors/Actar/LT.dat";
    //string LT_FileName = "./configs/LT.dat";
    ifstream LTConfigFile;
    LTConfigFile.open(LT_FileName.c_str());
    if(!LTConfigFile.is_open()){
        cout << "No Lookup Table in " << LT_FileName << " found!" << endl;
        return;
    }
    else{
        cout << "/// Using LookupTable from: " << LT_FileName << " ///" << endl;
        int co, as, ag, ch;
        int pX, pY;
        for(int i=0;i<NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel;i++){
            LTConfigFile >> co >> as >> ag >> ch >> pX >> pY;
            if(pX!=-1 && pY !=-1){
                m_PadToCobo[pX][pY] = co;
                m_PadToAsad[pX][pY] = as;
                m_PadToAGET[pX][pY] = ag;
                m_PadToChannel[pX][pY] = ch;
            }
        }
        //ReadingLookupTable = true;
    }
    LTConfigFile.close();

}

Actar::~Actar(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Actar::AddDetector(G4ThreeVector POS, string  Shape){
    // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4
    m_R.push_back(POS.mag());
    m_Theta.push_back(POS.theta());
    m_Phi.push_back(POS.phi());
    m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Actar::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
    m_R.push_back(R);
    m_Theta.push_back(Theta);
    m_Phi.push_back(Phi);
    m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Actar::BuildDetector(){

    G4Material* Cu= MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
    G4Material* Si= MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    G4Material* Al= MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* Mylar= MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");
    G4Material* MaterialCsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");

    if(!m_SquareDetector){
        // Main volume
        G4Box* sChamber = new G4Box("Actar_Box",Actar_NS::ChamberWidth*0.5,
                                    Actar_NS::ChamberHeight*0.5,Actar_NS::ChamberThickness*0.5);

        //Nose volume
        G4Tubs* sNose = new G4Tubs("Actar_Nose",Actar_NS::Nose_Rmin, Actar_NS::Nose_Rmax,Actar_NS::Nose_Length*0.5,
                                   0*deg, 360*deg);

        //Mylar volume
        G4Tubs* sWindows = new G4Tubs("Actar_Windows",0, Actar_NS::Mylar_Rmax,Actar_NS::Mylar_Thickness*0.5,
                                      0*deg, 360*deg);

        // Cage volume
        G4Box* sCage = new G4Box("Actar_Gas",Actar_NS::XGazVolume*0.5,
                                 Actar_NS::YGazVolume*0.5,Actar_NS::ZGazVolume*0.5);

        // Pad
        G4Box* sPad = new G4Box("Actar_Pad",256*mm*0.5,
                                2*mm*0.5,256*mm*0.5);

        // Micromegas
        G4Box* sMicromegas = new G4Box("Actar_Micromegas",256*mm*0.5,
                                       220*micrometer*0.5,256*mm*0.5);

        // Cathode
        G4Box* sCathode = new G4Box("Actar_Cathode",26.5*cm*0.5,
                                    1*um*0.5,25.6*cm*0.5);



        unsigned const int NumberOfGasMix = m_GasMaterial.size();

        double density=0;
        double density_sum=0;
        vector<G4Material*> GasComponent;
        vector<double> FractionMass;

        for(unsigned int i=0; i<NumberOfGasMix; i++){
            GasComponent.push_back(MaterialManager::getInstance()->GetGasFromLibrary(m_GasMaterial[i],m_Pressure,m_Temperature) );
        }
        for(unsigned int i=0; i<NumberOfGasMix; i++){
            density     += ((double)m_GasFraction[i]/100)*GasComponent[i]->GetDensity();
            density_sum += GasComponent[i]->GetDensity();
        }
        //cout << "density = " << density*cm3/g << endl;

        for(unsigned int i=0; i<NumberOfGasMix; i++){
            FractionMass.push_back(GasComponent[i]->GetDensity()/density_sum);
        }

        G4Material* GasMaterial = new G4Material("GasMix", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);
        G4Material* DriftGasMaterial = new G4Material("DriftGasMix", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);

        for(unsigned int i=0; i<NumberOfGasMix; i++){
            GasMaterial->AddMaterial(GasComponent[i], FractionMass[i]);
            DriftGasMaterial->AddMaterial(GasComponent[i], FractionMass[i]);
            cout << GasComponent[i] << endl;
        }

        G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
        MPT->AddConstProperty("DE_PAIRENERGY",20*eV);
        MPT->AddConstProperty("DE_YIELD",3e-1);
        //MPT->AddConstProperty("DE_AMPLIFICATION",2);
        MPT->AddConstProperty("DE_ABSLENGTH",1*pc);
        MPT->AddConstProperty("DE_DRIFTSPEED",0.8*cm/microsecond);
        MPT->AddConstProperty("DE_TRANSVERSALSPREAD",2e-5*mm2/ns);
        MPT->AddConstProperty("DE_LONGITUDINALSPREAD",4e-6*mm2/ns);

        DriftGasMaterial->SetMaterialPropertiesTable(MPT);

        G4MaterialPropertiesTable* MPT2 = new G4MaterialPropertiesTable();
        MPT2->AddConstProperty("DE_AMPLIFICATION",1000);
        MPT2->AddConstProperty("DE_ABSLENGTH",1*pc);

        Al->SetMaterialPropertiesTable(MPT2);

        m_SquareDetector    = new G4LogicalVolume(sChamber,GasMaterial,"logic_Actar_Box",0,0,0);
        m_logicGas          = new G4LogicalVolume(sCage,DriftGasMaterial,"logic_Gas",0,0,0);
        G4LogicalVolume* logicPad = new G4LogicalVolume(sPad,Cu,"logic_Pad",0,0,0);
        G4LogicalVolume* logicMicromegas = new G4LogicalVolume(sMicromegas,Al,"logic_Micromegas",0,0,0);

        G4LogicalVolume* logicNose = new G4LogicalVolume(sNose,Cu,"logic_Nose",0,0,0);
        /*G4LogicalVolume* logicCathode = */ new G4LogicalVolume(sCathode,Cu,"logic_Cathode",0,0,0);
        G4LogicalVolume* logicWindows = new G4LogicalVolume(sWindows,Mylar,"logic_Windows",0,0,0);

        G4RotationMatrix* Rot = new G4RotationMatrix();
        //new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,-Actar_NS::ChamberThickness*0.5+Actar_NS::Nose_Length*0.5)),
        //                logicNose,
        //              "ActarNose",m_SquareDetector,false,0);

        //new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,-Actar_NS::ChamberThickness*0.5+Actar_NS::Mylar_Thickness*0.5+Actar_NS::Nose_Length)),
        //                logicWindows,
        //              "ActarEntranceWindows",m_SquareDetector,false,0);

        new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,0)),
                          m_logicGas,
                          "ActarGas",m_SquareDetector,false,0);

        //new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,Actar_NS::YGazVolume*0.5-0.3*cm,0)),
          //                logicMicromegas,
            //              "ActarMicromegas",m_logicGas,false,0);

        new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,Actar_NS::YGazVolume*0.5,0)),
                          logicPad,
                          "ActarPad",m_logicGas,false,0);

        /*int pad=0;
         m_PadToXRow.clear();
         m_PadToZColumn.clear();
         for(int i=0; i<Actar_NS::PadX; i++){
         for(int j=0; j<Actar_NS::PadZ; j++){
         m_PadToXRow[pad] = i;
         m_PadToZColumn[pad] = j;
         double X=(i-64)*2*mm;
         double Z=(j-64)*2*mm;
         new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(X,Actar_NS::YGazVolume*0.5,Z)),
         logicPad,
         "ActarPad",m_logicGas,false,pad+1);

         pad++;
         }
         }*/

        /*new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(3*cm-0.5*1*um,0,0)),
         logicCathode,
         "ActarCathode",m_logicGas,false,0);



         new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0,0,-6*cm+6*micrometer)),
         logicWindows,
         "ActarEntranceWindows",m_SquareDetector,false,0);*/

        G4ElectricField* field = new G4UniformElectricField(G4ThreeVector(0.0,-400*volt/cm,0.0));
        // Create an equation of motion for this field
        G4EqMagElectricField*  Equation = new G4EqMagElectricField(field);
        G4MagIntegratorStepper* Stepper = new G4ClassicalRK4( Equation, 8 );

        // Get the global field manager
        G4FieldManager* FieldManager= new G4FieldManager();
        // Set this field to the global field manager
        FieldManager->SetDetectorField(field);
        m_logicGas->SetFieldManager(FieldManager,true);

        G4MagInt_Driver* IntgrDriver = new G4MagInt_Driver(0.1*mm,
                                                           Stepper,
                                                           Stepper->GetNumberOfVariables() );

        G4ChordFinder* ChordFinder = new G4ChordFinder(IntgrDriver);
        FieldManager->SetChordFinder( ChordFinder );


        logicPad->SetSensitiveDetector(m_ActarScorer);
        logicNose->SetVisAttributes(m_VisChamber);
        m_SquareDetector->SetVisAttributes(m_VisChamber);
        m_logicGas->SetVisAttributes(m_VisGas);
        logicWindows->SetVisAttributes(m_VisWindows);
        logicPad->SetVisAttributes(m_VisPads);
        logicMicromegas->SetVisAttributes(m_VisMicromegas);
        //m_SquareDetector->SetSensitiveDetector(m_ActarScorer);
    }


    ////////////////////////////////////////////////////
    ///////////////////// Beam Dump ////////////////////
    ////////////////////////////////////////////////////
    if(m_build_BeamDump){
        G4Tubs* sBeamDump = new G4Tubs("Actar_BeamDump",0*mm, Actar_NS::BeamDumpRadius,Actar_NS::BeamDumpThickness*0.5,
                                   0*deg, 360*deg);
        m_LogicBeamDump = new G4LogicalVolume(sBeamDump, Cu, "logicBeamDump",0,0,0);

        G4ThreeVector positionBeamDump = G4ThreeVector(0, 0, Actar_NS::BeamDump_PosZ);
        new G4PVPlacement(new G4RotationMatrix(0,0,0),
                          positionBeamDump,
                          m_LogicBeamDump,"BeamDump",
                          m_SquareDetector,false,0);

        m_LogicBeamDump->SetVisAttributes(m_BeamDumpVisAtt) ;

    }
    ///////////////////////////////////////////////////
    ///////////////////// Thin Si /////////////////////
    ///////////////////////////////////////////////////
    int SiliconNumber=0;
    if(m_build_Silicon){
        G4Box* solidSi = new G4Box("Si", 0.5*Actar_NS::SiliconWidth, 0.5*Actar_NS::SiliconHeight, 0.5*Actar_NS::SiliconThickness);	;
        m_LogicSilicon = new G4LogicalVolume(solidSi, Si, "logicSi", 0, 0, 0);

        for(int k=0;k<4; k++){
            for(int p=0; p<5; p++){
                double PosX;
                double PosY;
                if(k==0) PosY= -1.5*Actar_NS::SiliconHeight-2*Actar_NS::DistInterSi;
                if(k==1) PosY= -0.5*Actar_NS::SiliconHeight-1*Actar_NS::DistInterSi;
                if(k==2) PosY= 0.5*Actar_NS::SiliconHeight+1*Actar_NS::DistInterSi;
                if(k==3) PosY= 1.5*Actar_NS::SiliconHeight+2*Actar_NS::DistInterSi;
                if(p==0) PosX= -2*Actar_NS::SiliconWidth-2*Actar_NS::DistInterSi;
                if(p==1) PosX= -1*Actar_NS::SiliconWidth-1*Actar_NS::DistInterSi;
                if(p==2) PosX= 0;
                if(p==3) PosX= 1*Actar_NS::SiliconWidth+2*Actar_NS::DistInterSi;
                if(p==4) PosX= 2*Actar_NS::SiliconWidth+1*Actar_NS::DistInterSi;

                G4ThreeVector positionSi = G4ThreeVector(PosX, PosY, Actar_NS::Si_PosZ);
                new G4PVPlacement(new G4RotationMatrix(0,0,0),
                                  positionSi,
                                  m_LogicSilicon,"Si",
                                  m_SquareDetector,false,SiliconNumber);
                SiliconNumber++;
            }
        }

        // Set Si sensible
        m_LogicSilicon->SetSensitiveDetector(m_SiliconScorer);

        // Visualisation of ThinSi
        m_LogicSilicon->SetVisAttributes(m_SiliconVisAtt);
    }

    ////////////////////////////////////////////////////
    ///////////////////// Vamos Si /////////////////////
    ////////////////////////////////////////////////////
    if(m_build_Vamos_Silicon){
        G4Box* solidSi = new G4Box("Si", 0.5*Actar_NS::VamosSiliconWidth, 0.5*Actar_NS::VamosSiliconHeight, 0.5*Actar_NS::VamosSiliconThickness);	;
        m_LogicVamosSilicon = new G4LogicalVolume(solidSi, Si, "logicVamosSi", 0, 0, 0);


        int VamosSiliconNumber=0;
        for(int k=0;k<4; k++){
            for(int p=0; p<3; p++){
                double PosX;
                double PosY;
                if(k==0) PosX= -35*mm -0.5*Actar_NS::VamosSiliconWidth - Actar_NS::VamosSiliconDistanInterSi;
                if(k==1) PosX= -35*mm -1.5*Actar_NS::VamosSiliconWidth - Actar_NS::VamosSiliconDistanInterSi;
                if(k==2) PosX= 35*mm + 0.5*Actar_NS::VamosSiliconWidth + Actar_NS::VamosSiliconDistanInterSi;
                if(k==3) PosX= 35*mm + 1.5*Actar_NS::VamosSiliconWidth + Actar_NS::VamosSiliconDistanInterSi;

                if(p==0) PosY= -Actar_NS::VamosSiliconHeight-Actar_NS::VamosSiliconDistanInterSi;
                if(p==1) PosY= 0;
                if(p==2) PosY= Actar_NS::VamosSiliconHeight+Actar_NS::VamosSiliconDistanInterSi;

                G4ThreeVector positionSi = G4ThreeVector(PosX, PosY, Actar_NS::VamosSilicon_PosZ);
                new G4PVPlacement(new G4RotationMatrix(0,0,0),
                                  positionSi,
                                  m_LogicVamosSilicon,"Si",
                                  m_SquareDetector,false,SiliconNumber);
                VamosSiliconNumber++;
                SiliconNumber++;
            }
        }

        // Set Si sensible
        m_LogicVamosSilicon->SetSensitiveDetector(m_SiliconScorer);

        // Visualisation of ThinSi
        m_LogicVamosSilicon->SetVisAttributes(m_SiliconVisAtt);
    }

    ///////////////////////////////////////////////
    ///////////////////// CsI /////////////////////
    ///////////////////////////////////////////////
    if(m_build_CsI){
        G4Box* solidCsI = new G4Box("Si", 0.5*Actar_NS::CsIWidth, 0.5*Actar_NS::CsIHeight, 0.5*Actar_NS::CsIThickness);	;
        m_LogicCsICrystal = new G4LogicalVolume(solidCsI, MaterialCsI, "logicCsI", 0, 0, 0);

        int CsINumber=0;
        for(int k=0;k<8; k++){
            for(int p=0; p<10; p++){
                double PosX;
                double PosY;
                if(k<4) PosY= -0.5*Actar_NS::CsIHeight-0.5*Actar_NS::DistInterCsI+(k-3)*(Actar_NS::CsIHeight+Actar_NS::DistInterCsI);
                if(k>3) PosY= 0.5*Actar_NS::CsIHeight+0.5*Actar_NS::DistInterCsI+(k-4)*(Actar_NS::CsIHeight+Actar_NS::DistInterCsI);

                if(p<5) PosX= 0.5*Actar_NS::CsIWidth+0.5*Actar_NS::DistInterCsI+(4-p)*(Actar_NS::CsIWidth+Actar_NS::DistInterCsI);
                if(p>4) PosX= -0.5*Actar_NS::CsIWidth-0.5*Actar_NS::DistInterCsI+(5-p)*(Actar_NS::CsIWidth+Actar_NS::DistInterCsI);

                G4ThreeVector positionCsI = G4ThreeVector(PosX, PosY, Actar_NS::CsI_PosZ);
                new G4PVPlacement(new G4RotationMatrix(0,0,0),
                                  positionCsI,
                                  m_LogicCsICrystal,"CsI",
                                  m_SquareDetector,false,CsINumber);
                CsINumber++;
            }
        }

        // Set Si sensible
        m_LogicCsICrystal->SetSensitiveDetector(m_CsIScorer);

        // Visualisation of ThinSi
        m_LogicCsICrystal->SetVisAttributes(m_CsIVisAtt) ;
    }

    return m_SquareDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Actar::ReadConfiguration(NPL::InputParser parser){
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Actar");
    if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    vector<string> cart = {"POS","Shape","GasMaterial","GasFraction","Temperature","Pressure","Si","VamosSi","CsI","BeamDump"};
    //vector<string> sphe = {"R","Theta","Phi","Shape","GasMaterial","GasFraction","Temperature","Pressure","Si","CsI","BeamDump"};

    for(unsigned int i = 0 ; i < blocks.size() ; i++){
        if(blocks[i]->HasTokenList(cart)){
            if(NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  Actar " << i+1 <<  endl;

            G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
            string Shape = blocks[i]->GetString("Shape");
            vector<string> GasName = blocks[i]->GetVectorString("GasMaterial");
            vector<int> GasFraction = blocks[i]->GetVectorInt("GasFraction");
            for(unsigned int j=0; j<GasName.size(); j++){
                m_GasMaterial.push_back(GasName[j]);
                m_GasFraction.push_back(GasFraction[j]);
            }
            m_Temperature = blocks[i]->GetDouble("Temperature","kelvin");
            m_Pressure = blocks[i]->GetDouble("Pressure","bar");
            m_build_Silicon =  blocks[i]->GetInt("Si");
            m_build_Vamos_Silicon = blocks[i]->GetInt("VamosSi");
            m_build_CsI     =  blocks[i]->GetInt("CsI");
            m_build_BeamDump = blocks[i]->GetInt("BeamDump");

            AddDetector(Pos,Shape);
        }
        /*else if(blocks[i]->HasTokenList(sphe)){
            if(NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  Actar " << i+1 <<  endl;
            double R = blocks[i]->GetDouble("R","mm");
            double Theta = blocks[i]->GetDouble("Theta","deg");
            double Phi = blocks[i]->GetDouble("Phi","deg");
            string Shape = blocks[i]->GetString("Shape");
            vector<string> GasName = blocks[i]->GetVectorString("GasMaterial");
            vector<int> GasFraction = blocks[i]->GetVectorInt("GasFraction");
            for(unsigned int j=0; j<GasName.size(); j++){
                m_GasMaterial.push_back(GasName[j]);
                m_GasFraction.push_back(GasFraction[j]);
            }
            m_Temperature = blocks[i]->GetDouble("Temperature","kelvin");
            m_Pressure = blocks[i]->GetDouble("Pressure","bar");
            m_build_Silicon =  blocks[i]->GetInt("Si");
            m_build_CsI     =  blocks[i]->GetInt("CsI");

            AddDetector(R,Theta,Phi,Shape);
        }*/
        else{
            cout << "ERROR: check your input file formatting " << endl;
            exit(1);
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Actar::ConstructDetector(G4LogicalVolume* world){
    for (unsigned short i = 0 ; i < m_R.size() ; i++) {

        G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
        G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
        G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
        G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
        // So the face of the detector is at R instead of the middle
        Det_pos+=Det_pos.unit()*Actar_NS::ChamberThickness*0.5;
        // Building Detector reference frame
        G4double ii = cos(m_Theta[i]) * cos(m_Phi[i]);
        G4double jj = cos(m_Theta[i]) * sin(m_Phi[i]);
        G4double kk = -sin(m_Theta[i]);
        G4ThreeVector Y(ii,jj,kk);
        G4ThreeVector w = Det_pos.unit();
        G4ThreeVector u = w.cross(Y);
        G4ThreeVector v = w.cross(u);
        v = v.unit();
        u = u.unit();

        G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);

        new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
                          BuildDetector(),
                          "Actar",world,false,i+1);
    }
    if(!m_ReactionRegion){
        G4ProductionCuts* ecut = new G4ProductionCuts();
        ecut->SetProductionCut(1000,"e-");

        m_ReactionRegion= new G4Region("NPSimulationProcess");
        m_ReactionRegion->SetProductionCuts(ecut);
        m_ReactionRegion->AddRootLogicalVolume(m_logicGas);
        m_ReactionRegion->SetUserLimits(new G4UserLimits(1.2*mm));

        G4Region* Region_cut = new G4Region("RegionCut");
        Region_cut->SetProductionCuts(ecut);
        Region_cut->AddRootLogicalVolume(m_SquareDetector);
    }
    G4FastSimulationManager* mng = m_ReactionRegion->GetFastSimulationManager();
    unsigned int size = m_ReactionModel.size();
    for(unsigned int i = 0 ; i < size ; i++){
        mng->RemoveFastSimulationModel(m_ReactionModel[i]);
    }
    m_ReactionModel.clear();
    G4VFastSimulationModel* fsm;
    fsm = new NPS::BeamReaction("BeamReaction",m_ReactionRegion);
    m_ReactionModel.push_back(fsm);
    fsm = new NPS::Decay("Decay",m_ReactionRegion);
    m_ReactionModel.push_back(fsm);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Actar::InitializeRootOutput(){
    RootOutput *pAnalysis = RootOutput::getInstance();
    TTree *pTree = pAnalysis->GetTree();
    if(!pTree->FindBranch("data")){
        //pTree->Branch("Actar", "TActarData", &m_Event) ;
        pTree->Branch("data", "MEventReduced", &m_EventReduced) ;
    }
    //pTree->SetBranchAddress("Actar", &m_Event) ;
    pTree->SetBranchAddress("data", &m_EventReduced) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Actar::ReadSensitive(const G4Event*){
    m_EventReduced->CoboAsad.clear();
    static ReducedData DataReduced;


    //////////////
    // Pad scorer
    TPCScorers::PS_TPCCathode* PadScorer= (TPCScorers::PS_TPCCathode*) m_ActarScorer->GetPrimitive(0);
    unsigned int size = PadScorer->GetMult();
    // Loop on the Pad map
    //TH1D* h = new TH1D("h","h",25000,0,25000);
    for (unsigned int i = 0 ; i < size ; i++){
        DataReduced.clear();
        double Count = RandGauss::shoot(PadScorer->GetCharge(i),Actar_NS::ResoCharge*PadScorer->GetCharge(i));
        double Time =  PadScorer->GetTime(i);//RandGauss::shoot(Info[1],Actar_NS::ResoTime);
        //cout << "Time= " << Time << endl;
        //int iTime = ((int) Time*20/512)+1;
        //int PadNbr = Info[7];
        int Pad_X = PadScorer->GetPadX(i);//m_PadToXRow[PadNbr];
        int Pad_Y = PadScorer->GetPadY(i);//m_PadToZColumn[PadNbr];

        int co = m_PadToCobo[Pad_X][Pad_Y];
        int as = m_PadToAsad[Pad_X][Pad_Y];
        int ag = m_PadToAGET[Pad_X][Pad_Y];
        int ch = m_PadToChannel[Pad_X][Pad_Y];
        
        if(Count>Actar_NS::ChargeThreshold){
            DataReduced.globalchannelid = ch+(ag<<7)+(as<<9)+(co<<11);
            DataReduced.peakheight.push_back(Count);
            DataReduced.peaktime.push_back(Time);
        }
        m_EventReduced->CoboAsad.push_back(DataReduced);
    }
    
    /*vector<double> Q, T;
     for(int i=0; i<h->GetNbinsX(); i++){
     double count = h->GetBinContent(i);
     double time = h->GetBinCenter(i);
     if(count){
     Q.push_back(count);
     T.push_back(time+500);

     }
     }
     // clear map for next event
     SimulateDigitizer(Q,T,1.40*microsecond,0,8750,25,5);
     delete h;*/

    // Silicon //
    if(m_build_Silicon){
      DSSDScorers::PS_Rectangle* SiScorer= (DSSDScorers::PS_Rectangle*) m_SiliconScorer->GetPrimitive(0);
      unsigned int sizeSi = SiScorer->GetLengthMult();
      // Loop on the ThinSi map
      for(unsigned int i = 0 ; i < sizeSi ; i++){
            DataReduced.clear();
            double E_Si = RandGauss::shoot(SiScorer->GetEnergyLength(i),Actar_NS::ResoSilicon);

            int co = 31;
            int as = 0;
            int ag = 0;
            int ch = SiScorer->GetDetectorLength(i);

            if(E_Si>Actar_NS::EnergyThreshold){
                DataReduced.globalchannelid = ch+(ag<<7)+(as<<9)+(co<<11);
                DataReduced.peaktime.push_back(ch);
                DataReduced.peakheight.push_back(E_Si);
            }
            m_EventReduced->CoboAsad.push_back(DataReduced);
        }
        // Clear Map for next event
        SiScorer->clear();
    }

    // CsI //
    if(m_build_CsI){
       CalorimeterScorers::PS_Calorimeter* CsIScorer= (CalorimeterScorers::PS_Calorimeter*) m_CsIScorer->GetPrimitive(0);
      unsigned int sizeCsI = CsIScorer->GetMult();
      // Loop on the ThinSi map
      for(unsigned int i = 0 ; i < sizeCsI ; i++){
            DataReduced.clear();
            vector<unsigned int> level = CsIScorer->GetLevel(i);
            //double E_CsI = RandGauss::shoot(CsIScorer->GetEnergy(i),Actar_NS::ResoCsI);

            /*if(E_CsI>Actar_NS::EnergyThreshold){
             m_Event->SetCsIEnergy(E_CsI);
             m_Event->SetCsICrystalNumber(level[0]);
             }
             m_EventReduced->CoboAsad.push_back(DataReduced);*/
        }
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*void Actar::SimulateDigitizer(vector<double> E, vector<double> T, double fallTime,double start,double stop, double step,double noise){

 static string formula;
 formula= "";
 static string Es,Ts,var,cond;
 static string fall;
 fall=std::to_string(fallTime);

 for(unsigned int i = 0 ; i < E.size() ; i++){
 if(E[i]!=0 && T[i]!=0){
 Es = std::to_string(E[i]);
 Ts = std::to_string(T[i]);
 cond = ")*(x>"+Ts+")+";
 var = "(x-"+Ts+")";
 formula += Es+"*-1*exp(-"+var+"/"+fall+cond;
 }
 }
 formula+="0";
 //cout << formula << endl;
 TF1* f = new TF1("f",formula.c_str(),start,stop);
 unsigned int size = (stop-start)/step;
 for(unsigned int i = 0 ; i < size ; i++){
 double time = start+i*step;
 double energy = f->Eval(time)+noise*(1-2*G4UniformRand());
 }

 delete f;
 }*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void Actar::InitializeScorers() {
    // This check is necessary in case the geometry is reloaded

    bool already_exist = false;
    vector<G4int> NestingLevel;
    NestingLevel.push_back(0);
    NestingLevel.push_back(2);

    m_ActarScorer   = CheckScorer("ActarScorer",already_exist) ;
    m_SiliconScorer = CheckScorer("SiliconScorer",already_exist);
    //m_VamosSiliconScorer = CheckScorer("VamosSiliconScorer",already_exist);
    m_CsIScorer     = CheckScorer("CsIScorer",already_exist);

    if(already_exist) return;

    G4VPrimitiveScorer* SiScorer = new DSSDScorers::PS_Rectangle("SiliconScorer",0,Actar_NS::SiliconHeight,Actar_NS::SiliconWidth,1,1);
    m_SiliconScorer->RegisterPrimitive(SiScorer);

    /*G4VPrimitiveScorer* VamosSiScorer = new SILICONSCORERS::PS_Silicon_Rectangle("VamosSiliconScorer",0,Actar_NS::VamosSiliconHeight,Actar_NS::VamosSiliconWidth,1,1);
    m_VamosSiliconScorer->RegisterPrimitive(VamosSiScorer);*/

    G4VPrimitiveScorer* CsIScorer= new CalorimeterScorers::PS_Calorimeter("CsI",NestingLevel);
    m_CsIScorer->RegisterPrimitive(CsIScorer);

    vector<int> level; level.push_back(0);
    G4VPrimitiveScorer* Actar_dig= new TPCScorers::PS_TPCCathode("Actar_dig",0) ;
    m_ActarScorer->RegisterPrimitive(Actar_dig);

    G4SDManager::GetSDMpointer()->AddNewDetector(m_ActarScorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_SiliconScorer);
    //G4SDManager::GetSDMpointer()->AddNewDetector(m_VamosSiliconScorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_CsIScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Actar::Construct(){
    return  (NPS::VDetector*) new Actar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
    class proxy_nps_Actar{
    public:
        proxy_nps_Actar(){
            NPS::DetectorFactory::getInstance()->AddToken("Actar","Actar");
            NPS::DetectorFactory::getInstance()->AddDetector("Actar",Actar::Construct);
        }
    };

    proxy_nps_Actar p_nps_Actar;
}
