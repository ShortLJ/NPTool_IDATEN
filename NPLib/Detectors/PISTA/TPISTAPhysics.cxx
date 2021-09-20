/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TPISTAPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

//   ROOT
#include "TChain.h"

ClassImp(TPISTAPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TPISTAPhysics::TPISTAPhysics(){
    EventMultiplicity = 0;
    m_EventData = new TPISTAData;
    m_PreTreatedData = new TPISTAData;
    m_EventPhysics = this;
    m_Spectra = NULL;
    m_E_RAW_Threshold = 0; // adc channels
    m_E_Threshold = 0;     // MeV
    m_NumberOfDetectors = 0;
    m_MaximumStripMultiplicityAllowed = 10;
    m_StripEnergyMatching = 0.050;
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TPISTAPhysics::AddDetector(TVector3){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::AddDetector(double R, double Theta, double Phi){
  m_NumberOfDetectors++;

  double Height = 61.8; // mm
  double Base = 78.1; // mm
  double NumberOfStripsX = 62;
  double NumberOfStripsY = 97;
  
  double StripPitchHeight = Height / NumberOfStripsY; // mm
  double StripPitchBase = Base / NumberOfStripsX; // mm


  // Vector U on detector face (parallel to Y strips) Y strips are along X axis
  TVector3 U;
  // Vector V on detector face (parallel to X strips)
  TVector3 V;
  // Vector W normal to detector face (pointing to the back)
  TVector3 W;
  // Vector C position of detector face center
  TVector3 C;
  C = TVector3(R*sin(Theta)*cos(Phi),
        R*sin(Theta)*sin(Phi),
        Height*0.5+R*cos(Theta));

  TVector3 P = TVector3(cos(Theta)*cos(Phi),
      cos(Theta)*sin(Phi),
      -sin(Theta));

  W = C.Unit();
  U = W.Cross(P);
  V = W.Cross(U);

  U = U.Unit();
  V = V.Unit();

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  double X, Y, Z;

  // Moving C to the 1.1 Corner;
  TVector3 Strip_1_1;
  Strip_1_1 = C - (0.5*Base*U + 0.5*Height*V) + U*(StripPitchBase / 2.) + V*(StripPitchHeight / 2.);

  TVector3 StripPos;
  for(int i=0; i<NumberOfStripsX; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<NumberOfStripsY; j++){
      StripPos = Strip_1_1 + i*U*StripPitchBase + j*V*StripPitchHeight;
      lineX.push_back(StripPos.X());
      lineY.push_back(StripPos.Y());
      lineZ.push_back(StripPos.Z());
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_StripPositionX.push_back(OneDetectorStripPositionX);
  m_StripPositionY.push_back(OneDetectorStripPositionY);
  m_StripPositionZ.push_back(OneDetectorStripPositionZ);
} 

///////////////////////////////////////////////////////////////////////////
TVector3 TPISTAPhysics::GetPositionOfInteraction(const int i){
  TVector3 Position = TVector3(GetStripPositionX(DetectorNumber[i], E_StripX[i], DE_StripY[i]),
      GetStripPositionY(DetectorNumber[i], E_StripX[i], DE_StripY[i]),
      GetStripPositionZ(DetectorNumber[i], E_StripX[i], DE_StripY[i]));
  
  /*TVector3 Position = TVector3(GetStripPositionX(DetectorNumber[i], DE_StripX[i], E_StripY[i]),
      GetStripPositionY(DetectorNumber[i], DE_StripX[i], E_StripY[i]),
      GetStripPositionZ(DetectorNumber[i], DE_StripX[i], E_StripY[i]));
*/

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TPISTAPhysics::GetDetectorNormal(const int i){
  TVector3 U = TVector3(GetStripPositionX(DetectorNumber[i],62,1),
      GetStripPositionY(DetectorNumber[i],62,1),
      GetStripPositionZ(DetectorNumber[i],62,1))

    -TVector3(GetStripPositionX(DetectorNumber[i],62,1),
      GetStripPositionY(DetectorNumber[i],62,1),
      GetStripPositionZ(DetectorNumber[i],62,1));

  TVector3 V = TVector3(GetStripPositionX(DetectorNumber[i],62,97),
      GetStripPositionY(DetectorNumber[i],62,97),
      GetStripPositionZ(DetectorNumber[i],62,97))

    -TVector3(GetStripPositionX(DetectorNumber[i],62,1),
      GetStripPositionY(DetectorNumber[i],62,1),
      GetStripPositionZ(DetectorNumber[i],62,1));

  TVector3 Normal = U.Cross(V);

  return (Normal.Unit());
}
///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  if(1 /*CheckEvent() == 1*/){
    //vector<TVector2> couple = Match_X_Y();
    //EventMultiplicity = couple.size();

    int FirstStageMult = m_PreTreatedData->GetFirstStageMultXEnergy();
    int SecondStageMult = m_PreTreatedData->GetSecondStageMultXEnergy();

    for(unsigned int i=0; i<FirstStageMult; i++){
      for(unsigned int j=0; j<SecondStageMult; j++){
        int DE_N = m_PreTreatedData->GetFirstStage_XE_DetectorNbr(i);
        int E_N = m_PreTreatedData->GetSecondStage_XE_DetectorNbr(j);
        if(DE_N==E_N){

          int DE_X = m_PreTreatedData->GetFirstStage_XE_StripNbr(i);
          int DE_Y = m_PreTreatedData->GetFirstStage_YE_StripNbr(i);
          double DE_XE = m_PreTreatedData->GetFirstStage_XE_Energy(i);
          double DE_YE = m_PreTreatedData->GetFirstStage_YE_Energy(i);

          int E_X = m_PreTreatedData->GetSecondStage_XE_StripNbr(j);
          int E_Y = m_PreTreatedData->GetSecondStage_YE_StripNbr(j);
          double E_XE = m_PreTreatedData->GetSecondStage_XE_Energy(j);
          double E_YE = m_PreTreatedData->GetSecondStage_YE_Energy(j);

          DetectorNumber.push_back(DE_N);
          DE_StripX.push_back(DE_X);
          DE_StripY.push_back(DE_Y);
          //DE.push_back(DE_YE);
          DE.push_back(DE_XE);
          E_StripX.push_back(E_X);
          E_StripY.push_back(E_Y);
          E.push_back(E_XE);

          PosX.push_back(GetPositionOfInteraction(i).x());
          PosY.push_back(GetPositionOfInteraction(i).y());
          PosZ.push_back(GetPositionOfInteraction(i).z());

        }
      }
    }
  }
  EventMultiplicity = DetectorNumber.size();
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TPISTAPhysics::Match_X_Y(){
  vector<TVector2> ArrayOfGoodCouple;

  static unsigned int m_XEMult, m_YEMult;
  m_XEMult = m_PreTreatedData->GetFirstStageMultXEnergy();
  m_YEMult = m_PreTreatedData->GetFirstStageMultYEnergy();

  if(m_XEMult>m_MaximumStripMultiplicityAllowed || m_YEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_XEMult; i++){
    for(unsigned int j=0; j<m_YEMult; j++){

      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetFirstStage_XE_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetFirstStage_YE_DetectorNbr(j);

      // if same detector check energy
      if(XDetNbr == YDetNbr){
        // Declaration of variable for clarity
        double XE = m_PreTreatedData->GetFirstStage_XE_Energy(i);
        double YE = m_PreTreatedData->GetFirstStage_YE_Energy(i);
        double XStripNbr = m_PreTreatedData->GetFirstStage_XE_StripNbr(i);
        double YStripNbr = m_PreTreatedData->GetFirstStage_YE_StripNbr(i);

        // look if energy matches
        if(abs(XE-YE)/2.<m_StripEnergyMatching){
          ArrayOfGoodCouple.push_back(TVector2(i,j));
        }
      }
    }
  }

  return ArrayOfGoodCouple;
}

///////////////////////////////////////////////////////////////////////////
int TPISTAPhysics::CheckEvent(){
  // Check the size of the different elements
  if(m_PreTreatedData->GetFirstStageMultXEnergy() == m_PreTreatedData->GetFirstStageMultYEnergy() )
    return 1;

  else
    return -1;
}


///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  //////
  // First Stage Energy
  unsigned int sizeFront = m_EventData->GetFirstStageMultXEnergy();
  for (UShort_t i = 0; i < sizeFront ; ++i) {
    if (m_EventData->GetFirstStage_XE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetFirstStage_XE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("PISTA/ENERGY"+NPL::itoa(m_EventData->GetFirstStage_XE_DetectorNbr(i)),m_EventData->GetFirstStage_XE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetFirstStageXE(m_EventData->GetFirstStage_XE_DetectorNbr(i), m_EventData->GetFirstStage_XE_StripNbr(i), Energy);
      }
    }
  }
  unsigned int sizeBack = m_EventData->GetFirstStageMultYEnergy();
  for (UShort_t i = 0; i < sizeBack ; ++i) {
    if (m_EventData->GetFirstStage_YE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetFirstStage_YE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("PISTA/ENERGY"+NPL::itoa(m_EventData->GetFirstStage_YE_DetectorNbr(i)),m_EventData->GetFirstStage_YE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetFirstStageYE(m_EventData->GetFirstStage_YE_DetectorNbr(i), m_EventData->GetFirstStage_YE_StripNbr(i), Energy);
      }
    }
  }
  // First Stage Time 
  unsigned int mysize = m_EventData->GetFirstStageMultXTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("PISTA/TIME"+NPL::itoa(m_EventData->GetFirstStage_XT_DetectorNbr(i)),m_EventData->GetFirstStage_XT_Time(i));
    m_PreTreatedData->SetFirstStageXT(m_EventData->GetFirstStage_XT_DetectorNbr(i), m_EventData->GetFirstStage_XT_StripNbr(i), Time);
  }

  //////
  // Second Stage Energy
  sizeFront = m_EventData->GetSecondStageMultXEnergy();
  for (UShort_t i = 0; i < sizeFront ; ++i) {
    if (m_EventData->GetSecondStage_XE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetSecondStage_XE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("PISTA/ENERGY"+NPL::itoa(m_EventData->GetSecondStage_XE_DetectorNbr(i)),m_EventData->GetSecondStage_XE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetSecondStageXE(m_EventData->GetSecondStage_XE_DetectorNbr(i), m_EventData->GetSecondStage_XE_StripNbr(i), Energy);
      }
    }
  }
  sizeBack = m_EventData->GetSecondStageMultYEnergy();
  for (UShort_t i = 0; i < sizeBack ; ++i) {
    if (m_EventData->GetSecondStage_YE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetSecondStage_YE_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("PISTA/ENERGY"+NPL::itoa(m_EventData->GetSecondStage_YE_DetectorNbr(i)),m_EventData->GetSecondStage_YE_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetSecondStageYE(m_EventData->GetSecondStage_YE_DetectorNbr(i), m_EventData->GetSecondStage_YE_StripNbr(i), Energy);
      }
    }
  }

}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigPISTA.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigPISTA.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigPISTA.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigPISTA.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigPISTA";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::Clear() {
  EventMultiplicity = 0;

  // Position Information
  PosX.clear();
  PosY.clear();
  PosZ.clear();

  // DSSD
  DetectorNumber.clear();
  E.clear();
  DE_StripX.clear();
  DE_StripY.clear();
  E_StripX.clear();
  E_StripY.clear();
  Time.clear();
  DE.clear();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("PISTA");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");

      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");

      AddDetector(R, Theta, Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitSpectra() {
  m_Spectra = new TPISTASpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TPISTAPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("PISTA", "D"+ NPL::itoa(i+1)+"_ENERGY","PISTA_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("PISTA", "D"+ NPL::itoa(i+1)+"_TIME","PISTA_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("PISTA",  true );
  inputChain->SetBranchAddress("PISTA", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("PISTA", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("PISTA", "TPISTAPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TPISTAPhysics::Construct() {
  return (NPL::VDetector*) new TPISTAPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_PISTA{
    public:
      proxy_PISTA(){
        NPL::DetectorFactory::getInstance()->AddToken("PISTA","PISTA");
        NPL::DetectorFactory::getInstance()->AddDetector("PISTA",TPISTAPhysics::Construct);
      }
  };

  proxy_PISTA p_PISTA;
}

