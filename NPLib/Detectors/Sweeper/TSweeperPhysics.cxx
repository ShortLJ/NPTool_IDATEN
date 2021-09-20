/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: B. Monteagudo  contact address: monteagu@frib.msu.edu                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Sweeper Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSweeperPhysics.h"

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
#include "NPSystemOfUnits.h"

//   ROOT
#include "TChain.h"

ClassImp(TSweeperPhysics)

using namespace NPUNITS;

///////////////////////////////////////////////////////////////////////////
TSweeperPhysics::TSweeperPhysics()
   : m_EventData(new TSweeperData),
     m_PreTreatedData(new TSweeperData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_DriftSpeed(0), //cm/us
     m_DistDC(0), //m
     m_NumberOfDetectors(0){
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSweeperPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,shape);
} 
  
///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // Energy
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  for (UShort_t e = 0; e < mysizeE ; e++) {
       DetectorENumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
       Energy.push_back(m_PreTreatedData->Get_Energy(e));
  }
  //Time
  unsigned int mysizeT = m_PreTreatedData->GetMultTime();
  for (UShort_t t = 0; t < mysizeT ; t++) {
    DetectorTNumber.push_back(m_PreTreatedData->GetT_DetectorNbr(t));
    Time.push_back(m_PreTreatedData->Get_Time(t));
  }
  //Position
  unsigned int mysizeDC = m_PreTreatedData->GetMultPosition();
  for (UShort_t d = 0; d < mysizeDC ; d++) {
    DetectorDCNumber.push_back(m_PreTreatedData->GetDC_DetectorNbr(d));
    Drift_X.push_back(m_PreTreatedData->Get_X(d));
    Drift_Ypos.push_back(m_PreTreatedData->Get_Y(d));
    
    // Calculate Y from Drift Speed
    double DriftTime = m_PreTreatedData->Get_DriftTime(d)*microsecond;
    double Y = DriftTime*m_DriftSpeed;
    double HalfTime = 5*cm/m_DriftSpeed;
    
    if(Y<HalfTime) Drift_Y.push_back(-Y);
    else Drift_Y.push_back(Y);
  }
  //Calculate Theta and Phi
  if(mysizeDC==2) {
    TVector3 r(Drift_X[1]-Drift_X[0],Drift_Y[1]-Drift_Y[0], m_DistDC*m);
    
    Drift_ThetaX.push_back(atan(r.Px()/r.Pz())); 
    Drift_ThetaY.push_back(atan(r.Py()/r.Pz()));
    
  }
}

///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Energy
  unsigned int mysize = m_EventData->GetMultEnergy();
  for (UShort_t i = 0; i < mysize ; ++i) {
    if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = Cal->ApplyCalibration("Sweeper/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("Sweeper/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }
  // Position 
  mysize = m_EventData->GetMultPosition();
  for (UShort_t i = 0; i < mysize; ++i) {
    //Double_t Y= Cal->ApplyCalibration("Sweeper/DriftSpeed"+NPL::itoa(m_EventData->GetDC_DetectorNbr(i)),m_EventData->Get_DriftTime(i));
    m_PreTreatedData->SetDriftPosition(m_EventData->GetDC_DetectorNbr(i),m_EventData->Get_DriftTime(i), m_EventData->Get_X(i) );
    m_PreTreatedData->SetPosition(m_EventData->GetDC_DetectorNbr(i),m_EventData->Get_X(i), m_EventData->Get_Y(i) );
    
  }  
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSweeper.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSweeper.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSweeper.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSweeper.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSweeper";
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
      else if (whatToDo=="DriftSpeed") {
        AnalysisConfigFile >> DataBuffer;
        m_DriftSpeed = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_DriftSpeed << endl;
      }
      else if (whatToDo=="DistCRDCs") {
        AnalysisConfigFile >> DataBuffer;
        m_DistDC = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_DistDC << endl;
      }
   
      else {
        ReadingStatus = false;
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::Clear() {
  DetectorENumber.clear();
  DetectorTNumber.clear();
  DetectorDCNumber.clear();
  Energy.clear();
  Time.clear();
  Drift_X.clear();
  Drift_Y.clear();
  Drift_ThetaX.clear();
  Drift_ThetaY.clear();
  
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Sweeper");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Sweeper " << i+1 <<  endl;
    
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Sweeper " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R,Theta,Phi,Shape);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::InitSpectra() {
  m_Spectra = new TSweeperSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TSweeperPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Sweeper", "D"+ NPL::itoa(i+1)+"_ENERGY","Sweeper_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("Sweeper", "D"+ NPL::itoa(i+1)+"_TIME","Sweeper_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Sweeper",  true );
  inputChain->SetBranchAddress("Sweeper", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Sweeper", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSweeperPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Sweeper", "TSweeperPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSweeperPhysics::Construct() {
  return (NPL::VDetector*) new TSweeperPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Sweeper{
  public:
    proxy_Sweeper(){
      NPL::DetectorFactory::getInstance()->AddToken("Sweeper","Sweeper");
      NPL::DetectorFactory::getInstance()->AddDetector("Sweeper",TSweeperPhysics::Construct);
    }
};

proxy_Sweeper p_Sweeper;
}

