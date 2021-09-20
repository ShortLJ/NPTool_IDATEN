/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Catana Treated  data                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TCatanaPhysics.h"

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
using namespace NPUNITS;
#include "NPTrackingUtility.h"

//   ROOT
#include "TChain.h"

ClassImp(TCatanaPhysics)


///////////////////////////////////////////////////////////////////////////
TCatanaPhysics::TCatanaPhysics()
  : m_EventData(new TCatanaData),
  m_PreTreatedData(new TCatanaData),
  m_EventPhysics(this),
  m_Spectra(0),
  m_E_RAW_Threshold(0), // adc channels
  m_E_Threshold(0),     // MeV
  m_NumberOfDetectors(0) {
  }


///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::AddDetector(double X, double Y, double Z, double Theta, double Phi, int ID, int Type){
  m_NumberOfDetectors++;
  TVector3 Pos(X,Y,Z);
  m_Position[ID]=Pos+m_Ref;
  m_Theta[ID]=Theta;
  m_Phi[ID]=Phi;
  m_Type[ID]=Type;
}

///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}


///////////////////////////////////////////////////////////////////////////
TVector3 TCatanaPhysics::GetPositionOfInteraction(int& i){
  return m_Position[DetectorNumber[i]];  
}
///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::ReadCSV(string path){
  std::ifstream csv(path); 
  if(!csv.is_open()){
    std::ostringstream message;
    message << "Catana csv file " << path << " not found " << std::endl;
  }

  int ID, type,layer;
  double X,Y,Z,Theta,Phi;
  string buffer;
  // ignore first line
  getline(csv,buffer);
  while(csv >> ID >> buffer >> type >> buffer >> layer >> buffer >> X >> buffer >> Y >> buffer >> Z >> buffer >> Theta >> buffer >> Phi){
    if(type<6)
      AddDetector(X,Y,Z,Theta*deg,Phi*deg,ID,type);
    else{
      // ignore other type for which I don't have the geometry
    }
  }

  return;
}
///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  unsigned int mysizeT = m_PreTreatedData->GetMultTime();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    for (UShort_t t = 0; t < mysizeT ; t++) {
      if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t)) {
        DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
        Energy.push_back(m_PreTreatedData->Get_Energy(e));
        Time.push_back(m_PreTreatedData->Get_Time(t));
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
unsigned int TCatanaPhysics::FindClosestHitToLine(const TVector3& v1, const TVector3& v2,double& d){

  d = 1e32;
  unsigned result = 0; 
  unsigned int size = DetectorNumber.size();
  for(unsigned int i = 0 ; i < size ; i++){
    double current_d = NPL::MinimumDistancePointLine(v1,v2,m_Position[DetectorNumber[i]]) ;
    if(current_d < d){
      d=current_d;
      result=i; 
    }
  }

  if(d==1e32)
    d=-1000;

  return result;
}
///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::PreTreat() {
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
      Double_t Energy = Cal->ApplyCalibration("Catana/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("Catana/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigCatana.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigCatana.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigCatana.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigCatana.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigCatana";
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
void TCatanaPhysics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::ReadConfiguration(NPL::InputParser parser){
  // CSV config
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("Catana","CSV");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " CSV block found " << endl; 

  vector<string> token = {"Path","Pos","Rshift"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      string path = blocks[i]->GetString("Path");
      //double Rshift = blocks[i]->GetDouble("Rshift","micrometer");
      // Reference position of the whole array
      m_Ref = blocks[i]->GetTVector3("Pos","mm");
      ReadCSV(path);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  // Type 1
  blocks = parser.GetAllBlocksWithTokenAndValue("Catana","Detector");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  token = {"X","Y","Z","Theta","Phi","ID","Type"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Catana " << i+1 <<  endl;
      double X = blocks[i]->GetDouble("X","mm");
      double Y = blocks[i]->GetDouble("Y","mm");
      double Z = blocks[i]->GetDouble("Z","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      int    ID  = blocks[i]->GetInt("ID");
      int    Type =  blocks[i]->GetInt("Type"); 
      AddDetector(X,Y,Z,Theta,Phi,ID,Type);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::InitSpectra() {
  m_Spectra = new TCatanaSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TCatanaPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Catana", "D"+ NPL::itoa(i+1)+"_ENERGY","Catana_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("Catana", "D"+ NPL::itoa(i+1)+"_TIME","Catana_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Catana",  true );
  inputChain->SetBranchAddress("Catana", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Catana", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TCatanaPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Catana", "TCatanaPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TCatanaPhysics::Construct() {
  return (NPL::VDetector*) new TCatanaPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_Catana{
    public:
      proxy_Catana(){
        NPL::DetectorFactory::getInstance()->AddToken("Catana","Catana");
        NPL::DetectorFactory::getInstance()->AddDetector("Catana",TCatanaPhysics::Construct);
      }
  };

  proxy_Catana p_Catana;
}

