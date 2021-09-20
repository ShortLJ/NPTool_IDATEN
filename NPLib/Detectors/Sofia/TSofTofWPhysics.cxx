/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : November 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SofTofW Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSofTofWPhysics.h"

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

ClassImp(TSofTofWPhysics)


  ///////////////////////////////////////////////////////////////////////////
TSofTofWPhysics::TSofTofWPhysics()
  : m_EventData(new TSofTofWData),
  m_PreTreatedData(new TSofTofWData),
  m_EventPhysics(this),
  m_E_RAW_Threshold(0), // adc channels
  m_E_Threshold(0),     // MeV
  m_NumberOfPlastics(28),
  m_StartTime(-1),
  m_TofAlignedValue(0), // ns
  m_NumberOfDetectors(0) {
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSofTofWPhysics::AddDetector(TVector3 ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  if(m_StartTime == -1)
    return;

  PreTreat();

  double T1[28][10], T2[28][10];
  int mult1[28], mult2[28];
  for(int i=0; i<28; i++){
    mult1[i] = 0;
    mult2[i] = 0;
    for(int j=0; j<10; j++){
      T1[i][j] = 0;
      T2[i][j] = 0;
    }
  }

  unsigned int mysizeE = m_PreTreatedData->GetMultiplicity();
  for (UShort_t i = 0; i < mysizeE ; i++) {
    int plastic = m_PreTreatedData->GetPlasticNbr(i);
    int pmt     = m_PreTreatedData->GetPmt(i);
    int FT      = m_PreTreatedData->GetFineTime(i);
    int CT      = m_PreTreatedData->GetCoarseTime(i);

    double T = CalculateTimeNs(plastic, pmt, FT, CT);

    if(pmt==1){
      T1[plastic-1][mult1[plastic-1]] = T;
      mult1[plastic-1]++;
    }
    else if(pmt==2){
      T2[plastic-1][mult2[plastic-1]] = T;    
      mult2[plastic-1]++;
    }
  }

  static CalibrationManager* Cal = CalibrationManager::getInstance();
  for(int p=0; p<m_NumberOfPlastics; p++){
    if(mult1[p]==1 && mult2[p]==1){
      double time_ns = 0.5*(T1[p][0] + T2[p][0]);
      double rawpos = T1[p][0] - T2[p][0];
      double calpos = Cal->ApplyCalibration("SofTofW/TOFW"+NPL::itoa(p+1)+"_POSPAR",rawpos);
      double rawtof = time_ns - m_StartTime;
      double caltof = Cal->ApplyCalibration("SofTofW/TOFW"+NPL::itoa(p+1)+"_TOFPAR",rawtof) + m_TofAlignedValue;

      PlasticNbr.push_back(p+1);
      TimeNs.push_back(time_ns);
      RawPosY.push_back(rawpos);
      CalPosY.push_back(calpos);
      RawTof.push_back(rawtof);
      CalTof.push_back(caltof);
    }
  }
  m_StartTime = -1;
}

///////////////////////////////////////////////////////////////////////////
double TSofTofWPhysics::CalculateTimeNs(int det, int pmt, int ft, int ct){

  static CalibrationManager* Cal = CalibrationManager::getInstance();
  double par = Cal->GetValue("SofTofW/TOFW"+NPL::itoa(det)+"_PMT"+NPL::itoa(pmt)+"_TIME",ft);
  double r = (double)rand.Rndm()-0.5;
  double ClockOffset = Cal->GetValue("SofTofW/TOFW"+NPL::itoa(det)+"_PMT"+NPL::itoa(pmt)+"_CLOCKOFFSET",0);
  double ict_ns = ((double)ct - ClockOffset) * 5.; 
  double ift_ns;

  if(r<0){
    double par_prev = Cal->GetValue("SofTofW/TOFW"+NPL::itoa(det)+"_PMT"+NPL::itoa(pmt)+"_TIME",ft-1);
    ift_ns = par + r*(par - par_prev);
  }

  else{
    double par_next = Cal->GetValue("SofTofW/TOFW"+NPL::itoa(det)+"_PMT"+NPL::itoa(pmt)+"_TIME",ft+1);
    ift_ns = par + r*(par_next - par);
  }

  double time_ns = (double)ict_ns - ift_ns;

  return time_ns;
}

///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Energy
  unsigned int mysize = m_EventData->GetMultiplicity();
  for (UShort_t i = 0; i < mysize ; ++i) {
    if(m_EventData->GetWhichFlag(i)==1){
      m_PreTreatedData->SetPlasticNbr(m_EventData->GetPlasticNbr(i));
      m_PreTreatedData->SetPmt(m_EventData->GetPmt(i));
      m_PreTreatedData->SetCoarseTime(m_EventData->GetCoarseTime(i));
      m_PreTreatedData->SetFineTime(m_EventData->GetFineTime(i));
      m_PreTreatedData->SetWhichFlag(m_EventData->GetWhichFlag(i));
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSofTofW.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSofTofW.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSofTofW.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSofTofW.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSofTofW";
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
void TSofTofWPhysics::Clear() {
  PlasticNbr.clear();
  TimeNs.clear();
  RawPosY.clear();
  CalPosY.clear();
  RawTof.clear();
  CalTof.clear();
}



///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofTofW");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTofW " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTofW " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetector(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}


///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfPlastics; ++i) {
    Cal->AddParameter("SofTofW", "TOFW"+ NPL::itoa(i+1)+"_POSPAR","SofTofW_TOFW"+ NPL::itoa(i+1)+"_POSPAR");
    Cal->AddParameter("SofTofW", "TOFW"+ NPL::itoa(i+1)+"_TOFPAR","SofTofW_TOFW"+ NPL::itoa(i+1)+"_TOFPAR");

    for(int j = 0; j < 2; j++){
      Cal->AddParameter("SofTofW", "TOFW"+ NPL::itoa(i+1)+"_PMT"+NPL::itoa(j+1)+"_TIME","SofTofW_TOFW"+ NPL::itoa(i+1)+"_PMT"+NPL::itoa(j+1)+"_TIME");
      Cal->AddParameter("SofTofW", "TOFW"+ NPL::itoa(i+1)+"_PMT"+NPL::itoa(j+1)+"_CLOCKOFFSET","SofTofW_TOFW"+ NPL::itoa(i+1)+"_PMT"+NPL::itoa(j+1)+"_CLOCKOFFSET");
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SofTofW",  true );
  inputChain->SetBranchAddress("SofTofW", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SofTofW", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSofTofWPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("SofTofW", "TSofTofWPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSofTofWPhysics::Construct() {
  return (NPL::VDetector*) new TSofTofWPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_SofTofW{
    public:
      proxy_SofTofW(){
        NPL::DetectorFactory::getInstance()->AddToken("SofTofW","SofTofW");
        NPL::DetectorFactory::getInstance()->AddDetector("SofTofW",TSofTofWPhysics::Construct);
      }
  };

  proxy_SofTofW p_SofTofW;
}

