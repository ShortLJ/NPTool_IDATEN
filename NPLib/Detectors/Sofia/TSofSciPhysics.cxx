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
 *  This class hold SofSci Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSofSciPhysics.h"

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
#include "NPPhysicalConstants.h"
#include "NPGlobalSystemOfUnits.h"

//   ROOT
#include "TChain.h"

ClassImp(TSofSciPhysics)


  ///////////////////////////////////////////////////////////////////////////
TSofSciPhysics::TSofSciPhysics()
  : m_EventData(new TSofSciData),
  m_PreTreatedData(new TSofSciData),
  m_EventPhysics(this),
  m_DET1_PosNs_Min(-20),
  m_DET1_PosNs_Max(20),
  m_DET2_PosNs_Min(-20),
  m_DET2_PosNs_Max(20),
  m_RawTof_Min(0),
  m_RawTof_Max(2000),
  m_NumberOfSignals(3),
  m_NumberOfPmts(2),
  m_NumberOfDetectors(0){ 
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSofSciPhysics::AddDetector(TVector3 ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::BuildPhysicalEvent() {
  Clear();
  // apply thresholds and calibration
  //PreTreat();

  const int N = m_NumberOfDetectors;

  vector<double> S2_pmtL;
  vector<double> S2_pmtR;
  vector<double> S2_pmtTref;
  vector<double> S2_pos;
  vector<double> S2_time;

  vector<double> CC_pmtL;
  vector<double> CC_pmtR;
  vector<double> CC_pmtTref;
  vector<double> CC_pos;
  vector<double> CC_time;

  unsigned int mysizeE = m_EventData->GetMultiplicity();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    int det = m_EventData->GetDetectorNbr(e);
    int pmt = m_EventData->GetPmt(e);
    int FT  = m_EventData->GetFineTime(e);
    int CT  = m_EventData->GetCoarseTime(e);

    double T = CalculateTimeNs(det, pmt, FT, CT);

    if(m_NumberOfDetectors==2){
      if(det==1){
        if(pmt==1) S2_pmtR.push_back(T);
        else if(pmt==2) S2_pmtL.push_back(T);
        else if(pmt==3) S2_pmtTref.push_back(T);
      }

      else if(det==2){
        if(pmt==1) CC_pmtR.push_back(T);
        else if(pmt==2) CC_pmtL.push_back(T);
        else if(pmt==3) CC_pmtTref.push_back(T);
      }
    }

    else if(m_NumberOfDetectors==1){
      if(pmt==1) CC_pmtR.push_back(T);
      else if(pmt==2) CC_pmtL.push_back(T);
      else if(pmt==3) CC_pmtTref.push_back(T);
    }
  }

  if(m_NumberOfDetectors==2){
    double CC_rawpos;
    double S2_rawpos;
    double CC_calpos;
    double S2_calpos;
    double CC_rawtime;
    double S2_rawtime;
    double rawtof;

    multS2_R = S2_pmtR.size();
    multCC_R = CC_pmtR.size();
    multS2_L = S2_pmtL.size();
    multCC_L = CC_pmtL.size();


    static CalibrationManager* Cal = CalibrationManager::getInstance();
    if(S2_pmtTref.size()==1 && CC_pmtTref.size()==1){
      for(unsigned int i=0; i<CC_pmtR.size(); i++){
        for(unsigned int j=0; j<CC_pmtL.size(); j++){
          CC_rawpos = CC_pmtR[i] - CC_pmtL[j];
          CC_calpos = Cal->ApplyCalibration("SofSci/DET2_POSPAR", CC_rawpos);         
          CC_rawtime = 0.5*(CC_pmtR[i] + CC_pmtL[j]);
          if(CC_rawpos<m_DET2_PosNs_Min || CC_rawpos>m_DET2_PosNs_Max)
            continue;
          for(int p=0; p<S2_pmtR.size(); p++){
            for(int k=0; k<S2_pmtL.size(); k++){
              S2_rawpos = S2_pmtR[p] - S2_pmtL[k];
              S2_calpos = Cal->ApplyCalibration("SofSci/DET1_POSPAR", S2_rawpos);         
              S2_rawtime = 0.5*(S2_pmtR[p] + S2_pmtL[k]);

              if(S2_rawpos<m_DET1_PosNs_Min || S2_rawpos>m_DET1_PosNs_Max)
                continue;

              /*CC_pos.push_back(CC_rawpos);
                CC_time.push_back(CC_rawtime);
                S2_pos.push_back(S2_rawpos); 
                S2_time.push_back(S2_rawtime); */

              rawtof = CC_rawtime - CC_pmtTref[0] - (S2_rawtime - S2_pmtTref[0]);
              
              if(rawtof<m_RawTof_Min || rawtof>m_RawTof_Max)
                continue;

              double velocity_mns;
              double caltof;
              double betaS2;

              /*cout << "*** Printing physics calibration parameter:" << endl;
              cout << Cal->GetValue("SofSci/TOF2INV_V",0) << endl;
              cout << Cal->GetValue("SofSci/TOF2INV_V",1) << endl;
              cout << Cal->GetValue("SofSci/LENGTH_S2",0) << endl;*/

              velocity_mns = 1./Cal->ApplyCalibration("SofSci/TOF2INV_V",rawtof);
              caltof       = Cal->GetValue("SofSci/LENGTH_S2",0) / velocity_mns;
              betaS2       = velocity_mns * m/ns / NPUNITS::c_light;
              
              // Filling ouput Tree;
              DetectorNbr.push_back(1);
              TimeNs.push_back(S2_rawtime);
              PosNs.push_back(S2_rawpos);
              PosMm.push_back(S2_calpos);

              DetectorNbr.push_back(2);
              TimeNs.push_back(CC_rawtime);
              PosNs.push_back(CC_rawpos);
              PosMm.push_back(CC_calpos);

              RawTof.push_back(rawtof);
              CalTof.push_back(caltof);
              VelocityMNs.push_back(velocity_mns);
              Beta.push_back(betaS2);
            }
          }
        }
      }
    }
  }
  else if(m_NumberOfDetectors==1){
    double CC_rawpos;
    double CC_rawtime;
    if(CC_pmtTref.size()==1){
      for(unsigned int i=0; i<CC_pmtR.size(); i++){
        for(unsigned int j=0; j<CC_pmtL.size(); j++){
          CC_rawpos = CC_pmtR[i] - CC_pmtL[j];
          CC_rawtime = 0.5*(CC_pmtR[i] + CC_pmtL[j]);

          if(CC_rawpos<m_DET2_PosNs_Min || CC_rawpos>m_DET2_PosNs_Max)
            continue;

          //Filling output Tree;
          DetectorNbr.push_back(1);
          TimeNs.push_back(CC_rawtime);
          PosNs.push_back(CC_rawpos);
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
double TSofSciPhysics::CalculateTimeNs(int det, int pmt, int ft, int ct){

  static CalibrationManager* Cal = CalibrationManager::getInstance();
  //ft = ft+1; 
  double par = Cal->GetValue("SofSci/DET"+NPL::itoa(det)+"_SIGNAL"+NPL::itoa(pmt)+"_TIME",ft);
  double r = (double)rand.Rndm()-0.5;
  double ClockOffset = Cal->GetValue("SofSci/DET"+NPL::itoa(det)+"_SIGNAL"+NPL::itoa(pmt)+"_CLOCKOFFSET",0);
  double ict_ns = ((double)ct - ClockOffset) * 5.; 
  double ift_ns;

  if(r<0){
    double par_prev = Cal->GetValue("SofSci/DET"+NPL::itoa(det)+"_SIGNAL"+NPL::itoa(pmt)+"_TIME",ft-1);
    ift_ns = par + r*(par - par_prev);
  }

  else{
    double par_next = Cal->GetValue("SofSci/DET"+NPL::itoa(det)+"_SIGNAL"+NPL::itoa(pmt)+"_TIME",ft+1);
    ift_ns = par + r*(par_next - par);
  }

  double time_ns = (double)ict_ns - ift_ns;

  return time_ns;
}

///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysize = m_EventData->GetMultiplicity();
  for (unsigned int i = 0; i < mysize ; ++i) {

    m_PreTreatedData->SetDetectorNbr(m_EventData->GetDetectorNbr(i));
    m_PreTreatedData->SetPmt(m_EventData->GetPmt(i));
    m_PreTreatedData->SetCoarseTime(m_EventData->GetCoarseTime(i));
    m_PreTreatedData->SetFineTime(m_EventData->GetFineTime(i));
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSofSci.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSofSci.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSofSci.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSofSci.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSofSci";
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

      else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
      }
      else if (whatToDo=="DET1_POSNS_MIN") {
        AnalysisConfigFile >> DataBuffer;
        m_DET1_PosNs_Min = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_DET1_PosNs_Min << endl;
      }
      else if (whatToDo=="DET1_POSNS_MAX") {
        AnalysisConfigFile >> DataBuffer;
        m_DET1_PosNs_Max = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_DET1_PosNs_Max << endl;
      }
      else if (whatToDo=="DET2_POSNS_MIN") {
        AnalysisConfigFile >> DataBuffer;
        m_DET2_PosNs_Min = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_DET2_PosNs_Min << endl;
      }
      else if (whatToDo=="DET2_POSNS_MAX") {
        AnalysisConfigFile >> DataBuffer;
        m_DET2_PosNs_Max = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_DET2_PosNs_Max << endl;
      }
      else if (whatToDo=="RAWTOF_MIN") {
        AnalysisConfigFile >> DataBuffer;
        m_RawTof_Min = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_RawTof_Min << endl;
      }
      else if (whatToDo=="RAWTOF_MAX") {
        AnalysisConfigFile >> DataBuffer;
        m_RawTof_Max = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_RawTof_Max << endl;
      }



      else {
        ReadingStatus = false;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::Clear() {
  DetectorNbr.clear();
  TimeNs.clear();
  PosNs.clear();
  PosMm.clear();
  RawTof.clear();
  CalTof.clear();
  VelocityMNs.clear();
  Beta.clear();

  multS2_R=0;
  multS2_L=0;
  multCC_R=0;
  multCC_L=0;
}



///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofSci");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofSci " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      m_X = Pos.Z();
      m_Y = Pos.Z();
      m_Z = Pos.Z();
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofSci " << i+1 <<  endl;
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

  ReadAnalysisConfig();
}


///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  Cal->AddParameter("SofSci","TOF2INV_V","SofSci_TOF2INV_V");
  Cal->AddParameter("SofSci","LENGTH_S2","SofSci_LENGTH_S2");
  
  for(int d = 0; d < m_NumberOfDetectors; d++){
    Cal->AddParameter("SofSci","DET"+NPL::itoa(d+1)+"_POSPAR","SofSci_DET"+NPL::itoa(d+1)+"_POSPAR");
    for(int s = 0; s < m_NumberOfSignals; s++){
      Cal->AddParameter("SofSci","DET"+NPL::itoa(d+1)+"_SIGNAL"+NPL::itoa(s+1)+"_TIME","SofSci_DET"+NPL::itoa(d+1)+"_SIGNAL"+NPL::itoa(s+1)+"_TIME");
      Cal->AddParameter("SofSci","DET"+NPL::itoa(d+1)+"_SIGNAL"+NPL::itoa(s+1)+"_CLOCKOFFSET","SofSci_DET"+NPL::itoa(d+1)+"_SIGNAL"+NPL::itoa(s+1)+"_CLOCKOFFSET");
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SofSci",  true );
  inputChain->SetBranchAddress("SofSci", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SofSci", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSofSciPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("SofSci", "TSofSciPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSofSciPhysics::Construct() {
  return (NPL::VDetector*) new TSofSciPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_SofSci{
    public:
      proxy_SofSci(){
        NPL::DetectorFactory::getInstance()->AddToken("SofSci","SofSci");
        NPL::DetectorFactory::getInstance()->AddDetector("SofSci",TSofSciPhysics::Construct);
      }
  };

  proxy_SofSci p_SofSci;
}

