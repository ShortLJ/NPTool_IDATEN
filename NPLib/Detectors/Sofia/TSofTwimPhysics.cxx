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
 *  This class hold SofTwim Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSofTwimPhysics.h"

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

ClassImp(TSofTwimPhysics)


  ///////////////////////////////////////////////////////////////////////////
TSofTwimPhysics::TSofTwimPhysics()
  : m_EventData(new TSofTwimData),
  m_PreTreatedData(new TSofTwimData),
  m_EventPhysics(this),
  m_NumberOfDetectors(0), 
  m_Beta(-1), 
  m_BetaNorm(0.838), 
  m_NumberOfSections(4), 
  m_NumberOfAnodesPerSection(16) {
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSofTwimPhysics::AddDetector(TVector3 ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::BuildSimplePhysicalEvent() {
  if(m_Beta>0)
    BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::BuildPhysicalEvent() {

  //if(m_Beta<0)
  //  return;

  // apply thresholds and calibration
  PreTreat();
  /*if(m_PreTreatedData->GetMultiplicity() != 32){
  //m_Beta = -1;
  return;
  }*/

  vector<double> anode_energy_sec1;
  vector<double> anode_energy_sec2;
  vector<double> anode_energy_sec3;
  vector<double> anode_energy_sec4;
  vector<double> anode_dt_sec1;
  vector<double> anode_dt_sec2;
  vector<double> anode_dt_sec3;
  vector<double> anode_dt_sec4;

  unsigned int mysizeE = m_PreTreatedData->GetMultiplicity();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    int SectionNbr  = m_PreTreatedData->GetSectionNbr(e);
    int AnodeNumber = m_PreTreatedData->GetAnodeNbr(e);
    double Energy   = m_PreTreatedData->GetEnergy(e);
    double DT       = m_PreTreatedData->GetDriftTime(e);

    AnodeSecNbr.push_back(SectionNbr);
    AnodeNbr.push_back(AnodeNumber);
    AnodeEnergy.push_back(Energy);
    AnodeDT.push_back(DT);

    if(SectionNbr==1){
      anode_energy_sec1.push_back(Energy);
      anode_dt_sec1.push_back(DT);
    }
    if(SectionNbr==2){
      anode_energy_sec2.push_back(Energy);
      anode_dt_sec2.push_back(DT);
    }
    if(SectionNbr==3){
      anode_energy_sec3.push_back(Energy);
      anode_dt_sec3.push_back(DT);
    }
    if(SectionNbr==4){
      anode_energy_sec4.push_back(Energy);
      anode_dt_sec4.push_back(DT);
    }
  }

  static CalibrationManager* Cal = CalibrationManager::getInstance();
  double Esec1=0;
  double Esec2=0;
  double Esec3=0;
  double Esec4=0;
  double DTsec1=0;
  double DTsec2=0;
  double DTsec3=0;
  double DTsec4=0;
  for(int i=0; i<anode_energy_sec1.size(); i++){
    Esec1 += anode_energy_sec1[i];
    DTsec1 += anode_dt_sec1[i];
  }
  for(int i=0; i<anode_energy_sec2.size(); i++){
    Esec2 += anode_energy_sec2[i];
    DTsec2 += anode_dt_sec2[i];
  }
  for(int i=0; i<anode_energy_sec3.size(); i++){
    Esec3 += anode_energy_sec3[i];
    DTsec3 += anode_dt_sec3[i];
  }
  for(int i=0; i<anode_energy_sec4.size(); i++){
    Esec4 += anode_energy_sec4[i];
    DTsec4 += anode_dt_sec4[i];
  }

  if(Esec1>0 && anode_energy_sec1.size()>2){
    Esec1 = Esec1 / anode_energy_sec1.size();
    Esec1 = Cal->ApplyCalibration("SofTwim/SEC1_ALIGN",Esec1);

    DTsec1 = DTsec1 / anode_dt_sec1.size();
    EnergySection.push_back(Esec1);
    DriftTime.push_back(DTsec1);
    SectionNbr.push_back(1);
  }
  if(Esec2>0 && anode_energy_sec2.size()>2){
    Esec2 = Esec2 / anode_energy_sec2.size();
    Esec2 = Cal->ApplyCalibration("SofTwim/SEC2_ALIGN",Esec2);

    DTsec2 = DTsec2 / anode_dt_sec2.size();
    EnergySection.push_back(Esec2);
    DriftTime.push_back(DTsec2);
    SectionNbr.push_back(2);
  }
  if(Esec3>0 && anode_energy_sec3.size()>2){
    Esec3 = Esec3 / anode_energy_sec3.size();
    Esec3 = Cal->ApplyCalibration("SofTwim/SEC3_ALIGN",Esec3);

    DTsec3 = DTsec3 / anode_dt_sec3.size();
    EnergySection.push_back(Esec3);
    DriftTime.push_back(DTsec3);
    SectionNbr.push_back(3);
  }
  if(Esec4>0 && anode_energy_sec4.size()>2){
    Esec4 = Esec4 / anode_energy_sec4.size();
    Esec4 = Cal->ApplyCalibration("SofTwim/SEC4_ALIGN",Esec4);

    DTsec4 = DTsec4 / anode_dt_sec4.size();
    EnergySection.push_back(Esec4);
    DriftTime.push_back(DTsec4);
    SectionNbr.push_back(4);
  }

  m_Beta = -1;
}

///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysize = m_EventData->GetMultiplicity();
  for (unsigned int i = 0; i < mysize ; ++i) {
    if(m_EventData->GetPileUp(i) != 1 && m_EventData->GetOverflow(i) != 1){
      double Energy = Cal->ApplyCalibration("SofTwim/SEC"+NPL::itoa(m_EventData->GetSectionNbr(i))+"_ANODE"+NPL::itoa(m_EventData->GetAnodeNbr(i))+"_ENERGY",m_EventData->GetEnergy(i));
      double DT = Cal->ApplyCalibration("SofTwim/SEC"+NPL::itoa(m_EventData->GetSectionNbr(i))+"_ANODE"+NPL::itoa(m_EventData->GetAnodeNbr(i))+"_TIME",m_EventData->GetDriftTime(i));

      m_PreTreatedData->SetSectionNbr(m_EventData->GetSectionNbr(i));
      m_PreTreatedData->SetAnodeNbr(m_EventData->GetAnodeNbr(i));
      m_PreTreatedData->SetEnergy(Energy);
      m_PreTreatedData->SetDriftTime(DT);
      m_PreTreatedData->SetPileUp(m_EventData->GetPileUp(i));
      m_PreTreatedData->SetOverflow(m_EventData->GetOverflow(i));
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSofTwim.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSofTwim.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSofTwim.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSofTwim.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSofTwim";
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

      else if (whatToDo=="SPLINE_SECTION_BETA_PATH") {
        AnalysisConfigFile >> DataBuffer;
        m_SPLINE_SECTION_BETA_PATH = DataBuffer;
        cout << "*** Loading Spline for Beta correction per section ***" << endl;
        LoadSplineBeta();
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
void TSofTwimPhysics::LoadSplineBeta(){
  TString filename = m_SPLINE_SECTION_BETA_PATH;
  TFile* ifile = new TFile(filename,"read");

  if(ifile->IsOpen()){
    cout << "Loading splines..." << endl;
    for(int s=0; s<m_NumberOfSections; s++){
      TString splinename = Form("spline_sec%i",s+1);
      fcorr_beta_sec[s] = (TSpline3*) ifile->FindObjectAny(splinename);
      cout << fcorr_beta_sec[s]->GetName() << endl;
    }
  }
  else
    cout << "File " << filename << " not found!" << endl;
  ifile->Close();
}


///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::Clear() {
  SectionNbr.clear();
  EnergySection.clear();
  DriftTime.clear();
  AnodeNbr.clear();
  AnodeSecNbr.clear();
  AnodeEnergy.clear();
  AnodeDT.clear();
}



///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofTwim");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTwim " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTwim " << i+1 <<  endl;
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
void TSofTwimPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for(int sec = 0; sec < m_NumberOfSections; sec++){
    Cal->AddParameter("SofTwim","SEC"+NPL::itoa(sec+1)+"_ALIGN","SofTwim_SEC"+NPL::itoa(sec+1)+"_ALIGN");

    for(int anode = 0; anode < m_NumberOfAnodesPerSection; anode++){
      Cal->AddParameter("SofTwim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ENERGY","SofTwim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ENERGY");
      Cal->AddParameter("SofTwim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_TIME","SofTwim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_TIME");

    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SofTwim",  true );
  inputChain->SetBranchAddress("SofTwim", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SofTwim", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSofTwimPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("SofTwim", "TSofTwimPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSofTwimPhysics::Construct() {
  return (NPL::VDetector*) new TSofTwimPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_SofTwim{
    public:
      proxy_SofTwim(){
        NPL::DetectorFactory::getInstance()->AddToken("SofTwim","SofTwim");
        NPL::DetectorFactory::getInstance()->AddDetector("SofTwim",TSofTwimPhysics::Construct);
      }
  };

  proxy_SofTwim p_SofTwim;
}

