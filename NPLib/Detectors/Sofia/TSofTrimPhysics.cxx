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
 *  This class hold SofTrim Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSofTrimPhysics.h"

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

ClassImp(TSofTrimPhysics)


  ///////////////////////////////////////////////////////////////////////////
TSofTrimPhysics::TSofTrimPhysics()
  : m_EventData(new TSofTrimData),
  m_PreTreatedData(new TSofTrimData),
  m_EventPhysics(this),
  m_NumberOfDetectors(0), 
  m_Beta(-1),
  m_IsSplineSectionDriftTime(false),
  m_IsSplineSectionAngle(false),
  m_BetaNorm(0.838), 
  m_NumberOfSections(3), 
  m_NumberOfAnodesPaired(3),
  m_NumberOfAnodesPerSection(6) {
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSofTrimPhysics::AddDetector(TVector3 ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::BuildSimplePhysicalEvent() {
  if(m_Beta>0)
    BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::BuildPhysicalEvent() {

  if(m_Beta<0)
    return;

  // apply thresholds and calibration
  PreTreat();
  if(m_PreTreatedData->GetMultiplicity() != 18){
    m_Beta = -1;
    return;
  }

  double Ep1[3], DTp1[3];
  double Ep2[3], DTp2[3];
  double Ep3[3], DTp3[3];
  double Esec[3];
  int mult_p1[3];
  int mult_p2[3];
  int mult_p3[3];
  for(int i=0; i<m_NumberOfSections; i++){
    Ep1[i] = 0;
    Ep2[i] = 0;
    Ep3[i] = 0;
    DTp1[i] = 0;
    DTp2[i] = 0;
    DTp3[i] = 0;
    Esec[i] = 0;
    mult_p1[i] = 0;
    mult_p2[i] = 0;
    mult_p3[i] = 0;
  }

  unsigned int mysizeE = m_PreTreatedData->GetMultiplicity();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    int SectionNbr = m_PreTreatedData->GetSectionNbr(e);
    int AnodeNbr   = m_PreTreatedData->GetAnodeNbr(e);
    double Energy  = m_PreTreatedData->GetEnergy(e);
    double DT      = m_PreTreatedData->GetDriftTime(e);

    //cout << SectionNbr << " " << AnodeNbr << " " << Energy << " " << m_Beta << endl;

    if(AnodeNbr==1 || AnodeNbr==2){ 
      Ep1[SectionNbr-1] += Energy;
      DTp1[SectionNbr-1] += DT;
      mult_p1[SectionNbr-1]++;
    }
    if(AnodeNbr==3 || AnodeNbr==4){
      Ep2[SectionNbr-1] += Energy;
      DTp2[SectionNbr-1] += DT;
      mult_p2[SectionNbr-1]++;
    }
    if(AnodeNbr==5 || AnodeNbr==6){ 
      Ep3[SectionNbr-1] += Energy;
      DTp3[SectionNbr-1] += DT;
      mult_p3[SectionNbr-1]++;
    }
  }

  static CalibrationManager* Cal = CalibrationManager::getInstance();
  for(int i=0; i<m_NumberOfSections; i++){
    if(mult_p1[i] == 2){
      Ep1[i]  = 0.5*Ep1[i];
      DTp1[i] = 0.5*DTp1[i];
    }
    else{
      Ep1[i]  = -1;
      DTp1[i] = -1e5;
    }
    if(mult_p2[i] == 2){
      Ep2[i]  = 0.5*Ep2[i];
      DTp2[i] = 0.5*DTp2[i];
    }
    else{
      Ep2[i]  = -1;
      DTp2[i] = -1e5;
    }
    if(mult_p3[i] == 2){
      Ep3[i]  = 0.5*Ep3[i];
      DTp3[i] = 0.5*DTp3[i];
    }
    else{
      Ep3[i]  = -1;
      DTp3[i] = -1e5;
    }

  }

  double Ddt = DTp2[2] - DTp2[0];
  double p0_1[3], p0_2[3], p0_3[3];
  double p1_1[3], p1_2[3], p1_3[3];

  for(int i=0; i<m_NumberOfSections; i++){
    // Energy Alignement of pairs per section 
    Ep1[i] = Cal->ApplyCalibration("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE1_ALIGN",Ep1[i]);
    Ep2[i] = Cal->ApplyCalibration("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE2_ALIGN",Ep2[i]);
    Ep3[i] = Cal->ApplyCalibration("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE3_ALIGN",Ep3[i]);

    // Beta correction per pair: DE = [0] + [1]*pow(Beta, -5./3);
    p0_1[i] = Cal->GetValue("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE1_BETA",0);
    p0_2[i] = Cal->GetValue("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE2_BETA",0);
    p0_3[i] = Cal->GetValue("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE3_BETA",0);
    p1_1[i] = Cal->GetValue("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE1_BETA",1);
    p1_2[i] = Cal->GetValue("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE2_BETA",1);
    p1_3[i] = Cal->GetValue("SofTrim/SEC"+NPL::itoa(i+1)+"_ANODE3_BETA",1);

    double norm1 = p0_1[i] + p1_1[i]*TMath::Power(m_BetaNorm, -5./3.);
    double norm2 = p0_2[i] + p1_2[i]*TMath::Power(m_BetaNorm, -5./3.);
    double norm3 = p0_3[i] + p1_3[i]*TMath::Power(m_BetaNorm, -5./3.);
    Ep1[i] = norm1 * Ep1[i] / (p0_1[i] + p1_1[i]*TMath::Power(m_Beta, -5./3.));
    Ep2[i] = norm2 * Ep2[i] / (p0_2[i] + p1_2[i]*TMath::Power(m_Beta, -5./3.));
    Ep3[i] = norm3 * Ep3[i] / (p0_3[i] + p1_3[i]*TMath::Power(m_Beta, -5./3.));

    // Angle correction per pair: spline
    /*Ep1[i] = Ep1[i] / fcorr_EvsA[i][0]->Eval(Ddt) * fcorr_EvsA[i][0]->Eval(0);
    Ep2[i] = Ep2[i] / fcorr_EvsA[i][1]->Eval(Ddt) * fcorr_EvsA[i][1]->Eval(0);
    Ep3[i] = Ep3[i] / fcorr_EvsA[i][2]->Eval(Ddt) * fcorr_EvsA[i][2]->Eval(0);

    // DT correction per pair: spline
    Ep1[i] = Ep1[i] / fcorr_EvsDT[i][0]->Eval(DTp1[i]) * fcorr_EvsDT[i][0]->Eval(3000);
    Ep2[i] = Ep2[i] / fcorr_EvsDT[i][1]->Eval(DTp2[i]) * fcorr_EvsDT[i][1]->Eval(3000);
    Ep3[i] = Ep3[i] / fcorr_EvsDT[i][2]->Eval(DTp3[i]) * fcorr_EvsDT[i][2]->Eval(3000);
*/
    // Summing up Anode Energy per section 
    Esec[i] = (Ep1[i] + Ep2[i] + Ep3[i])/3;

    // Angle correction per section: spline   
    if(m_IsSplineSectionAngle)
      Esec[i] = Esec[i] / fcorr_sec_angle[i]->Eval(Ddt) * fcorr_sec_angle[i]->Eval(0);

    // 2nd DT correction per section: spline
    if(m_IsSplineSectionDriftTime)
      Esec[i] = Esec[i] / fcorr_sec_dt[i]->Eval(DTp2[i]) * fcorr_sec_dt[i]->Eval(0);

    // Section ALignement
    Esec[i] = Cal->ApplyCalibration("SofTrim/SEC"+NPL::itoa(i+1)+"_ALIGN",Esec[i]);

    // Filling Output Tree //
    if(DTp2[i]!=0){
      SectionNbr.push_back(i+1);
      EnergyPair1.push_back(Ep1[i]);
      EnergyPair2.push_back(Ep2[i]);
      EnergyPair3.push_back(Ep3[i]);
      DriftTimePair1.push_back(DTp1[i]);
      DriftTimePair2.push_back(DTp2[i]);
      DriftTimePair3.push_back(DTp3[i]);
      EnergySection.push_back(Esec[i]);
      Theta.push_back(DTp2[2]-DTp2[0]);
    }
  }
  m_Beta = -1;
}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysize = m_EventData->GetMultiplicity();
  for (unsigned int i = 0; i < mysize ; ++i) {
    double Energy = Cal->ApplyCalibration("SofTrim/SEC"+NPL::itoa(m_EventData->GetSectionNbr(i))+"_ANODE"+NPL::itoa(m_EventData->GetAnodeNbr(i))+"_ENERGY",m_EventData->GetEnergy(i));
    double DT = Cal->ApplyCalibration("SofTrim/SEC"+NPL::itoa(m_EventData->GetSectionNbr(i))+"_ANODE"+NPL::itoa(m_EventData->GetAnodeNbr(i))+"_TIME",m_EventData->GetDriftTime(i));

    m_PreTreatedData->SetSectionNbr(m_EventData->GetSectionNbr(i));
    m_PreTreatedData->SetAnodeNbr(m_EventData->GetAnodeNbr(i));
    m_PreTreatedData->SetEnergy(Energy);
    m_PreTreatedData->SetDriftTime(DT);
    m_PreTreatedData->SetPileUp(m_EventData->GetPileUp(i));
    m_PreTreatedData->SetOverflow(m_EventData->GetOverflow(i));
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSofTrim.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSofTrim.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSofTrim.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSofTrim.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSofTrim";
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

      else if (whatToDo=="SPLINE_PAIR_ANGLE_PATH") {
        AnalysisConfigFile >> DataBuffer;
        m_SPLINE_PAIR_ANGLE_PATH = DataBuffer;
        cout << "*** Loading Spline for Angle correction per pair ***" << endl;
        LoadSplinePairAngle();
      }

      else if (whatToDo=="SPLINE_PAIR_DT_PATH") {
        AnalysisConfigFile >> DataBuffer;
        m_SPLINE_PAIR_DT_PATH = DataBuffer;
        cout << "*** Loading Spline for Drtft Time correction per pair ***" << endl;
        LoadSplinePairDriftTime();
      }

      else if (whatToDo=="SPLINE_SECTION_DT_PATH") {
        AnalysisConfigFile >> DataBuffer;
        m_SPLINE_SECTION_DT_PATH = DataBuffer;
        cout << "*** Loading Spline for Drift Time correction per section ***" << endl;
        LoadSplineSectionDriftTime();
      }

      else if (whatToDo=="SPLINE_SECTION_ANGLE_PATH") {
        AnalysisConfigFile >> DataBuffer;
        m_SPLINE_SECTION_ANGLE_PATH = DataBuffer;
        cout << "*** Loading Spline for Angle correction per section ***" << endl;
        LoadSplineSectionAngle();
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
double TSofTrimPhysics::GetMaxEnergySection(){
  double Emax=-1;

  if(EnergySection.size() != 3){
    cout << "WARNING! Size of EnergySection different than 3, size= " << EnergySection.size() << endl;
    return Emax;
  }

  Emax = *max_element(EnergySection.begin(), EnergySection.end());

  return Emax;
}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::LoadSplinePairAngle(){
  TString filename = m_SPLINE_PAIR_ANGLE_PATH;
  TFile* ifile = new TFile(filename,"read");

  if(ifile->IsOpen()){
    cout << "Loading splines..." << endl;
    for(int s=0; s<m_NumberOfSections; s++){
      for(int a=0; a<3; a++){
        TString splinename = Form("spline_EvsA_sec%i_anode%i",s+1,a+1);
        fcorr_EvsA[s][a] = (TSpline3*) ifile->FindObjectAny(splinename);
        cout << fcorr_EvsA[s][a]->GetName() << endl;
      }
    }
  }
  else
    cout << "File " << filename << " not found!" << endl;
  ifile->Close();
}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::LoadSplinePairDriftTime(){
  TString filename = m_SPLINE_PAIR_DT_PATH;
  TFile* ifile = new TFile(filename,"read");

  if(ifile->IsOpen()){
    cout << "Loading splines..." << endl;
    for(int s=0; s<m_NumberOfSections; s++){
      for(int a=0; a<3; a++){
        TString splinename = Form("spline_EvsDT_sec%i_anode%i",s+1,a+1);
        fcorr_EvsDT[s][a] = (TSpline3*) ifile->FindObjectAny(splinename);
        cout << fcorr_EvsDT[s][a]->GetName() << endl;
      }
    }
  }
  else
    cout << "File " << filename << " not found!" << endl;
  ifile->Close();
}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::LoadSplineSectionDriftTime(){
  TString filename    = m_SPLINE_SECTION_DT_PATH;
  TFile* ifile    = new TFile(filename,"read");

  if(ifile->IsOpen()){
    cout << "Loading splines..." << endl;
    for(int s=0; s<m_NumberOfSections; s++){
      TString splinename = Form("spline_sec%i",s+1);
      fcorr_sec_dt[s] = (TSpline3*) ifile->FindObjectAny(splinename);
      cout << fcorr_sec_dt[s]->GetName() << endl;
    }
    m_IsSplineSectionDriftTime = true;
  }
  else
    cout << "File " << filename << " not found!" << endl;
  ifile->Close();
}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::LoadSplineSectionAngle(){
  TString filename = m_SPLINE_SECTION_ANGLE_PATH;
  TFile* ifile = new TFile(filename,"read");

  if(ifile->IsOpen()){
    cout << "Loading splines..." << endl;
    for(int s=0; s<m_NumberOfSections; s++){
      TString splinename = Form("spline_sec%i",s+1);
      
      fcorr_sec_angle[s] = (TSpline3*) ifile->FindObjectAny(splinename);
      cout << fcorr_sec_angle[s]->GetName() << endl;
    }
    m_IsSplineSectionAngle = true;
  }
  else
    cout << "File " << filename << " not found!" << endl;
  ifile->Close();
}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::Clear() {
  SectionNbr.clear();
  EnergyPair1.clear();
  EnergyPair2.clear();
  EnergyPair3.clear();
  DriftTimePair1.clear();
  DriftTimePair2.clear();
  DriftTimePair3.clear();
  EnergySection.clear();
  Theta.clear();
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofTrim");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTrim " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTrim " << i+1 <<  endl;
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
void TSofTrimPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for(int sec = 0; sec < m_NumberOfSections; sec++){
    for(int anode = 0; anode < m_NumberOfAnodesPerSection; anode++){
      Cal->AddParameter("SofTrim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ENERGY","SofTrim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ENERGY");
      Cal->AddParameter("SofTrim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_TIME","SofTrim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_TIME");

    }
  }
  for(int sec = 0; sec < m_NumberOfSections; sec++){
    Cal->AddParameter("SofTrim","SEC"+NPL::itoa(sec+1)+"_ALIGN","SofTrim_SEC"+NPL::itoa(sec+1)+"_ALIGN");

    for(int anode = 0; anode < m_NumberOfAnodesPaired; anode++){ 
      Cal->AddParameter("SofTrim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ALIGN","SofTrim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_ALIGN");
      Cal->AddParameter("SofTrim","SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_BETA","SofTrim_SEC"+NPL::itoa(sec+1)+"_ANODE"+NPL::itoa(anode+1)+"_BETA");
    }
  }

}

///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SofTrim",  true );
  inputChain->SetBranchAddress("SofTrim", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SofTrim", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSofTrimPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("SofTrim", "TSofTrimPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSofTrimPhysics::Construct() {
  return (NPL::VDetector*) new TSofTrimPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_SofTrim{
    public:
      proxy_SofTrim(){
        NPL::DetectorFactory::getInstance()->AddToken("SofTrim","SofTrim");
        NPL::DetectorFactory::getInstance()->AddDetector("SofTrim",TSofTrimPhysics::Construct);
      }
  };

  proxy_SofTrim p_SofTrim;
}

