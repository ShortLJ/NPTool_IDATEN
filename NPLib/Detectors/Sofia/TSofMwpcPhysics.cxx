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
 *  This class hold SofMwpc Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSofMwpcPhysics.h"

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

ClassImp(TSofMwpcPhysics)


  ///////////////////////////////////////////////////////////////////////////
TSofMwpcPhysics::TSofMwpcPhysics()
  : m_EventData(new TSofMwpcData),
  m_PreTreatedData(new TSofMwpcData),
  m_EventPhysics(this),
  m_NumberOfDetectors(0){
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSofMwpcPhysics::AddDetector(TVector3 Pos){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  DetPosX.push_back(Pos.X());
  DetPosY.push_back(Pos.Y());
  DetPosZ.push_back(Pos.Z());

  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::BuildPhysicalEvent() {

  PreTreat();

  unsigned int mysizeE = m_PreTreatedData->GetMultiplicity();

  vector<double> ChargeArray;
  ChargeArray.resize(288,0);

  int StripMaxX1[4]={0,0,0,0};
  int StripMaxX2[4]={0,0,0,0};
  int StripMaxY[4] ={0,0,0,0};
  double ChargeMaxX1[4]={0,0,0,0};
  double ChargeMaxX2[4]={0,0,0,0};
  double ChargeMaxY[4] ={0,0,0,0};


  for(int i=0; i<m_NumberOfDetectors; i++){
    Buffer_X1_Q.push_back(ChargeArray);
    Buffer_X2_Q.push_back(ChargeArray);
    Buffer_Y_Q.push_back(ChargeArray);
  }

  for (UShort_t e = 0; e < mysizeE ; e++) {
    int det = m_PreTreatedData->GetDetectorNbr(e);
    int plane = m_PreTreatedData->GetPlane(e);
    int strip = m_PreTreatedData->GetPad(e);// starts at 0
    double charge = m_PreTreatedData->GetCharge(e);

    // Xup for MWPC2 and 3 and X for MWPC1 and 4
    if(plane==1){
      Buffer_X1_Q[det-1][strip] = charge;

      if(charge>ChargeMaxX1[det-1]){
        ChargeMaxX1[det-1] = charge;
        StripMaxX1[det-1] = strip;
      }
    }
    // Xdown for MWPC2 and 3 and nothing for MWPC1 and 4
    else if(plane==2){
      Buffer_X2_Q[det-1][strip] = charge;

      if(charge>ChargeMaxX2[det-1]){
        ChargeMaxX2[det-1] = charge;
        StripMaxX2[det-1] = strip;
      }
    }
    // Y for all MWPCx
    if(plane==3){
      Buffer_Y_Q[det-1][strip] = charge;

      if(charge>ChargeMaxY[det-1]){
        ChargeMaxY[det-1] = charge;
        StripMaxY[det-1] = strip;
      }
    }
  }

  double X1 = -100;
  double X2 = -100;
  double Y  = -100;
  for(int i=0; i<m_NumberOfDetectors; i++){
    double qleft  = Buffer_X1_Q[i][StripMaxX1[i]-1];
    double qright = Buffer_X1_Q[i][StripMaxX1[i]+1];
    if(qleft>0 && qright>0){
      X1 = GetPositionX(ChargeMaxX1[i], StripMaxX1[i], qleft, qright);
      X1 = X1 - DetPosX[i];
    }
    qleft  = Buffer_X2_Q[i][StripMaxX2[i]-1];
    qright = Buffer_X2_Q[i][StripMaxX2[i]+1];
    if(qleft>0 && qright>0){
      X2 = GetPositionX(ChargeMaxX2[i], StripMaxX2[i], qleft, qright);
      X2 = X2 - DetPosX[i];
    }
    double qdown  = Buffer_Y_Q[i][StripMaxY[i]-1];
    double qup = Buffer_Y_Q[i][StripMaxY[i]+1];
    if(qdown>0 && qup>0){
      Y = GetPositionY(ChargeMaxY[i], StripMaxY[i], qdown, qup);
      Y = Y - DetPosY[i];
    }

    DetectorNbr.push_back(i+1);
    PositionX1.push_back(X1);
    PositionX2.push_back(X2);
    PositionY.push_back(Y);
  }
}

///////////////////////////////////////////////////////////////////////////
double TSofMwpcPhysics::GetPositionX(double qmax, int padmax, double qleft, double qright){
  double fwx = 3.125;
  double fSize = 200.0;

  double a3 = TMath::Pi() * fwx / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qleft) + TMath::Sqrt(qmax / qright))));

  Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qleft) - TMath::Sqrt(qmax / qright)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwx / a3)));

  return (-1. * padmax * fwx + (fSize / 2) - (fwx / 2) - a2); // Left is positive and right negative
}

///////////////////////////////////////////////////////////////////////////
double TSofMwpcPhysics::GetPositionY(double qmax, int padmax, double qdown, double qup){
  double fwy = 3.125;
  double fSize = 200.0;
  double a3 = TMath::Pi() * fwy / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qdown) + TMath::Sqrt(qmax / qup))));

  Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qdown) - TMath::Sqrt(qmax / qup)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwy / a3)));

  return (padmax * fwy - (fSize / 2) + (fwy / 2) + a2);
}

///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysize = m_EventData->GetMultiplicity();
  for (unsigned int i = 0; i < mysize ; ++i) {
    m_PreTreatedData->SetDetectorNbr(m_EventData->GetDetectorNbr(i));
    m_PreTreatedData->SetPlane(m_EventData->GetPlane(i));
    m_PreTreatedData->SetPad(m_EventData->GetPad(i));
    m_PreTreatedData->SetCharge(m_EventData->GetCharge(i));
  }
}



///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSofMwpc.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSofMwpc.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSofMwpc.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSofMwpc.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSofMwpc";
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

      else {
        ReadingStatus = false;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::Clear() {
  DetectorNbr.clear();
  PositionX1.clear();
  PositionX2.clear();
  PositionY.clear();

  Buffer_X1_Q.clear();
  Buffer_X2_Q.clear();
  Buffer_Y_Q.clear();
}



///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofMwpc");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofMwpc " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofMwpc " << i+1 <<  endl;
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
void TSofMwpcPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for(int sec = 0; sec < m_NumberOfDetectors; sec++){
    Cal->AddParameter("SofMwpc","SEC"+NPL::itoa(sec+1)+"_ALIGN","SofMwpc_SEC"+NPL::itoa(sec+1)+"_ALIGN");
  }
}

///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SofMwpc",  true );
  inputChain->SetBranchAddress("SofMwpc", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SofMwpc", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSofMwpcPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("SofMwpc", "TSofMwpcPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSofMwpcPhysics::Construct() {
  return (NPL::VDetector*) new TSofMwpcPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_SofMwpc{
    public:
      proxy_SofMwpc(){
        NPL::DetectorFactory::getInstance()->AddToken("SofMwpc","SofMwpc");
        NPL::DetectorFactory::getInstance()->AddDetector("SofMwpc",TSofMwpcPhysics::Construct);
      }
  };

  proxy_SofMwpc p_SofMwpc;
}

