/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Valerian Alcindor  contact address: *
 *                                                                           *
 * Creation Date  : October 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold GAGG Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "TGAGGPhysics.h"

//   STL
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdlib.h>
using namespace std;

//   NPL
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "RootInput.h"
#include "RootOutput.h"

//   ROOT
#include "TChain.h"

ClassImp(TGAGGPhysics)

    ///////////////////////////////////////////////////////////////////////////
    TGAGGPhysics::TGAGGPhysics()
    : m_EventData(new TGAGGData), m_PreTreatedData(new TGAGGData),
      m_EventPhysics(this), m_Spectra(0), m_E_RAW_Threshold(0), // adc channels
      m_E_Threshold(0), // MeV
      m_NumberOfDetectors(0) {
  m_Position.clear();
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TGAGGPhysics::AddDetector(TVector3 Pos, string shape) {
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant
  // positions (stripped silicon) or angles (gamma array)
  m_Position.push_back(Pos);
  m_NumberOfDetectors++;
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::AddDetector(double R, double Theta, double Phi,
                               string shape) {
  // Compute the TVector3 corresponding
  TVector3 Pos(R * sin(Theta) * cos(Phi), R * sin(Theta) * sin(Phi),
               R * cos(Theta));
  // Call the cartesian method
  AddDetector(Pos, shape);
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  unsigned int mysizeT = m_PreTreatedData->GetMultTime();
  for (UShort_t e = 0; e < mysizeE; e++) {
    for (UShort_t t = 0; t < mysizeT; t++) {
      if (m_PreTreatedData->GetE_DetectorNbr(e)
          == m_PreTreatedData->GetT_DetectorNbr(t)) {
        DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
        Energy.push_back(m_PreTreatedData->Get_Energy(e));
        Time.push_back(m_PreTreatedData->Get_Time(t));
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
TVector3 TGAGGPhysics::GetPositionOfInteraction(int& i) {
  return m_Position[DetectorNumber[i]];
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Energy
  unsigned int mysize = m_EventData->GetMultEnergy();
  for (UShort_t i = 0; i < mysize; ++i) {
    if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = Cal->ApplyCalibration(
          "GAGG/ENERGY" + NPL::itoa(m_EventData->GetE_DetectorNbr(i)),
          m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time = Cal->ApplyCalibration(
        "GAGG/TIME" + NPL::itoa(m_EventData->GetT_DetectorNbr(i)),
        m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigGAGG.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigGAGG.dat found: Default parameter loaded for Analayis "
         << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigGAGG.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig
      = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigGAGG.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer, DataBuffer, whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigGAGG";
    if (LineBuffer.compare(0, name.length(), name) == 0)
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {
      whatToDo = "";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }

      else if (whatToDo == "E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
      }

      else if (whatToDo == "E_THRESHOLD") {
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
void TGAGGPhysics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("GAGG");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> cart = {"POS", "Shape"};
  vector<string> sphe = {"R", "Theta", "Phi", "Shape"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  GAGG " << i + 1 << endl;

      TVector3 Pos   = blocks[i]->GetTVector3("POS", "mm");
      string   Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos, Shape);
    } else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  GAGG " << i + 1 << endl;
      double R     = blocks[i]->GetDouble("R", "mm");
      double Theta = blocks[i]->GetDouble("Theta", "deg");
      double Phi   = blocks[i]->GetDouble("Phi", "deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R, Theta, Phi, Shape);
    } else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::InitSpectra() {
  m_Spectra = new TGAGGSpectra(m_NumberOfDetectors);
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::FillSpectra() {
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::CheckSpectra() { m_Spectra->CheckSpectra(); }

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::ClearSpectra() {
  // To be done
}

///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TGAGGPhysics::GetSpectra() {
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::WriteSpectra() { m_Spectra->WriteSpectra(); }

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("GAGG", "D" + NPL::itoa(i + 1) + "_ENERGY",
                      "GAGG_D" + NPL::itoa(i + 1) + "_ENERGY");
    Cal->AddParameter("GAGG", "D" + NPL::itoa(i + 1) + "_TIME",
                      "GAGG_D" + NPL::itoa(i + 1) + "_TIME");
  }
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("GAGG", true);
  inputChain->SetBranchAddress("GAGG", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("GAGG", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TGAGGPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("GAGG", "TGAGGPhysics", &m_EventPhysics);
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TGAGGPhysics::Construct() {
  return (NPL::VDetector*)new TGAGGPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_GAGG {
public:
  proxy_GAGG() {
    NPL::DetectorFactory::getInstance()->AddToken("GAGG", "GAGG");
    NPL::DetectorFactory::getInstance()->AddDetector("GAGG",
                                                     TGAGGPhysics::Construct);
  }
};

proxy_GAGG p_GAGG;
}
