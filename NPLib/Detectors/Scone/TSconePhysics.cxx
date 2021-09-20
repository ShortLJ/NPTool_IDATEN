/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : March 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Scone Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSconePhysics.h"

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

ClassImp(TSconePhysics)


///////////////////////////////////////////////////////////////////////////
TSconePhysics::TSconePhysics()
   : m_EventData(new TSconeData),
     m_PreTreatedData(new TSconeData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSconePhysics::AddDetector(TVector3 ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

  
///////////////////////////////////////////////////////////////////////////
void TSconePhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  double Energy_det[40];
  double Time_det[40];
  int iTime[40];

  for(int i=0; i<40; i++){
    Energy_det[i] = 0;
    Time_det[i] = 0;
    iTime[i] = 0;
  }
  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  unsigned int mysizeT = m_PreTreatedData->GetMultTime();
  vector<double> tmp_energy;
  vector<int> tmp_det;
  vector<int> tmp_plastic;
  tmp_energy.clear();
  tmp_det.clear();
  tmp_plastic.clear();
  //for (UShort_t e = 0; e < mysizeE ; e++) {
  if (mysizeE == mysizeT) {
    for (UShort_t t = 0; t < mysizeT ; t++) {
      if (m_PreTreatedData->GetE_DetectorNbr(t) == m_PreTreatedData->GetT_DetectorNbr(t)) {
        int det = m_PreTreatedData->GetE_DetectorNbr(t);
        
        tmp_energy.push_back(m_PreTreatedData->Get_Energy(t));
        tmp_det.push_back(m_PreTreatedData->GetE_DetectorNbr(t));
        tmp_plastic.push_back(m_PreTreatedData->GetE_PlasticNbr(t));

        Energy_det[det-1] += m_PreTreatedData->Get_Energy(t);
        
        if(m_PreTreatedData->Get_Time(t)>Time_det[det-1]){
          Time_det[det-1] = m_PreTreatedData->Get_Time(t);
        }
        iTime[det-1]++;
        //DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
        //Energy.push_back(m_PreTreatedData->Get_Energy(e));
        //Time.push_back(m_PreTreatedData->Get_Time(t));
      }
    }
  }

  for(int i=0; i<40; i++){
    //Time_det[i] = Time_det[i]/iTime[i];
    if(Energy_det[i]>m_E_Threshold && Time_det[i]>0){
      DetectorNumber.push_back(i+1);
      Energy.push_back(Energy_det[i]);
      Time.push_back(Time_det[i]);

     /*if(Energy_det[i]>15){
        cout << "**********************************" << endl;
        cout << "event with E<15 MeV !!!!" << endl;
        double Esum=0;
        for(unsigned int k=0; k<tmp_energy.size(); k++){
          if(tmp_det[k]==i+1) {
            cout << i+1 << " / " <<tmp_det[k] << " / " << tmp_plastic[k] << " / " << tmp_energy[k] << endl;
            Esum += tmp_energy[k];
          }
        }
        cout << "final-> " << Energy_det[i] << " / " << Esum << endl;
      }*/
    }
  }

  // Capture Flag
  HasCaptured = m_PreTreatedData->GetCapture(0);

  for(UShort_t i=0; i< m_PreTreatedData->GetGammaMult(); i++){
    GammaEnergy.push_back(m_PreTreatedData->GetGammaEnergy(i));
  }
  for(UShort_t i=0; i< m_PreTreatedData->GetProtonMult(); i++){
    ProtonEnergy.push_back(m_PreTreatedData->GetProtonEnergy(i));
    ProtonTime.push_back(m_PreTreatedData->GetProtonTime(i));
  }

}

///////////////////////////////////////////////////////////////////////////
void TSconePhysics::PreTreat() {
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
      Double_t Energy = Cal->ApplyCalibration("Scone/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), m_EventData->GetE_PlasticNbr(i),Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("Scone/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), m_EventData->GetT_PlasticNbr(i), Time);
  }

  // Capture Flag;
  m_PreTreatedData->SetCapture(m_EventData->GetCapture(0));
  for(UShort_t i=0; i< m_EventData->GetGammaMult(); i++){
    m_PreTreatedData->SetGammaEnergy(m_EventData->GetGammaEnergy(i));
  }
  for(UShort_t i=0; i< m_EventData->GetProtonMult(); i++){
    m_PreTreatedData->SetProtonEnergy(m_EventData->GetProtonEnergy(i));
    m_PreTreatedData->SetProtonTime(m_EventData->GetProtonTime(i));
  }

}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigScone.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigScone.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigScone.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigScone.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigScone";
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
void TSconePhysics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
  HasCaptured = 0;
  GammaEnergy.clear();
  ProtonEnergy.clear();
  ProtonTime.clear();
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Scone");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Ring1","Ring2"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Scone " << i+1 <<  endl;
    
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      m_BuildRing1 = blocks[i]->GetInt("Ring1");
      m_BuildRing2 = blocks[i]->GetInt("Ring2");
      AddDetector(Pos);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSconePhysics::InitSpectra() {
  m_Spectra = new TSconeSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TSconePhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TSconePhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Scone", "D"+ NPL::itoa(i+1)+"_ENERGY","Scone_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("Scone", "D"+ NPL::itoa(i+1)+"_TIME","Scone_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Scone",  true );
  inputChain->SetBranchAddress("Scone", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Scone", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSconePhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Scone", "TSconePhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSconePhysics::Construct() {
  return (NPL::VDetector*) new TSconePhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Scone{
  public:
    proxy_Scone(){
      NPL::DetectorFactory::getInstance()->AddToken("Scone","Scone");
      NPL::DetectorFactory::getInstance()->AddDetector("Scone",TSconePhysics::Construct);
    }
};

proxy_Scone p_Scone;
}

