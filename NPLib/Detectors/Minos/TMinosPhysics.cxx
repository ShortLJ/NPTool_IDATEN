/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Developed by C. Santamaria / CEA Saclay                  *
 *                                                                           *
 * Creation Date  : 2014/11/24                                               *
 * Last update    : 2019/09 implemeted in NPTool by Cyril Lenain             *
 *                  lenain@lpccaen.in2p3.fr                                  *     
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TMinosPhysics.h"

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

TMinosPhysics* current_phy = 0;

///////////////////////////////////////////////////////////////////////////
TMinosPhysics::TMinosPhysics(){
  m_EventData=new TMinosData;
  m_EventPhysics=this;
  m_E_RAW_Threshold=0; // adc channels
  m_E_Threshold=0;    // MeV
  m_ransac.SetParameters(100,// max 100 iteration
      5,  // a point belong to a track if it is within 5 mm
      40, // max distance allowed between a pair of point to consider a track
      10);// minimum point to form a valid cluster 
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}
///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::BuildPhysicalEvent() {
  PreTreat();
  vector<NPL::LinearCluster3D> clusters = m_ransac.TreatEvent(X_Pad,Y_Pad,Z_Pad);
  unsigned int sizeC = clusters.size();
  for(unsigned int i = 0 ; i < sizeC ; i++){
    // Extract the Direction of the track 
    clusters[i].LinearFit();
    Tracks_P0.push_back(clusters[i].GetP0());
    Tracks_Dir.push_back(clusters[i].GetDir());
  }

  if(sizeC==2){
    static TVector3 Vertex,delta; 
    Delta_Vertex = 
      MinimumDistanceTwoLines(Tracks_P0[0],Tracks_P0[0]+Tracks_Dir[0], 
          Tracks_P0[1],Tracks_P0[1]+Tracks_Dir[1], 
          Vertex, delta);
      X_Vertex= Vertex.X();
      Y_Vertex= Vertex.Y();
      Z_Vertex= Vertex.Z();
      Theta_12=180.*Tracks_Dir[0].Angle(Tracks_Dir[1])/(M_PI);
  }
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::PreTreat() {
  // apply thresholds and calibration
  static unsigned int sizePad,sizeQ ;
  static unsigned short PadNumber;
  static double Q,T;
  static auto cal= CalibrationManager::getInstance();
  static string cal_v, cal_o;

  sizePad = m_EventData->GetPadMult();
  if(sizePad>20){
    for(unsigned int i = 0 ; i < sizePad ; i++){
      vector<unsigned short>* Charge = m_EventData->GetChargePtr(i);
      if(Charge->size()>40 ){
        vector<unsigned short>* Time   = m_EventData->GetTimePtr(i);
        m_utility.Calibrate(Time,Charge,i,T,Q);      

        if(T>0){
          PadNumber = m_EventData->GetPadNumber(i);
          double x_mm = m_X[PadNumber];
          double y_mm = m_Y[PadNumber];
          unsigned int ring = round((sqrt(x_mm*x_mm + y_mm*y_mm)-44.15)/2.1);
          cal_v="Minos/R"+NPL::itoa(ring)+"_VDRIFT";
          cal_o="Minos/R"+NPL::itoa(ring)+"_OFFSET";

          double z_mm = (T*m_TimeBin+cal->GetValue(cal_o,0))*cal->GetValue(cal_v,0);    

          TVector3 Pos=TVector3(x_mm+m_Position.X(),y_mm+m_Position.Y(),z_mm+m_Position.Z());
          Pos.RotateZ(m_ZRotation); 
          // Calibrate the Pad:
          X_Pad.push_back(Pos.X());
          Y_Pad.push_back(Pos.Y());
          Z_Pad.push_back(Pos.Z());    
          Ring_Pad.push_back(ring);
          Q_Pad.push_back(Q);    
          T_Pad.push_back(T);  
        }
      }

      else{
        /* X_Pad.push_back(-1);
           Y_Pad.push_back(-1);
           Z_Pad.push_back((-1);
           Ring_Pad.push_back(-1);
           Z_Pad.push_back(-1);    
           Q_Pad.push_back(-1);    
           T_Pad.push_back(-1);*/  
      }

    }

  }
}
///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;
  // path to file
  string FileName = "./configs/ConfigMinos.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigMinos.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigMinos.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigMinos.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigMinos";
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
void TMinosPhysics::Clear() {
  Ring_Pad.clear();  
  X_Pad.clear();  
  Y_Pad.clear();  
  Z_Pad.clear();
  Q_Pad.clear();
  T_Pad.clear();

  // Tracks information
  Tracks_P0.clear();
  Tracks_Dir.clear();

  // Vertex information
  X_Vertex=-1000;
  Y_Vertex=-1000;
  Z_Vertex=-1000;
  Theta_12=-1000;
  Delta_Vertex=-1000;
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Minos");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML","TimeBin","ShapingTime","Baseline","Sampling","Position","ZRotation"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
        
    if (blocks[i]->HasTokenList(token)) {
      cout << endl << "////  MINOS (" << i+1 << ")" << endl;
      m_ShapingTime = blocks[i]->GetDouble("ShapingTime","ns")/NPUNITS::ns;   
      m_TimeBin = blocks[i]->GetDouble("TimeBin","ns")/NPUNITS::ns;   
      m_Sampling= blocks[i]->GetInt("Sampling");   
      m_Baseline= blocks[i]->GetInt("BaseLine");   
      m_utility.SetParameters(m_TimeBin,m_ShapingTime,m_Baseline,m_Sampling);
      m_Position = blocks[i]->GetTVector3("Position","mm");   
      m_ZRotation= blocks[i]->GetDouble("ZRotation","deg");   
      string xmlpath = blocks[i]->GetString("XML");
      NPL::XmlParser xml;
      xml.LoadFile(xmlpath);
      ReadXML(xml);

    }
    else {
      cout << "ERROR: Missing token for Minos, check your input file"<< endl;
      cout << "Required token: "; 
      for(unsigned int i = 0 ; i < token.size() ; i++)
        cout << token[i] << " " ;
      cout << endl;
      exit(1);
    }

  }
}
///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::ReadXML(NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("MINOS");  
  unsigned int size = b.size();
  unsigned int valid_pad=0;
  for(unsigned int i = 0 ; i < size ; i++){

    unsigned short ID = b[i]->AsInt("ID"); 
    m_X[ID] = b[i]->AsDouble("X");  
    m_Y[ID] = b[i]->AsDouble("Y");  
    valid_pad++;
    if((m_X[ID]==0&&m_Y[ID]==0)||(m_X[ID]==-1&&m_Y[ID]==-1))
      valid_pad--;
    m_Period[ID] = b[i]->AsDouble("TPERIOD");  
    m_Offset[ID] = b[i]->AsDouble("TOFFSET");  
    m_Gain[ID]   = b[i]->AsDouble("GAIN");  
  }

  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << m_X.size() << " pads found with " << valid_pad << " valid pads" << endl; 

}



///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  unsigned int NbrRing = 18;
  vector<double> standardV={0.03475};
  vector<double> standardO={0};
  for(unsigned int i = 0 ; i < NbrRing ; i++){
    Cal->AddParameter("Minos", "R"+NPL::itoa(i+1)+"_VDRIFT","Minos_R"+NPL::itoa(i+1)+"_VDRIFT",standardV);
    Cal->AddParameter("Minos", "R"+NPL::itoa(i+1)+"_OFFSET","Minos_R"+NPL::itoa(i+1)+"_OFFSET",standardO);
  }
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Minos",  true );
  inputChain->SetBranchAddress("Minos", &m_EventData );
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Minos", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree("Minos");
  outputTree->Branch("Minos", "TMinosPhysics", &m_EventPhysics);
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TMinosPhysics::Construct() {
  return (NPL::VDetector*) new TMinosPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_Minos{
    public:
      proxy_Minos(){
        NPL::DetectorFactory::getInstance()->AddToken("Minos","Minos");
        NPL::DetectorFactory::getInstance()->AddDetector("Minos",TMinosPhysics::Construct);
      }
  };

  proxy_Minos p_Minos;
}
