/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Nebula Treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TNebulaPhysics.h"

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

ClassImp(TNebulaPhysics)


  ///////////////////////////////////////////////////////////////////////////
TNebulaPhysics::TNebulaPhysics()
  : m_EventData(new TNebulaData),
  m_EventPhysics(this),
  m_Spectra(0),
  m_Q_RAW_Threshold(0), // adc channels
  m_Q_Threshold(7),     // normal bars in MeV
  m_V_Threshold(1),     // veto bars in MeV
  m_NumberOfBars(0) {
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TNebulaPhysics::ReadXML(NPL::XmlParser xml){ 
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("NEBULA");  

  for(unsigned int i = 0 ; i < b.size() ; i++){
    m_NumberOfBars++;
    unsigned int id = b[i]->AsInt("ID");

    // position
    PositionX[id] = b[i]->AsDouble("PosX"); 
    PositionY[id] = b[i]->AsDouble("PosY"); 
    PositionZ[id] = b[i]->AsDouble("PosZ"); 
    // linear cal
    aQu[id] = b[i]->AsDouble("QUCal");
    bQu[id] = b[i]->AsDouble("QUPed");
    aQd[id] = b[i]->AsDouble("QDCal");
    bQd[id] = b[i]->AsDouble("QDPed");
    aTu[id] = b[i]->AsDouble("TUCal");
    bTu[id] = b[i]->AsDouble("TUOff");
    aTd[id] = b[i]->AsDouble("TDCal");
    bTd[id] = b[i]->AsDouble("TDOff");

    // T average offset
    avgT0[id] = b[i]->AsDouble("TAveOff");

    // slew correction T= tcal +slwT/sqrt(Qcal)
    slwTu[id] = b[i]->AsDouble("TUSlw");
    slwTd[id] = b[i]->AsDouble("TDSlw");

    // DT position cal
    DTa[id] = b[i]->AsDouble("DTCal");//!
    DTb[id] = b[i]->AsDouble("DTOff");//!


  } 
  cout << " -> " << m_NumberOfBars << " bars found" << endl;;
} 

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();
  static double rawQup,calQup,rawQdown,calQdown,rawTup,calTup,rawTdown,calTdown,calQ,calT,Y;
  static unsigned int ID;
  // All vector size 
  static unsigned int QUsize, QDsize, TUsize, TDsize ; 
  QUsize = m_EventData->GetChargeUpMult();
  QDsize = m_EventData->GetChargeDownMult();
  TUsize = m_EventData->GetTimeUpMult();
  TDsize = m_EventData->GetTimeDownMult();
  static double threshold;
  // loop on Qup
  for (unsigned int qup = 0; qup < QUsize ; qup++) {

    rawQup = m_EventData->GetChargeUp(qup);
    rawTup=-1;
    rawQdown=-1;
    rawTdown=-1;
    if (rawQup > m_Q_RAW_Threshold) {
      ID = m_EventData->GetChargeUpID(qup);
      if(ID<121)
        threshold=m_Q_Threshold;
      else
        threshold=m_V_Threshold;

      // look for associated Charge down
      for(unsigned int qdown = 0 ; qdown < QDsize ; qdown++){
        if(m_EventData->GetChargeDownID(qdown)==ID){
          rawQdown=m_EventData->GetChargeDown(qdown); 
          if(rawQdown > m_Q_RAW_Threshold){
            // Look for the associate time 
            for(unsigned int tdown = 0 ; tdown < TDsize; tdown++){
              if(m_EventData->GetTimeDownID(qdown)==ID) {
                rawTdown=m_EventData->GetTimeDown(qdown);
                break;
              }
            }// TDown
          }//if raw threshold down

          break;
        } //if match ID 

      }// Qdwown 

      if(rawTdown>0){ // Tdown is found, means Qdown as well
        // look for Tup  
        for(unsigned int tup = 0 ; tup < TUsize ; tup++){
          if(m_EventData->GetTimeUpID(tup)==ID){
            rawTup = m_EventData->GetTimeUp(tup);
            break;
          }
        }
      }
      // Got everything, do the math
      if(rawTup>0){
        // cal Q Up and Down
        calQup=aQu[ID]*(rawQup-bQu[ID]);
        calQdown=aQd[ID]*(rawQdown-bQd[ID]);
        
        // average value of Up and Down
        calQ=sqrt(calQup*calQdown); 

        // cal T  Up
        calTup=aTu[ID]*rawTup+bTu[ID];
        // slew correction
        calTup -= slwTu[ID]/sqrt(rawQup-bQu[ID]);

        // cal T Down
        calTdown=aTd[ID]*rawTdown+bTd[ID];
        // slew correction
        calTdown -= slwTd[ID]/sqrt(rawQdown-bQd[ID]);

        
        if(calQ>threshold){
          calT= (calTdown+calTup)*0.5+avgT0[ID]+Cal->GetPedestal("NEBULA_T_ID"+NPL::itoa(ID)); 
          Y=(calTdown-calTup)*DTa[ID]+DTb[ID]+Cal->GetPedestal("NEBULA_Y_ID"+NPL::itoa(ID));

          DetectorNumber.push_back(ID);
          Charge.push_back(calQ);
          TOF.push_back(calT);
          PosY.push_back(Y+PositionY[ID]);
          PosX.push_back(PositionX[ID]);
          PosZ.push_back(PositionZ[ID]);

          if(ID<121)
            IsVeto.push_back(0);
          else
            IsVeto.push_back(1);

        }
      }

    }// if raw threshold up
  } // Qup
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::PreTreat() {

}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::ReadAnalysisConfig() {
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::Clear() {
  DetectorNumber.clear();
  Charge.clear();
  TOF.clear();
  PosY.clear();
  PosX.clear();
  PosZ.clear();
  IsVeto.clear();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("NEBULA");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML","Offset","InvertX","InvertY"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      cout << endl << "////  Nebula (" << i+1 << ")" << endl;
      unsigned int det = std::atoi(blocks[i]->GetMainValue().c_str());
      string xmlpath = blocks[i]->GetString("XML");
      NPL::XmlParser xml;
      xml.LoadFile(xmlpath);
      ReadXML(xml);
      TVector3 offset = blocks[i]->GetTVector3("Offset","mm"); 
      bool invertX = blocks[i]->GetInt("InvertX"); 
      bool invertY = blocks[i]->GetInt("InvertY"); 
      m_offset[det] = offset;
      m_invertX[det] = invertX;
      m_invertY[det] = invertY;
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitSpectra() {
  m_Spectra = new TNebulaSpectra(m_NumberOfBars);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TNebulaPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  vector<double> standardO={0};
  for (int i = 0; i < m_NumberOfBars; ++i) {
    Cal->AddParameter("NEBULA_T_ID"+ NPL::itoa(i+1),standardO);
    Cal->AddParameter("NEBULA_Y_ID"+ NPL::itoa(i+1),standardO);
  }
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Nebula",  true );
  inputChain->SetBranchAddress("Nebula", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Nebula", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree("Nebula");
  outputTree->Branch("Nebula", "TNebulaPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TNebulaPhysics::Construct() {
  return (NPL::VDetector*) new TNebulaPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_Nebula{
    public:
      proxy_Nebula(){
        NPL::DetectorFactory::getInstance()->AddToken("NEBULA","Nebula");
        NPL::DetectorFactory::getInstance()->AddDetector("NEBULA",TNebulaPhysics::Construct);
      }
  };

  proxy_Nebula p_Nebula;
}

