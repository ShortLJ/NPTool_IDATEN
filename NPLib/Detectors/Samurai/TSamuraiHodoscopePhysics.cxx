/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : April 2021                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiHodoscope Treated data                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSamuraiHodoscopePhysics.h"

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

ClassImp(TSamuraiHodoscopePhysics)


  ///////////////////////////////////////////////////////////////////////////
  TSamuraiHodoscopePhysics::TSamuraiHodoscopePhysics(){
    m_EventData = new TSamuraiHodoscopeData;
    m_PreTreatedData = new TSamuraiHodoscopeData;
    m_EventPhysics = this;
    rawQ_low_threshold=0;
    Q_low_threshold=0;
    rawQ_high_threshold=4095;
    Q_high_threshold=300;
    rawT_low_threshold=0;
    T_low_threshold=0;
    rawT_high_threshold=4095;
    T_high_threshold=300;



    //m_Spectra = new TSamuraiHodoscopeSpectra;

  }

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  unsigned int sizeQUp = m_PreTreatedData->GetMultQUp();
  unsigned int sizeTUp = m_PreTreatedData->GetMultTUp();
  unsigned int sizeQDown = m_PreTreatedData->GetMultQDown();
  unsigned int sizeTDown = m_PreTreatedData->GetMultTDown();
  unsigned int IDqu,IDqd,IDtu,IDtd;
  double Qu,Qd,Tu,Td; 
  for(unsigned int qu = 0 ; qu < sizeQUp ; qu++){
    Qu = m_PreTreatedData->GetQUp_Charge(qu);
    Qd = Tu = Td = -1000;
    IDqu = m_PreTreatedData->GetQUp_ID(qu);

    // looking for match q down
    for(unsigned int qd = 0 ; qd < sizeQDown ; qd++){
      IDqd = m_PreTreatedData->GetQDown_ID(qd);
      if(IDqd==IDqu){
        Qd=m_PreTreatedData->GetQDown_Charge(qd); 
        break;
      }
    }

    // valid charge
    if(Qd>0){
      // Looking for associated time 
      for(unsigned int tu = 0 ; tu < sizeTUp ; tu++){
        IDtu =  m_PreTreatedData->GetTUp_ID(tu); 
        if(IDtu==IDqu){
          Tu=m_PreTreatedData->GetTUp_Time(tu);
          break;
        }
      }
      // Valide Time Up
      if(Tu>0){
        for(unsigned int td = 0 ; td < sizeTDown ; td++){
          IDtd =  m_PreTreatedData->GetTDown_ID(td); 
          if(IDtd==IDqu){
            Td=m_PreTreatedData->GetTDown_Time(td);
            break;
          }
        }
      }
      // Both Time are valid
      if(Td>0){
        Charge.push_back(sqrt(Qu*Qd));
        Time.push_back(0.5*(Tu+Td));
        ID.push_back(IDqu);
      } 
    }
  }  
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::PreTreat() {
  m_PreTreatedData->Clear();
  // UP // 
  // Q  // 
  unsigned int sizeQUp = m_EventData->GetMultQUp();
  for(unsigned int i = 0 ; i < sizeQUp ; i++){
    double rawQ = m_EventData->GetQUp_Charge(i)+rand.Uniform(0,1);    
    if(rawQ>rawQ_low_threshold && rawQ<rawQ_high_threshold){

      unsigned int oldID=m_EventData->GetQUp_ID(i);
      double Q = (rawQ-m_qup_offset[oldID])*m_qup_gain[oldID];    
      if(Q>Q_low_threshold)
        m_PreTreatedData->SetChargeUp(m_ID[oldID],Q); 
    }
  }
  // T  // 
  unsigned int sizeTUp = m_EventData->GetMultTUp();
  for(unsigned int i = 0 ; i < sizeTUp ; i++){
    unsigned int oldID=m_EventData->GetTUp_ID(i);
    double rawT = m_EventData->GetTUp_Time(i)+rand.Uniform(0,1); 

    if(rawT>rawT_low_threshold && rawT<rawT_high_threshold){
      double T = rawT*m_tup_gain[oldID]+m_tup_offset[oldID];    
      if(T>T_low_threshold && T<T_high_threshold)
        m_PreTreatedData->SetTimeUp(m_ID[oldID],T); 
    }
  }
  // Down // 
  // Q  // 
  unsigned int sizeQDown = m_EventData->GetMultQDown();
  for(unsigned int i = 0 ; i < sizeQDown ; i++){
    double rawQ = m_EventData->GetQDown_Charge(i)+rand.Uniform(0,1);    
    if(rawQ>rawQ_low_threshold && rawQ<rawQ_high_threshold){

      unsigned int oldID=m_EventData->GetQDown_ID(i);
      double Q = (rawQ-m_qdw_offset[oldID])*m_qdw_gain[oldID];    
      if(Q>Q_low_threshold)
        m_PreTreatedData->SetChargeDown(m_ID[oldID],Q); 
    }
  }
  // T  // 
  unsigned int sizeTDown = m_EventData->GetMultTDown();
  for(unsigned int i = 0 ; i < sizeTDown ; i++){
    unsigned int oldID=m_EventData->GetTDown_ID(i);
    double rawT=m_EventData->GetTDown_Time(i)+rand.Uniform(0,1);
    if(rawT>rawT_low_threshold && rawT<rawT_high_threshold){
      double T = rawT*m_tdw_gain[oldID]+m_tdw_offset[oldID];    
      if(T>T_low_threshold && T<T_high_threshold)
       m_PreTreatedData->SetTimeDown(m_ID[oldID],T); 
    }
  }


}

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::Clear() {
  ID.clear();
  Charge.clear();
  Time.clear();
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SAMURAIHOD");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> mandatory = {"XML",};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(mandatory)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SamuraiHodoscope " << i+1 <<  endl;

      string xmlpath= blocks[i]->GetString("XML");
      NPL::XmlParser xml;
      xml.LoadFile(xmlpath);
      ReadXML(xml);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  // fill the ID to ID map to reorder (probably specific to s034, need to be more modular)
  for(int i = 1 ; i<=24; i++) 
    m_ID[i] = 24-i+1;                                           
  for(int i = 27; i <= 40; i++)                                  
    m_ID[i] = 40-i+27;                                 
  m_ID[25]=26;                     
  m_ID[26]=25;                                          
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::ReadXML(NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("SAMURAIHOD");  
  unsigned int size = b.size();
  for(unsigned int i = 0 ; i < size ; i++){
    unsigned short ID = b[i]->AsInt("ID"); 
    m_qup_gain [ID]   = b[i]->AsDouble("qcal_up");  
    m_qup_offset [ID] = b[i]->AsDouble("qped_up");  
    m_qdw_gain [ID]   = b[i]->AsDouble("qcal_down");  
    m_qdw_offset [ID] = b[i]->AsDouble("qped_down");  
    m_tup_gain [ID]   = b[i]->AsDouble("tup_ch2ns");  
    m_tup_offset [ID] = b[i]->AsDouble("tup_offset");  
    m_tdw_gain [ID]   = b[i]->AsDouble("tdown_ch2ns");  
    m_tdw_offset [ID] = b[i]->AsDouble("tdown_offset");  
  }

}

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::InitSpectra() {
  //m_Spectra = new TSamuraiHodoscopeSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::FillSpectra() {
  //m_Spectra -> FillRawSpectra(m_EventData);
  //m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::CheckSpectra() {
  //m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TSamuraiHodoscopePhysics::GetSpectra() {
  //  if(m_Spectra)
  //return m_Spectra->GetMapHisto();
  // else{
  map< string , TH1*> empty;
  return empty;
  // }
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::WriteSpectra() {
  //  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::AddParameterToCalibrationManager() {
  /*  CalibrationManager* Cal = CalibrationManager::getInstance();
      for (int i = 0; i < m_NumberOfDetectors; ++i) {
      Cal->AddParameter("SamuraiHodoscope", "D"+ NPL::itoa(i+1)+"_ENERGY","SamuraiHodoscope_D"+ NPL::itoa(i+1)+"_ENERGY");
      Cal->AddParameter("SamuraiHodoscope", "D"+ NPL::itoa(i+1)+"_TIME","SamuraiHodoscope_D"+ NPL::itoa(i+1)+"_TIME");
      }*/
}



///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SamuraiHodoscope",  true );
  inputChain->SetBranchAddress("SamuraiHodoscope", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SamuraiHodoscope", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSamuraiHodoscopePhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree("SamuraiHodoscope");
  outputTree->Branch("SamuraiHodoscope", "TSamuraiHodoscopePhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSamuraiHodoscopePhysics::Construct() {
  return (NPL::VDetector*) new TSamuraiHodoscopePhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_SamuraiHodoscope{
    public:
      proxy_SamuraiHodoscope(){
        NPL::DetectorFactory::getInstance()->AddToken("SAMURAIHOD","Samurai");
        NPL::DetectorFactory::getInstance()->AddDetector("SAMURAIHOD",TSamuraiHodoscopePhysics::Construct);
      }
  };

  proxy_SamuraiHodoscope p_SamuraiHodoscope;
}

