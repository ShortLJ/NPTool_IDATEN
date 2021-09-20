/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F, Flavigny   contact address: flavigny@lpccaen.in2p3.fr *
 *                                                                           *
 * Creation Date  : April 2021                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold RIBF IC treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TBigRIPSICPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "NPOptionManager.h"
#include "NPDetectorFactory.h"
#include "NPSystemOfUnits.h"
//   ROOT
using namespace NPUNITS;
///////////////////////////////////////////////////////////////////////////

ClassImp(TBigRIPSICPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TBigRIPSICPhysics::TBigRIPSICPhysics(){
    m_EventData         = new TBigRIPSICData ;
    m_PreTreatedData    = new TBigRIPSICData ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
  }

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::BuildPhysicalEvent(){
  PreTreat();
  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::PreTreat(){

  ClearPreTreatedData();

  //Intermediate map associating a IC ID with its variables stored in an object
  static map<int, BigRIPSICVariables > Fdata ; 
  static map<int, BigRIPSICVariables >::iterator found;
  Fdata.clear() ; 

  int id,layer,ERaw;

  //E
  unsigned int sizeE = m_EventData->GetEMult();
  for(unsigned int i = 0 ; i < sizeE ; i++){
        id = m_EventData->GetEID(i);
        layer = m_EventData->GetELayer(i);
        ERaw = m_EventData->GetE(i);
        Fdata[id].FE.push_back(ERaw); 
        Fdata[id].FE_Layer.push_back(layer); 
  }

  int Nmulti;
  BigRIPSICVariables tmp;

  for(auto it = Fdata.begin();it!=Fdata.end();++it){
    id = it->first;
    ID.push_back(id);
    FP.push_back(IDtoFP[id]);
    tmp.Clear();
    tmp = it->second;
    //tmp.Print(); 
    
    //Calculate E
    double adc, Ecal;
    double AvSum=0;
    double SqSum=1;
    int NLFired=0;

    for (unsigned int j = 0; j < tmp.FE.size(); j++) {
        adc = tmp.FE[j] - pedestal[id][tmp.FE_Layer[j]];
/*
std::cout << ">> FE : " << tmp.FE[j] << std::endl;
std::cout << ">> FE_Layer : " << tmp.FE_Layer[j] << std::endl;
std::cout << ">> ped : " << pedestal[id][tmp.FE_Layer[j]] << std::endl;
std::cout << ">> adc : " << adc << std::endl;
*/
        if(adc>0){
            AvSum = AvSum + adc;
            SqSum = SqSum * adc;
            NLFired++;
            Ecal = ch2mev_0[id] + adc * ch2mev_1[id];
            E.push_back(Ecal);
            E_Layer.push_back(tmp.FE_Layer[j]);
//std::cout << ">> FE_Layer 22222222 : " << tmp.FE_Layer[j] << std::endl;
            E_ID.push_back(id);
        }

    }

    if(NLFired>0){ 
            AvSum = AvSum/NLFired;
            SqSum = std::pow(SqSum,1./NLFired);
            RawAvSum.push_back(AvSum);
            RawSqSum.push_back(SqSum);
            CalAvSum.push_back(ch2mev_0[id] + AvSum * ch2mev_1[id]);
            CalSqSum.push_back(ch2mev_0[id] + SqSum * ch2mev_1[id]);
            NLayerFired.push_back(NLFired);
    }else{
            RawAvSum.push_back(-9999);
            RawSqSum.push_back(-9999);
            CalAvSum.push_back(-9999);
            CalSqSum.push_back(-9999);
    }

  }//end of loop on FData (loop on all plastics)

  //Print(); 
  Fdata.clear();

  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::Clear(){
    ID.clear();
    FP.clear();
    RawAvSum.clear();
    RawSqSum.clear();
    CalAvSum.clear();
    CalSqSum.clear();
    E.clear();
    E_Layer.clear();
    E_ID.clear();
    NLayerFired.clear();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::Print(){
/*
    cout << "XXXXXXXXXXXXXXXXXXXXXXXX IC Physics Event XXXXXXXXXXXXXXXXX" << endl;
    cout << "TL_Mult = " << TL.size();
    for (UShort_t i = 0; i < TL.size(); i++){cout << "\tTL: " << TL[i] << endl;}
    cout << "TR_Mult = " << TR.size();
    for (UShort_t i = 0; i < TR.size(); i++){cout << "\tTR: " << TR[i] << endl;}
    cout << "QL_Mult = " << QL.size();
    for (UShort_t i = 0; i < QL.size(); i++){cout << "\tQL: " << QL[i] << endl;}
    cout << "QR_Mult = " << QR.size();
    for (UShort_t i = 0; i < QR.size(); i++){cout << "\tQR: " << QR[i] << endl;}
*/
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("BigRIPSIC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " IC detector file found " << endl; 

  vector<string> token= {"XML"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    cout << endl << "////  BigRIPSIC file (" << i+1 << ")" << endl;
    string xmlpath = blocks[i]->GetString("XML");
    NPL::XmlParser xml;
    xml.LoadFile(xmlpath);
    AddICs("BigRIPSIC",xml);
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::AddICs(string name, NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName(name);  

  if(name=="BigRIPSIC"){
    unsigned int size = b.size();
    for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = b[i]->AsInt("ID"); 
      IDtoFP[ID] = b[i]->AsInt("FPL"); 
      ch2mev_0[ID] = b[i]->AsDouble("ch2mev_0"); 
      ch2mev_1[ID] = b[i]->AsDouble("ch2mev_1"); 
      for(unsigned int j = 0 ; j < 10 ; j++){
        std::string str_ped = "pedestal" + std::to_string(j) ;  
        pedestal[ID][j] = b[i]->AsInt(str_ped); 
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::InitSpectra(){  
  //m_Spectra = new TBigRIPSICSpectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::FillSpectra(){  
  //  m_Spectra -> FillRawSpectra(m_EventData);
  //  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::CheckSpectra(){  
  //  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TBigRIPSICPhysics::GetSpectra() {
  /*  if(m_Spectra)
      return m_Spectra->GetMapHisto();
      else{
      map< string , TH1*> empty;
      return empty;
      }*/
  map< string , TH1*> empty;
  return empty;

} 

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::WriteSpectra(){
  // m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::AddParameterToCalibrationManager(){
  //CalibrationManager* Cal = CalibrationManager::getInstance();

  // each layer
  //for( int l = 0 ; l < 14 ; ++l){
  //  Cal->AddParameter("SamuraiFDC2", "L"+ NPL::itoa(l),"FDC2_L"+ NPL::itoa(l));
  //}

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "BigRIPSIC" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fIC_*"))
    inputChain->SetBranchStatus( "fIC_*",true);
  inputChain->SetBranchAddress( "BigRIPSIC" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::InitializeRootInputPhysics(){
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSICPhysics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch( "IC" , "TBigRIPSICPhysics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TBigRIPSICPhysics::Construct(){
  return (NPL::VDetector*) new TBigRIPSICPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_bigripsIC{
    public:
      proxy_bigripsIC(){
        NPL::DetectorFactory::getInstance()->AddToken("BigRIPSIC","BigRIPS");
        NPL::DetectorFactory::getInstance()->AddDetector("BigRIPSIC",TBigRIPSICPhysics::Construct);
      }
  };

  proxy_bigripsIC p_bigripsIC;
}

