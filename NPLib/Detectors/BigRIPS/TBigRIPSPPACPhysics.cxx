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
 *  This class hold RIBF PPAC treated data                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TBigRIPSPPACPhysics.h"

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

ClassImp(TBigRIPSPPACPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TBigRIPSPPACPhysics::TBigRIPSPPACPhysics(){
    m_EventData         = new TBigRIPSPPACData ;
    m_PreTreatedData    = new TBigRIPSPPACData ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
  }

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::BuildPhysicalEvent(){
  PreTreat();
  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::PreTreat(){

  ClearPreTreatedData();

  //Intermediate map associating a PPAC ID with its variables stored in an object
  // Condition for filling: Raw TDC value within [LowerLimit,UpperLimit]
  // Only first TDC hit for each variable is stored
  static map<int, BigRIPSPPACVariables > Fdata ; 
  static map<int, BigRIPSPPACVariables >::iterator found;
  Fdata.clear() ; 

  //pair of detector ID and variable type (TX1,TX2,TY1,TY2,TA)=(0,1,2,3,4)
  std::pair<unsigned int, double> pair_id_type; 
  int id, TimeRaw;


  //TX1
  unsigned int sizeTX1 = m_EventData->GetTX1Mult();
  for(unsigned int i = 0 ; i < sizeTX1 ; i++){
        id = m_EventData->GetTX1ID(i);
        TimeRaw = m_EventData->GetTX1(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
            if(Fdata[id].FTX1.size()==0){
                Fdata[id].FTX1.push_back(TimeRaw*ch2ns_TX1[id]); 
            }else {
                Fdata[id].FmultiHit[0]=1; //multiple value for same TDC ch.
            }
        }
  }

  //TX2
  unsigned int sizeTX2 = m_EventData->GetTX2Mult();
  for(unsigned int i = 0 ; i < sizeTX2 ; i++){
        id = m_EventData->GetTX2ID(i);
        TimeRaw = m_EventData->GetTX2(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
            if(Fdata[id].FTX2.size()==0){
                Fdata[id].FTX2.push_back(TimeRaw*ch2ns_TX2[id]); 
            }else {
                Fdata[id].FmultiHit[1]=1; //multiple value for same TDC ch.
            }
        }
  }

  //TY1
  unsigned int sizeTY1 = m_EventData->GetTY1Mult();
  for(unsigned int i = 0 ; i < sizeTY1 ; i++){
        id = m_EventData->GetTY1ID(i);
        TimeRaw = m_EventData->GetTY1(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
            if(Fdata[id].FTY1.size()==0){
                Fdata[id].FTY1.push_back(TimeRaw*ch2ns_TY1[id]); 
            }else {
                Fdata[id].FmultiHit[2]=1; //multiple value for same TDC ch.
            }
        }
  }

  //TY2
  unsigned int sizeTY2 = m_EventData->GetTY2Mult();
  for(unsigned int i = 0 ; i < sizeTY2 ; i++){
        id = m_EventData->GetTY2ID(i);
        TimeRaw = m_EventData->GetTY2(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
            if(Fdata[id].FTY2.size()==0){
                Fdata[id].FTY2.push_back(TimeRaw*ch2ns_TY2[id]); 
            }else {
                Fdata[id].FmultiHit[3]=1; //multiple value for same TDC ch.
            }
        }
  }

  //TA
  unsigned int sizeTA = m_EventData->GetTAMult();
  for(unsigned int i = 0 ; i < sizeTA ; i++){
        id = m_EventData->GetTAID(i);
        TimeRaw = m_EventData->GetTA(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
            if(Fdata[id].FTA.size()==0){
                Fdata[id].FTA.push_back(TimeRaw*ch2ns_TA[id]); 
            }else {
                Fdata[id].FmultiHit[4]=1; //multiple value for same TDC ch.
            }

        }
  }

  int Nmulti;
  BigRIPSPPACVariables tmp;

  for(auto it = Fdata.begin();it!=Fdata.end();++it){
    id = it->first;
    ID.push_back(id);
    FP.push_back(IDtoFP[id]);
    tmp.Clear();
    tmp = it->second;
    //tmp.Print(); 
    TX1=tmp.FTX1;
    TX2=tmp.FTX2;
    TY1=tmp.FTY1;
    TY2=tmp.FTY2;
    TA=tmp.FTA;
    
    //Calculate how many TDC variables got a multi hit
    Nmulti=0;
    for (int i=0; i<5; i++) if(tmp.FmultiHit[i]==1) Nmulti++;
    multiHit.push_back(Nmulti);

    //Calculate TSum, TDiff and X/Y
     double TSumX_tmp, TSumY_tmp;
     double TDiffX_tmp, TDiffY_tmp;
     double X_tmp, Y_tmp;
    
    if(tmp.HasTXs() && tmp.HasTA()){
        TSumX_tmp = TX1[0]+TX2[0]-2*TA[0];
        TSumX.push_back(TSumX_tmp);
    }else {
        TSumX_tmp = -99999;
        TSumX.push_back(TSumX_tmp);
    }

    if(tmp.HasTXs() && ( ignore_txsum_cut[id] || 
            (TSumX_tmp >= txsum_min[id] && TSumX_tmp <= txsum_max[id]) ) ){
        TDiffX_tmp = (TX1[0] - TX2[0] - xns_off[id]);
        TDiffX.push_back(TDiffX_tmp);
        X_tmp = -1*((TDiffX_tmp * x_ns2mm[id])/2. - x_offset[id] - xpos_offset[id]);
        X.push_back(X_tmp);
    }else{
        TDiffX.push_back(-99999);
        X.push_back(-99999);
    }

    if(tmp.HasTYs() && tmp.HasTA()) {
        TSumY_tmp = TY1[0]+TY2[0]-2*TA[0];
        TSumY.push_back(TSumY_tmp);
    } else {
        TSumY_tmp = -99999;
        TSumY.push_back(TSumY_tmp);
    }

    if(tmp.HasTYs() && (ignore_tysum_cut[id] || 
         (TSumY_tmp >= tysum_min[id] && TSumY_tmp <= tysum_max[id])) ){
        TDiffY_tmp = (TY1[0] - TY2[0] - yns_off[id]);
        TDiffY.push_back(TDiffY_tmp);
        Y_tmp = ((TDiffY_tmp * y_ns2mm[id])/2. - y_offset[id] - ypos_offset[id]);
        Y.push_back(Y_tmp);
/*
        if(id==19){
            std::cout << "_____________"<< std::endl;
            std::cout << "TY1:\t\t" << TY1[0] << std::endl;
            std::cout << "TY2:\t\t" << TY2[0] << std::endl;
            std::cout << "yns_off:\t" <<  yns_off[id]<< std::endl;
            std::cout << "TDiffY:\t\t" << TDiffY_tmp << std::endl;
            std::cout << "y_ns2mm:\t" <<  y_ns2mm[id]<< std::endl;
            std::cout << "y_off:\t\t" <<  y_offset[id]<< std::endl;
            std::cout << "ypos_off:\t" <<  ypos_offset[id]<< std::endl;
            std::cout << "Y:\t\t" <<  Y_tmp<< std::endl;
        }
*/
    }else {
        TDiffY.push_back(-99999);
        Y.push_back(-99999);
    }
    

  }
  //Print(); 
  Fdata.clear();

  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::Clear(){
    TX1.clear();
    TX2.clear();
    TY1.clear();
    TY2.clear();
    TA.clear();
    TSumX.clear();
    TDiffX.clear();
    X.clear();
    TSumY.clear();
    TDiffY.clear();
    Y.clear();
    ID.clear();
    FP.clear();
    multiHit.clear();
    //Data.clear();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::Print(){
    cout << "XXXXXXXXXXXXXXXXXXXXXXXX PPAC Physics Event XXXXXXXXXXXXXXXXX" << endl;
    cout << "TX1_Mult = " << TX1.size();
    for (UShort_t i = 0; i < TX1.size(); i++){cout << "\tTX1: " << TX1[i] << endl;}
    cout << "TX2_Mult = " << TX2.size();
    for (UShort_t i = 0; i < TX2.size(); i++){cout << "\tTX2: " << TX2[i] << endl;}
    cout << "TY1_Mult = " << TY1.size();
    for (UShort_t i = 0; i < TY1.size(); i++){cout << "\tTY1: " << TY1[i] << endl;}
    cout << "TY2_Mult = " << TY2.size();
    for (UShort_t i = 0; i < TY2.size(); i++){cout << "\tTY2: " << TY2[i] << endl;}
    cout << "TA_Mult = " << TA.size();
    for (UShort_t i = 0; i < TA.size(); i++){cout << "\tTA: " << TA[i] << endl;}
    cout << "TDiffX_Mult = " << TDiffX.size();
    for (UShort_t i = 0; i < TDiffX.size(); i++){cout << "\tTDiffX: " << TDiffX[i] << endl;}
    cout << "TDiffY_Mult = " << TDiffY.size();
    for (UShort_t i = 0; i < TDiffY.size(); i++){cout << "\tTDiffY: " << TDiffY[i] << endl;}
    cout << "TSumX_Mult = " << TSumX.size();
    for (UShort_t i = 0; i < TSumX.size(); i++){cout << "\tTSumX: " << TSumX[i] << endl;}
    cout << "TSumY_Mult = " << TSumY.size();
    for (UShort_t i = 0; i < TSumY.size(); i++){cout << "\tTSumY: " << TSumY[i] << endl;}
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("BigRIPSPPAC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " PPAC detector file found " << endl; 

  vector<string> token= {"XML"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    cout << endl << "////  BigRIPSPPAC file (" << i+1 << ")" << endl;
    string xmlpath = blocks[i]->GetString("XML");
    NPL::XmlParser xml;
    xml.LoadFile(xmlpath);
    AddPPACs("BigRIPSPPAC",xml);
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::AddPPACs(string name, NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName(name);  

  if(name=="BigRIPSPPAC"){
    unsigned int size = b.size();
    for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = b[i]->AsInt("ID"); 
      //string sDir = b[i]->AsString("anodedir");
      RawLowerLimit[ID] = b[i]->AsInt("tdc_underflow"); 
      RawUpperLimit[ID] = b[i]->AsInt("tdc_overflow"); 
      IDtoFP[ID] = b[i]->AsInt("FPL"); 
      ch2ns_TX1[ID] = b[i]->AsDouble("x1_ch2ns"); 
      ch2ns_TX2[ID] = b[i]->AsDouble("x2_ch2ns"); 
      ch2ns_TY1[ID] = b[i]->AsDouble("y1_ch2ns"); 
      ch2ns_TY2[ID] = b[i]->AsDouble("y2_ch2ns"); 
      ch2ns_TA[ID] = b[i]->AsDouble("a_ch2ns"); 
      xns_off[ID] = b[i]->AsDouble("xns_off"); 
      yns_off[ID] = b[i]->AsDouble("yns_off"); 
      x_ns2mm[ID] = b[i]->AsDouble("xfactor"); //prop. speed along delay line 
      y_ns2mm[ID] = b[i]->AsDouble("yfactor"); //prop. speed along delay line 
      x_offset[ID] = b[i]->AsDouble("xoffset"); 
      y_offset[ID] = b[i]->AsDouble("yoffset"); 
      xpos_offset[ID] = b[i]->AsDouble("xpos_off"); 
      ypos_offset[ID] = b[i]->AsDouble("ypos_off"); 
      txsum_min[ID] = b[i]->AsDouble("txsum_min"); 
      txsum_max[ID] = b[i]->AsDouble("txsum_max"); 
      tysum_min[ID] = b[i]->AsDouble("tysum_min"); 
      tysum_max[ID] = b[i]->AsDouble("tysum_max"); 
      ignore_txsum_cut[ID] = false; 
      if(txsum_min[ID] >= txsum_max[ID]) ignore_txsum_cut[ID] = true; 
      ignore_tysum_cut[ID] = false; 
      if(tysum_min[ID] >= tysum_max[ID]) ignore_tysum_cut[ID] = true; 
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitSpectra(){  
  //m_Spectra = new TBigRIPSPPACSpectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::FillSpectra(){  
  //  m_Spectra -> FillRawSpectra(m_EventData);
  //  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::CheckSpectra(){  
  //  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TBigRIPSPPACPhysics::GetSpectra() {
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
void TBigRIPSPPACPhysics::WriteSpectra(){
  // m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::AddParameterToCalibrationManager(){
  //CalibrationManager* Cal = CalibrationManager::getInstance();

  // each layer
  //for( int l = 0 ; l < 14 ; ++l){
  //  Cal->AddParameter("SamuraiFDC2", "L"+ NPL::itoa(l),"FDC2_L"+ NPL::itoa(l));
  //}

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "BigRIPSPPAC" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fPPAC_*"))
    inputChain->SetBranchStatus( "fPPAC_*",true);
  inputChain->SetBranchAddress( "BigRIPSPPAC" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitializeRootInputPhysics(){
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACPhysics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch( "PPAC" , "TBigRIPSPPACPhysics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TBigRIPSPPACPhysics::Construct(){
  return (NPL::VDetector*) new TBigRIPSPPACPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_bigripsPPAC{
    public:
      proxy_bigripsPPAC(){
        NPL::DetectorFactory::getInstance()->AddToken("BigRIPSPPAC","BigRIPS");
        NPL::DetectorFactory::getInstance()->AddDetector("BigRIPSPPAC",TBigRIPSPPACPhysics::Construct);
      }
  };

  proxy_bigripsPPAC p_bigripsPPAC;
}

