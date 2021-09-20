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
 *  This class hold RIBF Plastic treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TBigRIPSPlasticPhysics.h"

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

ClassImp(TBigRIPSPlasticPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TBigRIPSPlasticPhysics::TBigRIPSPlasticPhysics(){
    m_EventData         = new TBigRIPSPlasticData ;
    m_PreTreatedData    = new TBigRIPSPlasticData ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
  }

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::BuildPhysicalEvent(){
  PreTreat();
  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::PreTreat(){

  ClearPreTreatedData();

  //Intermediate map associating a Plastic ID with its variables stored in an object
  // Condition for filling: Raw TDC value within [LowerLimit,UpperLimit]
  // Only first TDC hit for each variable is stored
  static map<int, BigRIPSPlasticVariables > Fdata ; 
  static map<int, BigRIPSPlasticVariables >::iterator found;
  Fdata.clear() ; 

  //pair of detector ID and variable type (TL,TR,QR,QL)=(0,1,2,3)
  std::pair<unsigned int, double> pair_id_type; 
  int id, TimeRaw, QRaw;


  //TL
  unsigned int sizeTL = m_EventData->GetTLMult();
  for(unsigned int i = 0 ; i < sizeTL ; i++){
        id = m_EventData->GetTLID(i);
        TimeRaw = m_EventData->GetTL(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
            if(Fdata[id].FTL.size()==0){
                Fdata[id].FTL.push_back(TimeRaw); 

            }else {
                Fdata[id].FmultiHit[0]=1; //multiple value for same TDC ch.
            }
        }
  }

  //TR
  unsigned int sizeTR = m_EventData->GetTRMult();
  for(unsigned int i = 0 ; i < sizeTR ; i++){
        id = m_EventData->GetTRID(i);
        TimeRaw = m_EventData->GetTR(i);
        if(TimeRaw>RawLowerLimit[id] && TimeRaw<RawUpperLimit[id]){
            if(Fdata[id].FTR.size()==0){
                Fdata[id].FTR.push_back(TimeRaw); 
            }else {
                Fdata[id].FmultiHit[1]=1; //multiple value for same TDC ch.
            }
        }
  }

  //QL
  unsigned int sizeQL = m_EventData->GetQLMult();
  for(unsigned int i = 0 ; i < sizeQL ; i++){
        id = m_EventData->GetQLID(i);
        QRaw = m_EventData->GetQL(i);
            if(Fdata[id].FQL.size()==0){
                Fdata[id].FQL.push_back(QRaw); 
            }else {
                Fdata[id].FmultiHit[2]=1; //multiple value for same Q ch.
            }
  }

  //QR
  unsigned int sizeQR = m_EventData->GetQRMult();
  for(unsigned int i = 0 ; i < sizeQR ; i++){
        id = m_EventData->GetQRID(i);
        QRaw = m_EventData->GetQR(i);
            if(Fdata[id].FQR.size()==0){
                Fdata[id].FQR.push_back(QRaw); 
            }else {
                Fdata[id].FmultiHit[3]=1; //multiple value for same Q ch.
            }
  }

  int Nmulti;
  BigRIPSPlasticVariables tmp;

  for(auto it = Fdata.begin();it!=Fdata.end();++it){
    id = it->first;
    ID.push_back(id);
    FP.push_back(IDtoFP[id]);
    tmp.Clear();
    tmp = it->second;
    //tmp.Print(); 
    //TL=tmp.FTL * tcal_left[id];
    //TR=tmp.FTR * tcal_right[id];
    //QL=tmp.FQL;
    //QR=tmp.FQR;
    
    //Calculate how many TDC variables got a multi hit
    Nmulti=0;
    for (int i=0; i<4; i++) if(tmp.FmultiHit[i]==1) Nmulti++;
    multiHit.push_back(Nmulti);

    //Calculate TL, TLSlew
    double TL_cal, TLSlew_cal;
    if(tmp.FTL.size()==1){
        TL_cal = tmp.FTL[0] * tcal_left[id];
        TL.push_back(TL_cal);
        if(tmp.HasTLandQL()){            
            TLSlew_cal = tcal_left[id]*
                       (tmp.FTL[0] 
                       + tslew_left_a[id]/sqrt(tmp.FQL[0]) 
                       + tslew_left_b[id]);
            TLSlew.push_back(TLSlew_cal);
       }else{
            TLSlew.push_back(TL_cal);
       }
    }else{ 
       TL.push_back(-99999);
       TLSlew.push_back(-99999);
    }

    //Calculate TR, TRSlew
    double TR_cal, TRSlew_cal;
    if(tmp.FTR.size()==1){
        TR_cal = tmp.FTR[0] * tcal_right[id];
        TR.push_back(TR_cal);
        if(tmp.HasTRandQR()){            
            TRSlew_cal = tcal_right[id]*
                       (tmp.FTR[0] 
                       + tslew_right_a[id]/sqrt(tmp.FQR[0]) 
                       + tslew_right_b[id]);
            TRSlew.push_back(TRSlew_cal);
       }else{
            TRSlew.push_back(TR_cal);
       }
    }else{ 
       TR.push_back(-99999);
       TRSlew.push_back(-99999);
    }

    //Calculate T, TSlew
    if(tmp.FTL.size()==1 && tmp.FTR.size()==1){
       T.push_back((TR_cal+TL_cal)/2.);
       if(tmp.FQL.size()==1 && tmp.FQR.size()==1){
          TSlew.push_back((TRSlew_cal+TLSlew_cal)/2.);
       }else{
          TSlew.push_back(-99999);
       }
    }else{
       T.push_back(-99999);
       TSlew.push_back(-99999);
    }

    //Simply fill QL and QR
    if(tmp.FQL.size()==1){
       QL.push_back(tmp.FQL[0]);
    }else{
       QL.push_back(-99999);
    } 
    if(tmp.FQR.size()==1){
       QR.push_back(tmp.FQR[0]);            
    }else{
       QR.push_back(-99999);
    } 

    if(tmp.HasEverything()==1){
       fired.push_back(true);            
    }else{
       fired.push_back(false);
    } 

  }//end of loop on FData (loop on all plastics)

  //Print(); 
  Fdata.clear();

  return;
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::Clear(){
    ID.clear();
    FP.clear();
    T.clear();
    TL.clear();
    TR.clear();
    TSlew.clear();
    TLSlew.clear();
    TRSlew.clear();
    QL.clear();
    QR.clear();
    multiHit.clear();
    fired.clear();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::Print(){
    cout << "XXXXXXXXXXXXXXXXXXXXXXXX Plastic Physics Event XXXXXXXXXXXXXXXXX" << endl;
    cout << "TL_Mult = " << TL.size();
    for (UShort_t i = 0; i < TL.size(); i++){cout << "\tTL: " << TL[i] << endl;}
    cout << "TR_Mult = " << TR.size();
    for (UShort_t i = 0; i < TR.size(); i++){cout << "\tTR: " << TR[i] << endl;}
    cout << "QL_Mult = " << QL.size();
    for (UShort_t i = 0; i < QL.size(); i++){cout << "\tQL: " << QL[i] << endl;}
    cout << "QR_Mult = " << QR.size();
    for (UShort_t i = 0; i < QR.size(); i++){cout << "\tQR: " << QR[i] << endl;}
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("BigRIPSPlastic");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " Plastic detector file found " << endl; 

  vector<string> token= {"XML"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    cout << endl << "////  BigRIPSPlastic file (" << i+1 << ")" << endl;
    string xmlpath = blocks[i]->GetString("XML");
    NPL::XmlParser xml;
    xml.LoadFile(xmlpath);
    AddPlastics("BigRIPSPlastic",xml);
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::AddPlastics(string name, NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName(name);  

  if(name=="BigRIPSPlastic"){
    unsigned int size = b.size();
    for(unsigned int i = 0 ; i < size ; i++){
      unsigned int ID = b[i]->AsInt("ID"); 
      //string sDir = b[i]->AsString("anodedir");
      RawLowerLimit[ID] = b[i]->AsInt("tdc_underflow"); 
      RawUpperLimit[ID] = b[i]->AsInt("tdc_overflow"); 
      IDtoFP[ID] = b[i]->AsInt("FPL"); 
      tcal_left[ID] = b[i]->AsDouble("tcal_left"); 
      tcal_right[ID] = b[i]->AsDouble("tcal_right"); 
      tslew_left_a[ID] = b[i]->AsDouble("tslew_left_a"); 
      tslew_left_b[ID] = b[i]->AsDouble("tslew_left_b"); 
      tslew_right_a[ID] = b[i]->AsDouble("tslew_right_a"); 
      tslew_right_b[ID] = b[i]->AsDouble("tslew_right_b"); 
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::InitSpectra(){  
  //m_Spectra = new TBigRIPSPlasticSpectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::FillSpectra(){  
  //  m_Spectra -> FillRawSpectra(m_EventData);
  //  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::CheckSpectra(){  
  //  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TBigRIPSPlasticPhysics::GetSpectra() {
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
void TBigRIPSPlasticPhysics::WriteSpectra(){
  // m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::AddParameterToCalibrationManager(){
  //CalibrationManager* Cal = CalibrationManager::getInstance();

  // each layer
  //for( int l = 0 ; l < 14 ; ++l){
  //  Cal->AddParameter("SamuraiFDC2", "L"+ NPL::itoa(l),"FDC2_L"+ NPL::itoa(l));
  //}

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "BigRIPSPlastic" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fPlastic_*"))
    inputChain->SetBranchStatus( "fPlastic_*",true);
  inputChain->SetBranchAddress( "BigRIPSPlastic" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::InitializeRootInputPhysics(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "BigRIPSPlastic" , true );
  inputChain->SetBranchAddress( "BigRIPSPlastic" , &m_EventPhysics );
}

///////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticPhysics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch( "BigRIPSPlastic" , "TBigRIPSPlasticPhysics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TBigRIPSPlasticPhysics::Construct(){
  return (NPL::VDetector*) new TBigRIPSPlasticPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_bigripsPlastic{
    public:
      proxy_bigripsPlastic(){
        NPL::DetectorFactory::getInstance()->AddToken("BigRIPSPlastic","BigRIPS");
        NPL::DetectorFactory::getInstance()->AddDetector("BigRIPSPlastic",TBigRIPSPlasticPhysics::Construct);
      }
  };

  proxy_bigripsPlastic p_bigripsPlastic;
}

