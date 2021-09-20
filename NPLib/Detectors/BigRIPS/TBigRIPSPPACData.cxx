#include "TBigRIPSPPACData.h"
#include <iostream>

TBigRIPSPPACData::TBigRIPSPPACData(){};
TBigRIPSPPACData::~TBigRIPSPPACData(){};

////////////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACData::Clear(){
  fPPAC_TX1.clear(); 
  fPPAC_TX2.clear(); 
  fPPAC_TY1.clear(); 
  fPPAC_TY2.clear();
  fPPAC_TA.clear(); 
  fPPAC_QX1.clear(); 
  fPPAC_QX2.clear(); 
  fPPAC_QY1.clear(); 
  fPPAC_QY2.clear();
  fPPAC_QA.clear(); 

  fPPAC_TX1_ID.clear(); 
  fPPAC_TX2_ID.clear(); 
  fPPAC_TY1_ID.clear(); 
  fPPAC_TY2_ID.clear();
  fPPAC_TA_ID.clear(); 
  fPPAC_QX1_ID.clear(); 
  fPPAC_QX2_ID.clear(); 
  fPPAC_QY1_ID.clear(); 
  fPPAC_QY2_ID.clear();
  fPPAC_QA_ID.clear(); 
}

////////////////////////////////////////////////////////////////////////////////
void TBigRIPSPPACData::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  //cout << "   - Multiplicity: " << Mult() << endl;

}
ClassImp(TBigRIPSPPACData); 
