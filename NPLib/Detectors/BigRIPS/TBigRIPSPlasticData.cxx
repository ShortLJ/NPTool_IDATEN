#include "TBigRIPSPlasticData.h"
#include <iostream>

TBigRIPSPlasticData::TBigRIPSPlasticData(){};
TBigRIPSPlasticData::~TBigRIPSPlasticData(){};


////////////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticData::Clear(){
  fPlastic_TL.clear(); 
  fPlastic_TR.clear(); 
  fPlastic_QL.clear(); 
  fPlastic_QR.clear(); 

  fPlastic_TL_ID.clear(); 
  fPlastic_TR_ID.clear(); 
  fPlastic_QL_ID.clear(); 
  fPlastic_QR_ID.clear(); 
}

////////////////////////////////////////////////////////////////////////////////
void TBigRIPSPlasticData::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  //cout << "   - Multiplicity: " << Mult() << endl;

}
ClassImp(TBigRIPSPlasticData); 
