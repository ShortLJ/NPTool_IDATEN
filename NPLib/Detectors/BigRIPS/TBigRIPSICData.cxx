#include "TBigRIPSICData.h"
#include <iostream>

TBigRIPSICData::TBigRIPSICData(){};
TBigRIPSICData::~TBigRIPSICData(){};


////////////////////////////////////////////////////////////////////////////////
void TBigRIPSICData::Clear(){
  fIC_E.clear(); 
  fIC_E_ID.clear(); 
  fIC_E_Layer.clear(); 

  fIC_T.clear(); 
  fIC_T_ID.clear(); 
  fIC_T_Layer.clear(); 
}

////////////////////////////////////////////////////////////////////////////////
void TBigRIPSICData::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  //cout << "   - Multiplicity: " << Mult() << endl;

}
ClassImp(TBigRIPSICData); 
