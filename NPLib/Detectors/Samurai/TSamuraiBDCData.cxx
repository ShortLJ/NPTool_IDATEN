#include "TSamuraiBDCData.h"
#include <iostream>

TSamuraiBDCData::TSamuraiBDCData(){};
TSamuraiBDCData::~TSamuraiBDCData(){};

////////////////////////////////////////////////////////////////////////////////
void TSamuraiBDCData::SetData(const int& Det, const int& Layer, const int& Wire, const double& Time, const int& Edge){
  fBDC_DetectorNbr.push_back(Det);
  fBDC_LayerNbr.push_back(Layer); 
  fBDC_WireNbr.push_back(Wire); 
  fBDC_Time.push_back(Time); 
  fBDC_Edge.push_back(Edge); 
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiBDCData::Clear(){
  fBDC_DetectorNbr.clear();
  fBDC_LayerNbr.clear();
  fBDC_WireNbr.clear();
  fBDC_Time.clear();
  fBDC_Edge.clear();
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiBDCData::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  cout << "   - Multiplicity: " << Mult() << endl;

}
////////////////////////////////////////////////////////////////////////////////
unsigned int TSamuraiBDCData::MultLayer(unsigned int det , unsigned int layer, int edge){
  unsigned int mult=0;
  unsigned int size = fBDC_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fBDC_DetectorNbr[i]==det)
      if(fBDC_LayerNbr[i]==layer)
        if(fBDC_Edge[i]==edge && edge!=-1) // edge type is specified (0 or 1)
          mult++;
        else if(edge==-1)// edge type is not specified
          mult++;
  }
  return mult;

}
////////////////////////////////////////////////////////////////////////////////
std::vector<int> TSamuraiBDCData::GetWire(unsigned int det , unsigned int layer){
  std::vector<int> wires;
  unsigned int size = fBDC_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fBDC_DetectorNbr[i]==det)
      if(fBDC_LayerNbr[i]==layer)
        wires.push_back(fBDC_WireNbr[i]);
  }
  return wires;
}
ClassImp(TSamuraiBDCData); 
