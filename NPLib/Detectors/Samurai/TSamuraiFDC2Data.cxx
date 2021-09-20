#include "TSamuraiFDC2Data.h"
#include <iostream>

TSamuraiFDC2Data::TSamuraiFDC2Data(){};
TSamuraiFDC2Data::~TSamuraiFDC2Data(){};

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Data::SetData(const int& Det, const int& Layer, const int& Wire, const double& Time, const int& Edge){
  fFDC2_DetectorNbr.push_back(Det);
  fFDC2_LayerNbr.push_back(Layer); 
  fFDC2_WireNbr.push_back(Wire); 
  fFDC2_Time.push_back(Time); 
  fFDC2_Edge.push_back(Edge); 
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Data::Clear(){
  fFDC2_DetectorNbr.clear();
  fFDC2_LayerNbr.clear();
  fFDC2_WireNbr.clear();
  fFDC2_Time.clear();
  fFDC2_Edge.clear();
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Data::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  cout << "   - Multiplicity: " << Mult() << endl;

}
////////////////////////////////////////////////////////////////////////////////
unsigned int TSamuraiFDC2Data::MultLayer(unsigned int det , unsigned int layer, int edge){
  unsigned int mult=0;
  unsigned int size = fFDC2_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fFDC2_DetectorNbr[i]==det)
      if(fFDC2_LayerNbr[i]==layer)
        if(fFDC2_Edge[i]==edge && edge!=-1) // edge type is specified (0 or 1)
          mult++;
        else if(edge==-1)// edge type is not specified
          mult++;
  }
  return mult;

}
////////////////////////////////////////////////////////////////////////////////
std::vector<int> TSamuraiFDC2Data::GetWire(unsigned int det , unsigned int layer){
  std::vector<int> wires;
  unsigned int size = fFDC2_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fFDC2_DetectorNbr[i]==det)
      if(fFDC2_LayerNbr[i]==layer)
        wires.push_back(fFDC2_WireNbr[i]);
  }
  return wires;
}
ClassImp(TSamuraiFDC2Data); 
