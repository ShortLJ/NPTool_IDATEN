#include "TSamuraiFDC0Data.h"
#include <iostream>

TSamuraiFDC0Data::TSamuraiFDC0Data(){};
TSamuraiFDC0Data::~TSamuraiFDC0Data(){};

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC0Data::SetData(const int& Det, const int& Layer, const int& Wire, const double& Time, const int& Edge){
  fFDC0_DetectorNbr.push_back(Det);
  fFDC0_LayerNbr.push_back(Layer); 
  fFDC0_WireNbr.push_back(Wire); 
  fFDC0_Time.push_back(Time); 
  fFDC0_Edge.push_back(Edge); 
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC0Data::Clear(){
  fFDC0_DetectorNbr.clear();
  fFDC0_LayerNbr.clear();
  fFDC0_WireNbr.clear();
  fFDC0_Time.clear();
  fFDC0_Edge.clear();
}

////////////////////////////////////////////////////////////////////////////////
void TSamuraiFDC0Data::Print(){
  using namespace std;

  cout << " -- Event:" << endl;
  cout << "   - Multiplicity: " << Mult() << endl;

}
////////////////////////////////////////////////////////////////////////////////
unsigned int TSamuraiFDC0Data::MultLayer(unsigned int det , unsigned int layer, int edge){
  unsigned int mult=0;
  unsigned int size = fFDC0_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fFDC0_DetectorNbr[i]==det)
      if(fFDC0_LayerNbr[i]==layer)
        if(fFDC0_Edge[i]==edge && edge!=-1) // edge type is specified (0 or 1)
          mult++;
        else if(edge==-1)// edge type is not specified
          mult++;
  }
  return mult;

}
////////////////////////////////////////////////////////////////////////////////
std::vector<int> TSamuraiFDC0Data::GetWire(unsigned int det , unsigned int layer){
  std::vector<int> wires;
  unsigned int size = fFDC0_DetectorNbr.size();
  for(unsigned int i = 0 ; i< size ; i++ ){
    if(fFDC0_DetectorNbr[i]==det)
      if(fFDC0_LayerNbr[i]==layer)
        wires.push_back(fFDC0_WireNbr[i]);
  }
  return wires;
}
ClassImp(TSamuraiFDC0Data); 
