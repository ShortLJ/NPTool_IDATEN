#ifndef TSamuraiFDC0Data_H
#define TSamuraiFDC0Data_H
#include "TObject.h"
#include <vector>

class TSamuraiFDC0Data: public TObject{
  public:
    TSamuraiFDC0Data();
    ~TSamuraiFDC0Data();

  private:
    std::vector<int> fFDC0_DetectorNbr;
    std::vector<int> fFDC0_LayerNbr;
    std::vector<int> fFDC0_WireNbr;
    std::vector<double> fFDC0_Time;
    std::vector<int> fFDC0_Edge;
  
  public:
    void Clear();
    void Print();
    void Clear(const Option_t*) {};
    void Dump() const{};
  
  public:
    void SetData(const int& Det, const int& Layer, const int& Wire, const double& Time, const int& Edge);
    unsigned int Mult(){return fFDC0_DetectorNbr.size();};
    unsigned int MultLayer(unsigned int det , unsigned int layer, int edge=-1);
    std::vector<int> GetWire(unsigned int det , unsigned int layer);
    int const GetDetectorNbr(const unsigned int& i){return fFDC0_DetectorNbr[i];};
    int const GetLayerNbr(const unsigned int& i){return fFDC0_LayerNbr[i];};
    int const GetWireNbr(const unsigned int& i){return fFDC0_WireNbr[i];};
    double const GetTime(const unsigned int& i){return fFDC0_Time[i];};
    int const GetEdge(const unsigned int& i){return fFDC0_Edge[i];};

    ClassDef(TSamuraiFDC0Data,1); 
};

#endif
