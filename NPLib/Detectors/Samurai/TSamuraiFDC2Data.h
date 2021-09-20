#ifndef TSamuraiFDC2Data_H
#define TSamuraiFDC2Data_H
#include "TObject.h"
#include <vector>

class TSamuraiFDC2Data: public TObject{
  public:
    TSamuraiFDC2Data();
    ~TSamuraiFDC2Data();

  private:
    std::vector<int> fFDC2_DetectorNbr;
    std::vector<int> fFDC2_LayerNbr;
    std::vector<int> fFDC2_WireNbr;
    std::vector<double> fFDC2_Time;
    std::vector<int> fFDC2_Edge;
  
  public:
    void Clear();
    void Print();
    void Clear(const Option_t*) {};
    void Dump() const{};
  
  public:
    void SetData(const int& Det, const int& Layer, const int& Wire, const double& Time, const int& Edge);
    unsigned int Mult(){return fFDC2_DetectorNbr.size();};
    unsigned int MultLayer(unsigned int det , unsigned int layer, int edge=-1);
    std::vector<int> GetWire(unsigned int det , unsigned int layer);
    int const GetDetectorNbr(const unsigned int& i){return fFDC2_DetectorNbr[i];};
    int const GetLayerNbr(const unsigned int& i){return fFDC2_LayerNbr[i];};
    int const GetWireNbr(const unsigned int& i){return fFDC2_WireNbr[i];};
    double const GetTime(const unsigned int& i){return fFDC2_Time[i];};
    int const GetEdge(const unsigned int& i){return fFDC2_Edge[i];};

    ClassDef(TSamuraiFDC2Data,1); 
};

#endif
