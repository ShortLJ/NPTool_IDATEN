#ifndef TSamuraiBDCData_H
#define TSamuraiBDCData_H
#include "TObject.h"
#include <vector>

class TSamuraiBDCData: public TObject{
  public:
    TSamuraiBDCData();
    ~TSamuraiBDCData();

  private:
    std::vector<int> fBDC_DetectorNbr;
    std::vector<int> fBDC_LayerNbr;
    std::vector<int> fBDC_WireNbr;
    std::vector<double> fBDC_Time;
    std::vector<int> fBDC_Edge;
  
  public:
    void Clear();
    void Print();
    void Clear(const Option_t*) {};
    void Dump() const{};
  
  public:
    void SetData(const int& Det, const int& Layer, const int& Wire, const double& Time, const int& Edge);
    unsigned int Mult(){return fBDC_DetectorNbr.size();};
    unsigned int MultLayer(unsigned int det , unsigned int layer, int edge=-1);
    std::vector<int> GetWire(unsigned int det , unsigned int layer);
    int const GetDetectorNbr(const unsigned int& i){return fBDC_DetectorNbr[i];};
    int const GetLayerNbr(const unsigned int& i){return fBDC_LayerNbr[i];};
    int const GetWireNbr(const unsigned int& i){return fBDC_WireNbr[i];};
    double const GetTime(const unsigned int& i){return fBDC_Time[i];};
    int const GetEdge(const unsigned int& i){return fBDC_Edge[i];};

    ClassDef(TSamuraiBDCData,1); 
};

#endif
