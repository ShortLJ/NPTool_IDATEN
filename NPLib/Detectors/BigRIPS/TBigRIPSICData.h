#ifndef TBigRIPSICData_H
#define TBigRIPSICData_H
#include "TObject.h"
#include <vector>
class TBigRIPSICData: public TObject{
  public:
    TBigRIPSICData();
    ~TBigRIPSICData();

  private:

    std::vector<int> fIC_E;
    std::vector<int> fIC_E_ID;
    std::vector<int> fIC_E_Layer;

    std::vector<int> fIC_T;
    std::vector<int> fIC_T_ID;
    std::vector<int> fIC_T_Layer;

  public:
    void Clear();
    void Print();
    void Clear(const Option_t*) {};
    void Dump() const{};
  
  public:
    void SetE(const int& E, const int& ID, const int& Layer){fIC_E.push_back(E);fIC_E_ID.push_back(ID);fIC_E_Layer.push_back(Layer);}
    void SetT(const int& T, const int& ID, const int& Layer){fIC_T.push_back(T);fIC_T_ID.push_back(ID);fIC_T_Layer.push_back(Layer);}

    int GetEMult() const {return fIC_E_ID.size();}
    int const GetE(const unsigned int& i){return fIC_E[i];};
    int const GetEID(const unsigned int& i){return fIC_E_ID[i];};
    int const GetELayer(const unsigned int& i){return fIC_E_Layer[i];};

    int GetTMult() const {return fIC_T_ID.size();}
    int const GetT(const unsigned int& i){return fIC_T[i];};
    int const GetTID(const unsigned int& i){return fIC_T_ID[i];};
    int const GetTLayer(const unsigned int& i){return fIC_T_Layer[i];};

    ClassDef(TBigRIPSICData,1); 
};

#endif
