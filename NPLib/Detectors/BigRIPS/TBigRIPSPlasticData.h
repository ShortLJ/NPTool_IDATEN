#ifndef TBigRIPSPlasticData_H
#define TBigRIPSPlasticData_H
#include "TObject.h"
#include <vector>
class TBigRIPSPlasticData: public TObject{
  public:
    TBigRIPSPlasticData();
    ~TBigRIPSPlasticData();

  private:
    std::vector<int> fPlastic_TL;
    std::vector<int> fPlastic_TR;
    std::vector<int> fPlastic_QL;
    std::vector<int> fPlastic_QR;

    std::vector<int> fPlastic_TL_ID;
    std::vector<int> fPlastic_TR_ID;
    std::vector<int> fPlastic_QL_ID;
    std::vector<int> fPlastic_QR_ID;

  public:
    void Clear();
    void Print();
    void Clear(const Option_t*) {};
    void Dump() const{};
  
  public:
    void SetTL(const int& T, const int& ID){fPlastic_TL.push_back(T);fPlastic_TL_ID.push_back(ID);}
    void SetTR(const int& T, const int& ID){fPlastic_TR.push_back(T);fPlastic_TR_ID.push_back(ID);}

    void SetQL(const int& Q, const int& ID){fPlastic_QL.push_back(Q);fPlastic_QL_ID.push_back(ID);}
    void SetQR(const int& Q, const int& ID){fPlastic_QR.push_back(Q);fPlastic_QR_ID.push_back(ID);}

    int GetTLMult() const {return fPlastic_TL_ID.size();}
    int const GetTL(const unsigned int& i){return fPlastic_TL[i];};
    int const GetTLID(const unsigned int& i){return fPlastic_TL_ID[i];};

    int GetTRMult() const {return fPlastic_TR_ID.size();}
    int const GetTR(const unsigned int& i){return fPlastic_TR[i];};
    int const GetTRID(const unsigned int& i){return fPlastic_TR_ID[i];};

    int GetQLMult() const {return fPlastic_QL_ID.size();}
    int const GetQL(const unsigned int& i){return fPlastic_QL[i];};
    int const GetQLID(const unsigned int& i){return fPlastic_QL_ID[i];};

    int GetQRMult() const {return fPlastic_QR_ID.size();}
    int const GetQR(const unsigned int& i){return fPlastic_QR[i];};
    int const GetQRID(const unsigned int& i){return fPlastic_QR_ID[i];};

    ClassDef(TBigRIPSPlasticData,1); 
};

#endif
