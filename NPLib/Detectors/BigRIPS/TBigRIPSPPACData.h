#ifndef TBigRIPSPPACData_H
#define TBigRIPSPPACData_H
#include "TObject.h"
#include <vector>
class TBigRIPSPPACData: public TObject{
  public:
    TBigRIPSPPACData();
    ~TBigRIPSPPACData();

  private:
    //std::vector<int> fPPAC_FP; //focal plane
    //std::vector<int> fPPAC_ID; //specific ID
    //std::vector<int> fPPAC_Edge;

    std::vector<int> fPPAC_TX1;
    std::vector<int> fPPAC_TX2;
    std::vector<int> fPPAC_TY1;
    std::vector<int> fPPAC_TY2;
    std::vector<int> fPPAC_TA;
    std::vector<int> fPPAC_QX1;
    std::vector<int> fPPAC_QX2;
    std::vector<int> fPPAC_QY1;
    std::vector<int> fPPAC_QY2;
    std::vector<int> fPPAC_QA;

    std::vector<int> fPPAC_TX1_ID;
    std::vector<int> fPPAC_TX2_ID;
    std::vector<int> fPPAC_TY1_ID;
    std::vector<int> fPPAC_TY2_ID;
    std::vector<int> fPPAC_TA_ID;
    std::vector<int> fPPAC_QX1_ID;
    std::vector<int> fPPAC_QX2_ID;
    std::vector<int> fPPAC_QY1_ID;
    std::vector<int> fPPAC_QY2_ID;
    std::vector<int> fPPAC_QA_ID;

/*
    enum ChannelType{
        TX1,
        TX2,
        TY1,
        TY2,
        TA,
        QX1,
        QX2,
        QY1,
        QY2,
        QA
    };//should match with the Channel Type define in Unpacker 
  */

  public:
    void Clear();
    void Print();
    void Clear(const Option_t*) {};
    void Dump() const{};
  
  public:
    //void SetData(const int& FP, const int& ID, const int& Edge, const int& Value,const int& VariableType);

    //unsigned int Mult(){return fPPAC_DetectorNbr.size();};
    //unsigned int MultLayer(unsigned int det , unsigned int layer, int edge=-1);
    //int const GetID(const unsigned int& i){return fPPAC_ID[i];};
    //int const GetFP(const unsigned int& i){return fPPAC_FP[i];};
    //int const GetEdge(const unsigned int& i){return fPPAC_Edge[i];};
    /*
    double const GetTX2(const unsigned int& i){return fPPAC_TX2[i];};
    double const GetTY1(const unsigned int& i){return fPPAC_TY1[i];};
    double const GetTY2(const unsigned int& i){return fPPAC_TY2[i];};
    double const GetTA(const unsigned int& i){return fPPAC_TA[i];};
    double const GetQX1(const unsigned int& i){return fPPAC_QX1[i];};
    double const GetQX2(const unsigned int& i){return fPPAC_QX2[i];};
    double const GetQY1(const unsigned int& i){return fPPAC_QY1[i];};
    double const GetQY2(const unsigned int& i){return fPPAC_QY2[i];};
    double const GetQA(const unsigned int& i){return fPPAC_QA[i];};
*/

    //void SetID(const int& i){fPPAC_ID.push_back(i);}
    //void SetFP(const int& i){fPPAC_FP.push_back(i);}
    //void SetEdge(const int& i){fPPAC_Edge.push_back(i);}
    void SetTX1(const int& T, const int& ID){fPPAC_TX1.push_back(T);fPPAC_TX1_ID.push_back(ID);}
    void SetTX2(const int& T, const int& ID){fPPAC_TX2.push_back(T);fPPAC_TX2_ID.push_back(ID);}
    void SetTY1(const int& T, const int& ID){fPPAC_TY1.push_back(T);fPPAC_TY1_ID.push_back(ID);}
    void SetTY2(const int& T, const int& ID){fPPAC_TY2.push_back(T);fPPAC_TY2_ID.push_back(ID);}
    void SetTA(const int& T, const int& ID){fPPAC_TA.push_back(T);fPPAC_TA_ID.push_back(ID);}

    void SetQX1(const int& Q, const int& ID){fPPAC_QX1.push_back(Q);fPPAC_TX1_ID.push_back(ID);}
    void SetQX2(const int& Q, const int& ID){fPPAC_QX2.push_back(Q);fPPAC_TX2_ID.push_back(ID);}
    void SetQY1(const int& Q, const int& ID){fPPAC_QY1.push_back(Q);fPPAC_TY1_ID.push_back(ID);}
    void SetQY2(const int& Q, const int& ID){fPPAC_QY2.push_back(Q);fPPAC_TY2_ID.push_back(ID);}
    void SetQA(const int& Q, const int& ID){fPPAC_QA.push_back(Q);fPPAC_TA_ID.push_back(ID);}


    int GetTX1Mult() const {return fPPAC_TX1_ID.size();}
    int const GetTX1(const unsigned int& i){return fPPAC_TX1[i];};
    int const GetTX1ID(const unsigned int& i){return fPPAC_TX1_ID[i];};

    int GetTX2Mult() const {return fPPAC_TX2_ID.size();}
    int const GetTX2(const unsigned int& i){return fPPAC_TX2[i];};
    int const GetTX2ID(const unsigned int& i){return fPPAC_TX2_ID[i];};

    int GetTY1Mult() const {return fPPAC_TY1_ID.size();}
    int const GetTY1(const unsigned int& i){return fPPAC_TY1[i];};
    int const GetTY1ID(const unsigned int& i){return fPPAC_TY1_ID[i];};

    int GetTY2Mult() const {return fPPAC_TY2_ID.size();}
    int const GetTY2(const unsigned int& i){return fPPAC_TY2[i];};
    int const GetTY2ID(const unsigned int& i){return fPPAC_TY2_ID[i];};

    int GetTAMult() const {return fPPAC_TA_ID.size();}
    int const GetTA(const unsigned int& i){return fPPAC_TA[i];};
    int const GetTAID(const unsigned int& i){return fPPAC_TA_ID[i];};

    int GetQX1Mult() const {return fPPAC_QX1_ID.size();}
    int const GetQX1(const unsigned int& i){return fPPAC_QX1[i];};
    int const GetQX1ID(const unsigned int& i){return fPPAC_QX1_ID[i];};

    int GetQX2Mult() const {return fPPAC_QX2_ID.size();}
    int const GetQX2(const unsigned int& i){return fPPAC_QX2[i];};
    int const GetQX2ID(const unsigned int& i){return fPPAC_QX2_ID[i];};

    int GetQY1Mult() const {return fPPAC_QY1_ID.size();}
    int const GetQY1(const unsigned int& i){return fPPAC_QY1[i];};
    int const GetQY1ID(const unsigned int& i){return fPPAC_QY1_ID[i];};

    int GetQY2Mult() const {return fPPAC_QY2_ID.size();}
    int const GetQY2(const unsigned int& i){return fPPAC_QY2[i];};
    int const GetQY2ID(const unsigned int& i){return fPPAC_QY2_ID[i];};

    int GetQAMult() const {return fPPAC_QA_ID.size();}
    int const GetQA(const unsigned int& i){return fPPAC_QA[i];};
    int const GetQAID(const unsigned int& i){return fPPAC_QA_ID[i];};

    ClassDef(TBigRIPSPPACData,1); 
};

#endif
