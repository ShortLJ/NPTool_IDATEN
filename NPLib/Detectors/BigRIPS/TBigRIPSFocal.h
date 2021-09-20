#ifndef TBigRIPSFocal_H
#define TBigRIPSFocal_H
#include "TObject.h"
#include <vector>
#include <iostream>
class TBigRIPSFocal: public TObject{

  public:
    TBigRIPSFocal();
    ~TBigRIPSFocal();

  private:
    std::vector<std::vector<double>> FPTrack  ;

  public:
    void Clear();
    void Init(int);
    void Print(int);
  
  public:
    void SetFPTrack(int FPNbr, std::vector<double> track){
            FPTrack[FPNbr]=track;
    };
    void SetFPTrack(int FPNbr, double x, double a, double y, double b){
            FPTrack[FPNbr][0]=x;
            FPTrack[FPNbr][1]=a;
            FPTrack[FPNbr][2]=y;
            FPTrack[FPNbr][3]=b;
    };
    std::vector<double> GetFPTrack(int FPNbr) {return FPTrack[FPNbr];}

    ClassDef(TBigRIPSFocal,1); 
};

#endif
