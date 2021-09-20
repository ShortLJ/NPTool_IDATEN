#ifndef TBigRIPSReco_H
#define TBigRIPSReco_H
#include "TObject.h"
#include <vector>
#include <iostream>
#include "NPPhysicalConstants.h"

using namespace NPUNITS;

class TBigRIPSReco: public TObject{

  public:
    TBigRIPSReco();
    ~TBigRIPSReco();

  private:
    double aoq;
    double beta;
    double delta;
    double brho1; // upstream brho if twofold
    double brho2; // downstream brho if twofold
    double brho;
    double angle;
    double z;
  public:
    void Clear();
    void Init();
    void Print();
    void RecBrho(std::vector<double>,std::vector<double>,std::vector<std::vector<double>>,double);
    void RecAoqOneFold(double, double);
    void RecAoqTwoFold(double, double, double, int);
    void RecZet(double,double, std::vector<double>);
    //const double c_mm_ns = NPUNITS::c_light; //  299.779 mm/ns
    //const double amu_c2 = NPUNITS::amu_c2;  // 931.494  MeV
  
  public:
    void SetAoq(double value){aoq=value;}
    void SetBeta(double value){beta=value;}
    void SetDelta(double value){delta=value;}
    void SetBrho(double value){brho=value;}
    void SetBrho1(double value){brho1=value;}
    void SetBrho2(double value){brho2=value;}
    void SetZ(double value){z=value;}
    
    double GetAoq(){return aoq;}
    double GetBeta(){return beta;}
    double GetDelta(){return delta;}
    double GetBrho1(){return brho1;}
    double GetBrho2(){return brho2;}
    double GetBrho(){return brho;}
    double GetAngle(){return angle;}
    double GetZ(){return z;}

    ClassDef(TBigRIPSReco,1); 
};

#endif
