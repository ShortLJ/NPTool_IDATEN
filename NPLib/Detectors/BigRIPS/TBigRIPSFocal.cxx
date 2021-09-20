#include "TBigRIPSFocal.h"

TBigRIPSFocal::TBigRIPSFocal(){};
TBigRIPSFocal::~TBigRIPSFocal(){Clear();};


////////////////////////////////////////////////////////////////////////////////
void TBigRIPSFocal::Clear(){
    FPTrack.clear();
}
////////////////////////////////////////////////////////////////////////////////
void  TBigRIPSFocal::Init(int N_FP){
    std::vector<double> tmp;
    for(int j=0;j<4;j++) tmp.push_back(-9999);
    for(int i=0;i<N_FP;i++) FPTrack.push_back(tmp);
}
////////////////////////////////////////////////////////////////////////////////
void TBigRIPSFocal::Print(int FPNbr){
    std::cout << "----FP----:" << FPNbr <<std::endl;
    std::cout << "X: " << FPTrack[FPNbr][0] << std::endl;
    std::cout << "A: " << FPTrack[FPNbr][1] << std::endl;
    std::cout << "Y: " << FPTrack[FPNbr][2] << std::endl;
    std::cout << "B: " << FPTrack[FPNbr][3] << std::endl;

}
ClassImp(TBigRIPSFocal); 
