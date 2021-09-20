#include "TBigRIPSReco.h"
#include <cmath>
TBigRIPSReco::TBigRIPSReco(){Init();}
TBigRIPSReco::~TBigRIPSReco(){Clear();}


////////////////////////////////////////////////////////////////////////////////
void TBigRIPSReco::Clear(){
}
////////////////////////////////////////////////////////////////////////////////
void TBigRIPSReco::Init(){
    aoq=-9999;
    beta=-9999;
    delta=-9999;
    brho1=-9999;
    brho2=-9999;
    brho=-9999;
}

////////////////////////////////////////////////////////////////////////////////

void TBigRIPSReco::RecBrho(std::vector<double> RecFPUpstream,std::vector<double> RecFPDownstream, std::vector<std::vector<double>> matrix,double BrhoCentral){

    if(matrix.size()!=6) {
        std::cout << "MATRIX used for Brho Rec HAS NOT 6 rows" << std::endl;
        return;
    }else{
        for(int i=0;i<6;i++){
            if(matrix[i].size()!=5){
                 std::cout << "MATRIX used for Brho Rec. HAS NOT 5 columns" << std::endl;
                 std::cout << "but" << matrix[i].size()<<std::endl;
                 return;
            }
        }
    }

    if(RecFPUpstream[0]==-9999 || 
       RecFPUpstream[1]==-9999 || 
       RecFPDownstream[0]==-9999)
        return;

    double x_x = matrix[0][0];
    double a_x = matrix[0][1];
    double x_a = matrix[1][0];
    double x_d = matrix[5][0];
    double a_a = matrix[1][1];
    double a_d = matrix[5][1];

    double x_up = RecFPUpstream[0]; // downst. X
    double a_up = RecFPUpstream[1]; // downst. X
    double x_down = RecFPDownstream[0]; // downst. X
    double a_down = RecFPDownstream[1]; // downst. A

    delta = (a_a * x_down - x_a * a_down - x_up) / (x_d * a_a - x_a * a_d);
    angle = (a_a * x_down - x_a * a_down - x_up) / (x_d * a_a - x_a * a_d);

    brho = BrhoCentral*(1.0+delta*0.01);
}


////////////////////////////////////////////////////////////////////////////////
void TBigRIPSReco::RecAoqOneFold(double tof, double length){
    beta = length /(tof * c_light);
    double gamma = 1./sqrt(1 - beta*beta);
    aoq = (brho * c_light) / (amu_c2 * beta * gamma);
}

////////////////////////////////////////////////////////////////////////////////
void TBigRIPSReco::RecAoqTwoFold(double tof, double length1, double length2, int useBeta){

    if(!(length1>0 && length2>0)){
        std::cout << "Length1 or Length2 is not positive (or both)" << std::endl;
        return;
    }
    if(!(brho1>0 && brho2>0)){
        //std::cout << "brho1 or brho2 was not properly initialized/setup" << std::endl;
        //std::cout << "brho1:" <<brho1<< std::endl;
        //std::cout << "brho2:" <<brho2<< std::endl;
        return;
    }

    double alpha  = brho2 / brho1;
    double a1     = sqrt(alpha * alpha * c_squared * tof * tof
          + (pow(alpha,4) - alpha * alpha) * length1 * length1
          + (1 - alpha*alpha) * length2 * length2);

  
    double rbeta1 = ( a1 * length1 + length2 * c_light * tof ) / 
                    ( a1 * c_light * tof + (1 - alpha * alpha) * length1 * length2);
    double  gamma1 = 1 / sqrt(1 - pow(rbeta1,2));

    double rbeta2 = ( a1 * length1 + length2 * c_light * tof ) /
                    ( c_squared * tof * tof + (alpha * alpha - 1) * length1 * length1);
    double gammab = 1 / sqrt(1 - pow(rbeta2,2));

    aoq = brho1 * c_light / amu_c2 / rbeta1 / gamma1; // should be same as brho2/beta2/gamma2
      
    if(useBeta == 2 ) beta = rbeta2;
    else if(useBeta == 1) beta = rbeta1;

    brho = brho2; 
/*
    std::cout << "-------------" << std::endl;
    std::cout << "brho1:" <<brho1<< std::endl;
    std::cout << "brho2:" <<brho2<< std::endl;
    std::cout << "brho:" <<brho<< std::endl;
    std::cout << "alpha:" <<alpha<< std::endl;
    std::cout << "rbeta1:" <<rbeta1<< std::endl;
    std::cout << "rbeta2:" <<rbeta2<< std::endl;
    std::cout << "aoq:" <<aoq<< std::endl;
  */  
}

////////////////////////////////////////////////////////////////////////////////
void TBigRIPSReco::RecZet(double dE,double ionpair, std::vector<double> coeff){
   if(coeff.size()<2){
       std::cout << "Zcoeff for Z reco not defined" << std::endl;
       return;
   }
   double de_v = log(ionpair*beta*beta) - log((1-beta*beta)) - beta*beta;
   z = coeff[0] * sqrt(dE/de_v)*beta + coeff[1];   
}
////////////////////////////////////////////////////////////////////////////////
void TBigRIPSReco::Print(){
    std::cout << "aoq: "  << aoq << std::endl;
    std::cout << "z: "  << z << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "delta: "<< delta << std::endl;
    std::cout << "brho: " << brho << std::endl;
}

ClassImp(TBigRIPSReco); 
