/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta   contact address: matta@lpcaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : October 2020                                            *
 * Last update     : October 2020                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class have all the method needed to analyse Drift Chambers          *
 *****************************************************************************/

#include"NPDCReconstruction.h"
#include "Math/Factory.h"
#include "TError.h"
#include "TGraph.h"
using namespace std;
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
DCReconstruction::DCReconstruction(){
  m_min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  m_func=ROOT::Math::Functor(this,&NPL::DCReconstruction::SumD,2); 
  m_min->SetFunction(m_func);
  m_min->SetPrintLevel(0);

  // this avoid error 
  gErrorIgnoreLevel = kError;
  //m_min->SetMaxFunctionCalls(1000); 
  //m_min->SetMaxIterations(1000);
  //m_min->SetTolerance(1);
  //m_min->SetPrecision(1e-10);
}
////////////////////////////////////////////////////////////////////////////////
DCReconstruction::~DCReconstruction(){
  delete m_min;
}

////////////////////////////////////////////////////////////////////////////////
double DCReconstruction::BuildTrack2D(const vector<double>& X,const vector<double>& Z,const vector<double>& R,double& X0,double& X100,double& a, double& b ){
  fitX=&X;
  fitZ=&Z;
  fitR=&R;
  // assume all X,Z,R of same size
  sizeX = X.size();
  // Define the starting point of the fit: a straight line passing through the 
  // the first and last wire
  // z = ax+b -> x=(z-b)/a
  unsigned int i = 1;
  ai=1/0.;
  while(isinf(ai)&&i!=sizeX){
    ai = (Z[sizeX-i]-Z[0])/(X[sizeX-i]-X[0]);
    bi = Z[0]-ai*(X[0]);
    i++;
  }
  if(isinf(ai)){ // then there is no two point in different layer
    a=-10000;
    b=-10000;
    X0=-10000;
    X100=-10000;
    return 10000;
  }

  m_min->Clear(); 
  m_min->SetVariable(0,"a",ai,1);
  m_min->SetVariable(1,"b",bi,1);

  // Perform minimisation
  m_min->Minimize(); 
  //std::cout << "EDM:" <<  m_min->Edm() << std::endl;
  // access set of parameter that produce the minimum
  const double *xs = m_min->X();
  a=xs[0];
  b=xs[1];
  X0=-b/a;
  X100=(100-b)/a;
  return m_min->MinValue() ;
}

////////////////////////////////////////////////////////////////////////////////
void DCReconstruction::ResolvePlane(const TVector3& L,const double& ThetaU ,const TVector3& H, const double& ThetaV, TVector3& PosXY){
  // direction of U and V wire
  TVector3 u(0,1,0);
  u.RotateZ(ThetaU);

  TVector3 v(0,1,0);
  v.RotateZ(ThetaV);


  // Compute the coeff of the two line of vecotr u (v) going through H (L)
  // dv : y = av*x+bv
  av = v.Y()/v.X();
  bv = H.Y() - av*H.X();

  // du : y = au*x+bu
  au = u.Y()/u.X();
  bu = L.Y() - au*L.X();

  // We look for M(xM, yM) that intersect du and dv:
  if(!isinf(au) && !isinf(av)){ // au and av are not inf, i.e. not vertical line
    xM = (bv-bu)/(au-av);
    yM = au*xM+bu;
  }
  else if(isinf(av)){// av is inf, so v is along Y axis, H is direct measure of X
    xM = H.X();
    yM = au*xM+bu;
  }
  else if (isinf(au)){//au is inf, so u is along Y axis, L is direct measure of X
    xM = L.X();
    yM = av*xM+bv;
  }
  else{ // all is lost
    xM=-10000;
    yM=-10000;
  }
  PosXY.SetXYZ(xM,yM,0);
}


////////////////////////////////////////////////////////////////////////////////
double DCReconstruction::SumD(const double* parameter ){
  // Compute the sum P of the distance between the circle and the track
  P = 0;
  a = parameter[0];
  b = parameter[1];
  ab= a*b;
  a2=a*a;
  for(unsigned int i = 0 ; i < sizeX ; i++){
    c = (*fitX)[i];
    d = (*fitZ)[i];
    r = (*fitR)[i];
    x = (a*d-ab+c)/(1+a2);
    z = a*x+b;
    p= (x-c)*(x-c)+(z-d)*(z-d)-r*r;
    // numerical trick to have a smooth derivative instead of using abs
    P+= sqrt(p*p+0.1);
  }

  // return normalized power
  return P/sizeX;
}


////////////////////////////////////////////////////////////////////////////////
TGraph* DCReconstruction::Scan(double a, double b, int tovary, double minV, double maxV){
  vector<double> x,y;
  unsigned int sizeT=1000;
  double step = (maxV-minV)/sizeT;
  double p[2]={a,b};
  for(unsigned int i = 0 ; i < sizeT ; i++){
    if(!tovary){
      p[0]=minV+step*i; 
      x.push_back(p[0]);
      y.push_back(SumD(p));
    }

    else{
      p[1]=minV+step*i; 
      x.push_back(p[1]);
      y.push_back(SumD(p));
    }
  }

  TGraph* g = new TGraph(x.size(),&x[0],&y[0]);
  return g;
}
