#if __cplusplus > 199711L // require c++11 
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

#include"NPDCReconstructionMT.h"

// ROOT
#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "TError.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TROOT.h"
using namespace std;
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
DCReconstructionMT::DCReconstructionMT(unsigned int number_thread){
  ROOT::EnableThreadSafety();
  m_nbr_thread= number_thread;
  // force loading of the minimizer plugin ahead
  ROOT::Math::Minimizer* mini=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  delete mini;
}
////////////////////////////////////////////////////////////////////////////////
DCReconstructionMT::~DCReconstructionMT(){
  StopThread();
}

////////////////////////////////////////////////////////////////////////////////
double DCReconstructionMT::GetResults(unsigned int ID,double& X0,double& X100,double& a, double& b){
  if(m_X0.find(ID)!=m_X0.end()){
    X0=m_X0[ID]; 
    X100=m_X100[ID]; 
    a=m_a[ID];
    b=m_b[ID];
    return m_minimum[ID];
  }

  return -1;
}
////////////////////////////////////////////////////////////////////////////////
void DCReconstructionMT::AddPlan(unsigned int ID, const vector<double>& X, const vector<double>& Z, const vector<double>& R){
  // Select a free thread
  unsigned sizeR = m_Ready.size();
  unsigned int free_thread;
  bool found_thread=false;
  while(!found_thread){
    for(unsigned int i = 0 ; i < sizeR ; i++){
      if(!m_Ready[i]){
        free_thread=i;
        found_thread=true;
        break;
      }
    }
  }

  fitX[free_thread]=&X;
  fitZ[free_thread]=&Z;
  fitR[free_thread]=&R;
  m_uid[free_thread]=ID;
  // assume all X,Z,R of same size
  sizeX[free_thread] = X.size();
  m_Ready[free_thread]=true;
  return;
}
////////////////////////////////////////////////////////////////////////////////
void DCReconstructionMT::BuildTrack2D(){
  while(!IsDone())
    std::this_thread::yield();
  return;
}

////////////////////////////////////////////////////////////////////////////////
bool DCReconstructionMT::IsDone(){
  std::vector<bool>::iterator begin = m_Ready.begin(); 
  std::vector<bool>::iterator end = m_Ready.end();
  for(auto it = begin ; it!=end ; it++){
    if((*it))
      return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void DCReconstructionMT::ResolvePlane(const TVector3& L,const double& ThetaU ,const TVector3& H, const double& ThetaV, TVector3& PosXY){
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
double DCReconstructionMT::SumD(const double* parameter ){
  // Compute the sum P of the distance between the circle and the track
  double P = 0;
  double a = parameter[0];
  double b = parameter[1];
  double ab= a*b;
  double a2=a*a;
  unsigned int id = parameter[2];
  unsigned int size =  sizeX[id];
  const std::vector<double>* X=fitX[id];
  const std::vector<double>* Z=fitZ[id];
  const std::vector<double>* R=fitR[id];

  double c,d,r,x,z,p;
  for(unsigned int i = 0 ; i < size ; i++){
    c = (*X)[i];
    d = (*Z)[i];
    r = (*R)[i];
    x = (a*d-ab+c)/(1+a2);
    z = a*x+b;
    p= (x-c)*(x-c)+(z-d)*(z-d)-r*r;
    // numerical trick to have a smooth derivative instead of using abs
    P+= sqrt(p*p+0.1);
  }

  // return normalized power
  return P/size;
}


////////////////////////////////////////////////////////////////////////////////
TGraph* DCReconstructionMT::Scan(double a, double b, int tovary, double minV, double maxV){
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
////////////////////////////////////////////////////////////////////////////////
void NPL::DCReconstructionMT::InitThreadPool(){

  // this avoid error printout during fitting
  gErrorIgnoreLevel = kError;

  StopThread();
  m_ThreadPool.clear();
  m_Ready.clear();
  m_Ready.resize(m_nbr_thread,false);
  for (unsigned int i=0; i < m_nbr_thread; i++) { 
    // Register minimiser for futur deletion
    m_ThreadPool.push_back( std::thread(&NPL::DCReconstructionMT::StartThread,this,i) );
  }
  m_stop = false;
  for(auto& th: m_ThreadPool){
    th.detach();
  }
}
////////////////////////////////////////////////////////////////////////////////
void DCReconstructionMT::StartThread(unsigned int id){
  // usefull variable
  double ai,bi;
  unsigned int uid;
  const double* xs;
  // create the functor 
  // each threads needs its own or the minisation is not thread safe 
  ROOT::Math::Functor* func= new ROOT::Math::Functor(this,&NPL::DCReconstructionMT::SumD,3); 
  //Create the minimiser (deleted by the thread)
  m_mtx.lock();
  ROOT::Math::Minimizer* mini=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  m_mtx.unlock();

  mini->SetFunction(*func);
  mini->SetPrintLevel(0);

  // Let the main thread start
  std::this_thread::sleep_for(std::chrono::milliseconds(250));
  while(true){
    // Do the job if possible
    if(m_Ready[id]){
      // Define the starting point of the fit: a straight line passing through the 
      // the first and last wire
      // z = ax+b -> x=(z-b)/a
        ai = ((*fitZ[id])[sizeX[id]-1]-(*fitZ[id])[0])/((*fitX[id])[sizeX[id]-1]-(*fitX[id])[0]+(*fitR[id])[sizeX[id]-1]-(*fitR[id])[0]);
        bi = (*fitZ[id])[0]-ai*((*fitX[id])[0]+(*fitR[id])[0]);

      if(isinf(ai)){ // then there is no two point in different layer
        m_a[uid]=-10000;
        m_b[uid]=-10000;
        m_X0[uid]=-10000;
        m_X100[uid]=-10000;
        m_minimum[uid] = 10000;
      }

      else{
        mini->Clear(); 
        mini->SetVariable(0,"a",ai,1);
        mini->SetVariable(1,"b",bi,1);
        mini->SetFixedVariable(2,"id",id);
        // Perform minimisation
        mini->Minimize(); 

        // access set of parameter that produce the minimum
        xs = mini->X();
        uid = m_uid[id]; 
        m_a[uid]=xs[0];
        m_b[uid]=xs[1];
        m_X0[uid]=-m_b[uid]/m_a[uid];
        m_X100[uid]=(100-m_b[uid])/m_a[uid];
        m_minimum[uid] = mini->MinValue();
      }
      // notify main thread job is done
      m_mtx.lock();// make sure no other thread is reading/writing to the map
      m_Ready[id].flip();
      m_mtx.unlock();
      // Let other thread move up in the queu
      std::this_thread::yield();
    }
    else{
      std::this_thread::yield();
    }
    // Return if stopped
    if(m_stop){
      delete mini;
      delete func;
      return;
    }
  }   
}
////////////////////////////////////////////////////////////////////////////////
void DCReconstructionMT::StopThread(){
  // make sure the last thread are schedule before stopping;
  std::this_thread::yield();
  m_stop=true;
  std::this_thread::yield();
}
#endif
