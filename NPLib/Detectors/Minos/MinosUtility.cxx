#include "MinosUtility.h"
#include "Math/Factory.h"
#include "TROOT.h"
#include <algorithm>
#include <fstream>
#include <vector>
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
MinosUtility::MinosUtility(){
  ROOT::EnableThreadSafety();
  // Setting up for signal fitting 
  m_signal_min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  m_signal_func=ROOT::Math::Functor(this,&NPL::MinosUtility::signal_chi2,2); 
  m_signal_min->SetFunction(m_signal_func);
  //30 ns sample (33MHz FEM clock) and 333.9 shaping time on preamp, fit with 1/10 of the point
  SetParameters(30,333.9,250,10);
  m_counter=0;

}

////////////////////////////////////////////////////////////////////////////////
double MinosUtility::signal_shape(const double& x, const double* p){
  double a = (x-p[1])/m_Tau;
  if(x > p[1] && x < 512.) 
    return (p[0]*exp(-3.*a)*sin(a)*(a*a*a) + m_Baseline);
  else 
    return m_Baseline;
}

////////////////////////////////////////////////////////////////////////////////
double MinosUtility::signal_chi2(const double* p){
  double chi2=0;
  unsigned int counter=0;
  double diff,expected;

  for(auto it = m_vTQ.begin(); it!=m_vTQ.end() ; it++ ){
    expected = signal_shape(it->first,p);
    diff  = it->second - expected ; 
    chi2 += diff*diff; 
    counter++;
  }
  return chi2/counter;
}
////////////////////////////////////////////////////////////////////////////////
void MinosUtility::sample_signal(){
  m_vTQ.clear();
  for(unsigned int i = m_minbin; i < m_maxbin ; i+= m_Sampling ){
    m_vTQ.push_back(std::make_pair((*m_fitSignalT)[i],(*m_fitSignalQ)[i]));
  }
}
////////////////////////////////////////////////////////////////////////////////
void MinosUtility::find_maxQ(){
  m_qmax=0; m_qmaxbin=0;
  for(unsigned int i = 0 ; i < m_signal_size ; i++ ){
    if((*m_fitSignalQ)[i]>m_qmax){
      m_qmax=(*m_fitSignalQ)[i];
      m_qmaxbin=i;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
double MinosUtility::Calibrate(const std::vector<unsigned short>* T,const std::vector<unsigned short>* Q, const unsigned int& i , double& time, double& charge){
  // setting up for the minisation
  m_fitSignalT=T;
  m_fitSignalQ=Q;
  m_signal_size=m_fitSignalT->size();   
  find_maxQ();
  if(m_qmaxbin<20){
    time=-1;    
    charge=-1;
    return -10000;
  }

  sample_signal();
  m_guess_t0_bin = m_qmaxbin-20;
  m_minbin = m_guess_t0_bin;
  m_maxbin = std::min(m_guess_t0_bin+40,m_signal_size);
  m_signal_min->Clear();
  m_signal_min->SetLimitedVariable(0,"A",m_qmax*10,100,0,m_qmax*20);
  m_signal_min->SetLimitedVariable(1,"t0",(*m_fitSignalT)[m_guess_t0_bin],1,0,(*m_fitSignalT)[m_qmaxbin]);

  // Perform minimisation
  m_signal_min->Minimize(); 
  // access set of parameter that produce the minimum
  const double* xs = m_signal_min->X();
  charge= xs[0];
  time=xs[1];

  return m_signal_min->MinValue() ;

}
