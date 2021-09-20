//----------! Developed by C. Santamaria / CEA Saclay !----------
//----------!      Version date :: 2014/11/12         !----------
// Tracking Functions
#include "Math/Vector3D.h"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include "TClonesArray.h"
#include <vector>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraph.h"
/* #include "Vector3D.h" */
// #include "/opt/root/v5.34.14/include/Math/Vector3D.h"
/* #include "/opt/cern/root/root_v5.34.36/math/genvector/inc/Math/Vector3D.h" */
/* #include "TMinosResult.h" */
#include "Tracking.h"

using namespace std;

/* TClonesArray data_result; */
/* TMinosResult *minosdata_result; */

ClassImp(NPL::Tracking)

NPL::Tracking::Tracking() 
{
// Constructor
}

NPL::Tracking::~Tracking()
{
// Destructor
}

// Fit function for E(T)
/* double NPL::Tracking::conv_fit(double *x, double *p) { */
/*     double val; */
/*     if(!(x[0]<p[1] || x[0]>512.)) val = p[0] * exp(-3.*(x[0]-p[1])/p[2])  * sin((x[0]-p[1])/p[2]) * pow((x[0]-p[1])/p[2], 3) + 250; */
/*     //else val = p[3]; */
/*     else val = 250; */
/*     return(val); */
/* } */

int NPL::Tracking::Hough_modified(vector<double> *x,vector<double> *y,vector<double> *q, vector<double> *x_out,vector<double> *y_out, vector<double> *q_out, vector<int> *ringbool) {
  double bint1=2.;
  double bint2=2.;
  int maxt = 360.;
  int mint = 0.;
  int nt1=(maxt-mint)/bint1;
  int nt2=(maxt-mint)/bint1;

  double PI = TMath::Pi();
  double Rint = 45.2;
  double Rext = 45.2 + 18*2.1;
  int filter_result = 0;

  static TH2F hp_xy1 = TH2F("hp_xy1","hp_xy1",nt1,mint,maxt,nt2,mint,maxt);
  static TH2F hpDiag_xy = TH2F("hpDiag_xy","hpDiag_xy",nt1,mint,maxt,nt2,mint,maxt);
  hp_xy1.Reset();
  hpDiag_xy.Reset();
  
  double max_xy;
  static vector<double> xTemp, yTemp, qTemp;
  xTemp.clear(); yTemp.clear(); qTemp.clear();

  double theta1, theta2, xt1, yt1, xt2, yt2;
  double line0=0., line1=0.;
  double delta=0., AA=0., BB=0., CC=0.;
  double maxtheta1=0., maxtheta2=0., xmax1=0., ymax1=0., xmax2=0., ymax2=0.;
  double par0=0., par1=0.;
  double r_mm=0.;
  int ringsum=0;
  bool maxfound = false;

  for(unsigned int i=0;i<x->size();i++){
    xTemp.push_back(x->at(i));
    yTemp.push_back(y->at(i));
    qTemp.push_back(q->at(i));

    //Loop of indices
    for(int j=0;j<nt1;j++){

      theta1 = (j+0.5)*bint1 + mint;
      xt1 = Rint * TMath::Cos(theta1*PI/180.);
      yt1 = Rint * TMath::Sin(theta1*PI/180.);
      line1 = (yt1 - y->at(i))/(xt1 - x->at(i));
      line0 = yt1 - xt1 * line1;
      AA = 1 + line1*line1;
      BB = 2*line0*line1;
      /* CC = line0*line0 - Rext; // Warning likely Error : Rext and not Rext^2 (Cyril) */
      CC = line0*line0 - Rext*Rext; // Warning likely Error : Rext and not Rext^2 (Cyril)

      delta = BB*BB - 4*AA*CC;

      if(delta>=0){
        xt2 = (-BB - sqrt(delta))/(2*AA);
        yt2 = line0 + line1*xt2;
        if(xt2<=0) theta2=  180 - asin(yt2/Rext)*180/PI;
        else if(xt2>0){
          if(yt2>0)	theta2=  asin(yt2/Rext)*180/PI;
          else if(yt2<=0)	      theta2=  360 + asin(yt2/Rext)*180/PI;
        }
        if( (xt1*x->at(i) + yt1*y->at(i))>=0 && (xt2*x->at(i) + yt2*y->at(i))>=0  && (xt1*xt2+yt1*yt2)>=0){
          hp_xy1.Fill(theta1,theta2);
          if(abs(theta1-theta2)<=10) hpDiag_xy.Fill(theta1,theta2);
        }
        else{
          if(delta!=0){
            xt2 = (-BB + sqrt(delta))/(2*AA);
            yt2 = line0 + line1*xt2;
            if(xt2<=0) theta2=  180 - asin(yt2/Rext)*180/PI;
            else if(xt2>0){
              if(yt2>0)	theta2=  asin(yt2/Rext)*180/PI;
              else if(yt2<=0)	      theta2=  360 + asin(yt2/Rext)*180/PI;
            }
            if( (xt1*x->at(i) + yt1*y->at(i))>=0 && (xt2*x->at(i) + yt2*y->at(i))>=0  && (xt1*xt2+yt1*yt2)>=0){
              hp_xy1.Fill(theta1,theta2);
              if(abs(theta1-theta2)<=10) hpDiag_xy.Fill(theta1,theta2);
            }
          }
        }
      } // end if delta>=0
    } // end loop on indices
  } // end loop on pads

  x->clear(); 
  y->clear();
  q->clear();

  if(hpDiag_xy.GetMaximum()>=10) max_xy = hpDiag_xy.GetMaximum();
  //		cout << "Max taken in diag... withh value=" << max_xy << endl;
  else max_xy = hp_xy1.GetMaximum();

  for(int ii=0; ii<nt1; ii++){
    if(maxfound ==true) break;
    for(int jj=0; jj<nt2; jj++){
      if(hp_xy1.GetBinContent(ii+1, jj+1) == max_xy){
        maxtheta1 = (ii+0.5)*bint1 + mint;
        maxtheta2 = (jj+0.5)*bint2 + mint;
        maxfound = true;
        //cout << "xy: theta max are " << maxtheta1 << " , " << maxtheta2 << endl;
      }
      if(maxfound ==true) break;
    }
  }

  xmax1 = Rint * TMath::Cos(maxtheta1*PI/180.);
  ymax1 = Rint * TMath::Sin(maxtheta1*PI/180.);
  xmax2 = Rext * TMath::Cos(maxtheta2*PI/180.);
  ymax2 = Rext * TMath::Sin(maxtheta2*PI/180.);

  // xy PEAK
  par1 = (ymax2-ymax1)/(xmax2-xmax1);
  par0 = (ymax1 - xmax1*par1);
  //Selection of x,y points IN the maxmean+/-1 found in Obertelli transform of xy plane
  static int xTempSize;
  xTempSize = xTemp.size();
  for(unsigned int i=0;i<xTempSize;i++){
    if( (abs(par1*xTemp[i]-yTemp[i]+par0)/sqrt(1+par1*par1))<= 6 && ((xmax1*xTemp[i] + ymax1*yTemp[i]) >= 0) && ((xmax2*xTemp[i] + ymax2*yTemp[i]) >= 0) && ((xmax1*xmax2 + ymax1*ymax2) >= 0)){
      //			hcnew_xy->Fill(xTemp[i],yTemp[i],qTemp[i]);
      x_out->push_back(xTemp[i]); 
      y_out->push_back(yTemp[i]);
      q_out->push_back(qTemp[i]);
      filter_result++;
      /* r_mm = sqrt(xTemp[i]*xTemp[i]+yTemp[i]*yTemp[i]); */
      /* if(r_mm<(45.2+5*2.1)) ringsum++; */  // commented by Cyril   
    }
    else{	
      x->push_back(xTemp[i]);
      y->push_back(yTemp[i]);
      q->push_back(qTemp[i]);
    }
  }
  for(int ip=0; ip<filter_result; ip++){
    /* if(ringsum>2) ringbool->push_back(1); */
    /* else ringbool->push_back(0); */
    ringbool->push_back(1);
  }
  return filter_result;
}

double NPL::Tracking::FitFunction(double *x, double *p) {
  double val=p[0]+p[1]*x[0];
  return(val);
}
void NPL::Tracking::FindStart(double pStart[4], double chi[2],  int fitStatus[2],TGraph *grxz, TGraph *gryz) {
  static double par1D[2];
  static TF1 *myfit1 = new TF1("myfit1","[0]+[1]*x", -100,500);
  myfit1->Clear();
  myfit1->SetParameter(0,0);
  myfit1->SetParameter(1,10);
  fitStatus[0] =0;
  grxz->Fit(myfit1,"RQM");
  chi[0]=myfit1->GetChisquare();
  par1D[0]=myfit1->GetParameter(0); 
  par1D[1]=myfit1->GetParameter(1);
  pStart[0]=par1D[0];
  pStart[1]=par1D[1];
  fitStatus[1] =0;
  gryz->Fit(myfit1,"RQM");
  chi[1]=myfit1->GetChisquare();
  par1D[0]=myfit1->GetParameter(0);
  par1D[1]=myfit1->GetParameter(1);
  pStart[2]=par1D[0];
  pStart[3]=par1D[1];
  //AC 07/12/14
}

// Calculation of the distance line-point
double NPL::Tracking::distance2(double x,double y,double z, double *p) {
  // distance line point is D= | (xp-x0) cross  ux |
  // where ux is direction of line and x0 is a point in the line (like t = 0)
  ROOT::Math::XYZVector xp(x,y,z); //point of the track
  ROOT::Math:: XYZVector x0(p[0], p[2], 0. );
  ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); //line
  ROOT::Math::XYZVector u = (x1-x0).Unit();
  double d2 = ((xp-x0).Cross(u)) .Mag2();
  return d2;
}

//void Hough_3D(TCanvas *c1, vector<double> *x,vector<double> *y,vector<double> *z,vector<double> *q, vector<double> *x_out,vector<double> *y_out,vector<double> *z_out,vector<double> *q_out) {
void NPL::Tracking::Hough_3D(vector<double> *x,vector<double> *y,vector<double> *z,vector<double> *q,vector<double> *x_out,vector<double> *y_out,vector<double> *z_out,vector<double> *q_out) {

  int nt_xy=180;
  int nt_xz=180;
  int nt_yz=180;
  int nr_xy=45;
  int nr_xz=300;
  int nr_yz=300;
  double bint_xy=2.;
  double bint_xz=2.;
  double bint_yz=2.;
  double binr_xy=3.;
  double binr_xz=3.;
  double binr_yz=3.;
  int nt,nr;
  double PI = TMath::Pi();

  double rho_xy,rho_xz,rho_yz;
  double theta_xy,theta_xz,theta_yz;
  
  static TH2F hp_xy("hp_xy","hp_xy",nt_xy,0,180,nr_xy,-1*nr_xy,nr_xy);
  static TH2F hp_xz("hp_xz","hp_xz",nt_xz,0,180,nr_xz,-1*nr_xz,nr_xz);
  static TH2F hp_yz("hp_yz","hp_yz",nt_yz,0,180,nr_yz,-1*nr_yz,nr_yz);
  hp_xy.Reset();
  hp_xz.Reset();
  hp_yz.Reset();

  //	int npeaks_xy, npeaks_xz, npeaks_yz;
  vector<double> thetapeaks_xy, rpeaks_xy, thetapeaks_xz, rpeaks_xz, thetapeaks_yz, rpeaks_yz;
  double max_xy, max_xz, max_yz;
  double rmean_xy=0, thetamean_xy=0, rmean_xz=0, thetamean_xz=0, rmean_yz=0, thetamean_yz=0;

  double r0_xy=0., r0_xz=0., r0_yz=0., rmin_xy=0., rmin_xz=0., rmin_yz=0., rmax_xy=0., rmax_xz=0., rmax_yz=0.;
  double tmin=0., tmax=0.;
  double rinf=0., rsup=0.;

  nt=nt_xy;
  nr=nr_xy;
  if(nt<nt_xz)nt=nt_xz;
  if(nr<nr_xz)nr=nr_xz;
  if(nt<nt_yz)nt=nt_yz;
  if(nr<nr_yz)nr=nr_yz;

  for(unsigned int i=0;i<x->size();i++)
  {
    //Fill coordinate space histograms for plots
    //Loop of indices and fill Histograms
    for(int j=0;j<nt;j++)
    {
      //xy
      theta_xy = j*180./nt_xy;
      rho_xy = x->at(i)*TMath::Cos(theta_xy*PI/180.)+y->at(i)*TMath::Sin(theta_xy*PI/180.);
      if(abs(theta_xy)<180. && abs(rho_xy)<nr_xy)
      {
        hp_xy.Fill(theta_xy,rho_xy);
      }

      //xz
      theta_xz = j*180./nt_xz;
      rho_xz = z->at(i)*TMath::Cos(theta_xz*PI/180.)+x->at(i)*TMath::Sin(theta_xz*PI/180.);
      if(abs(theta_xz)<180. && abs(rho_xz)<nr_xz)
      {
        hp_xz.Fill(theta_xz,rho_xz);
      }

      //yz
      theta_yz = j*180./nt_yz;
      rho_yz = z->at(i)*TMath::Cos(theta_yz*PI/180.)+y->at(i)*TMath::Sin(theta_yz*PI/180.);
      if(abs(theta_yz)<180. && abs(rho_yz)<nr_yz)
      {
        hp_yz.Fill(theta_yz,rho_yz);
      }
    }
  }

  max_xy = hp_xy.GetMaximum();
  max_xz = hp_xz.GetMaximum();
  max_yz = hp_yz.GetMaximum();

  for(int ii=0; ii<nt; ii++)
  {
    for(int jj=0; jj<nr; jj++)
    {
      if(hp_xy.GetBinContent(ii+1, jj+1) == max_xy && jj<nr_xy)
      {
        thetapeaks_xy.push_back((ii+0.5)*nt_xy/nt);
        rpeaks_xy.push_back((jj+0.5)*2 - nr_xy);
        rmean_xy += rpeaks_xy.back();
        thetamean_xy += thetapeaks_xy.back();
      }
      if(hp_xz.GetBinContent(ii+1, jj+1) == max_xz)
      {
        thetapeaks_xz.push_back((ii+0.5)*nt_xz/nt);
        rpeaks_xz.push_back((jj+0.5)*2 - nr_xz);
        rmean_xz += rpeaks_xz.back();
        thetamean_xz += thetapeaks_xz.back();
      }
      if(hp_yz.GetBinContent(ii+1, jj+1) == max_yz)
      {
        thetapeaks_yz.push_back((ii+0.5)*nt_yz/nt);
        rpeaks_yz.push_back((jj+0.5)*2 - nr_yz);
        rmean_yz += rpeaks_yz.back();
        thetamean_yz += thetapeaks_yz.back();
      }
    }
  }

  // xy PEAK
  rmean_xy = rmean_xy / rpeaks_xy.size();
  thetamean_xy = thetamean_xy / thetapeaks_xy.size();

  // xz PEAK
  rmean_xz = rmean_xz / rpeaks_xz.size();
  thetamean_xz = thetamean_xz / thetapeaks_xz.size();

  // yz PEAK
  rmean_yz = rmean_yz / rpeaks_yz.size();
  thetamean_yz = thetamean_yz / thetapeaks_yz.size();

  rmean_xy = rpeaks_xy[0];
  thetamean_xy = thetapeaks_xy[0];
  rmean_xz = rpeaks_xz[0];
  thetamean_xz = thetapeaks_xz[0];
  rmean_yz = rpeaks_yz[0];
  thetamean_yz = thetapeaks_yz[0];

  //Selection of x,y,z points COMMON to the 3 maxmean+/-1 found in Hough spaces for xy, xz and yz spaces
  for(unsigned int i=0;i<x->size();i++)
  {
    r0_xy = x->at(i)*TMath::Cos(thetamean_xy*PI/180.)+y->at(i)*TMath::Sin(thetamean_xy*PI/180.);
    tmin = thetamean_xy-bint_xy;
    tmax = thetamean_xy+bint_xy;
    if((tmin)<0) tmin = tmin + 180.;
    if((tmax)>180) tmax = tmax - 180.;
    rmin_xy = x->at(i)*TMath::Cos(tmin*PI/180.)+y->at(i)*TMath::Sin(tmin*PI/180.);
    rmax_xy = x->at(i)*TMath::Cos(tmax*PI/180.)+y->at(i)*TMath::Sin(tmax*PI/180.);

    rinf = min( rmean_xy - binr_xy, rmean_xy + binr_xy);
    rsup = max( rmean_xy - binr_xy, rmean_xy + binr_xy);
    if((r0_xy>=rinf || rmin_xy>=rinf || rmax_xy>=rinf) && (r0_xy<=rsup || rmin_xy<=rsup || rmax_xy<=rsup))
    {
      r0_xz = z->at(i)*TMath::Cos(thetamean_xz*PI/180.)+x->at(i)*TMath::Sin(thetamean_xz*PI/180.);
      tmin = thetamean_xz-bint_xz;
      tmax = thetamean_xz+bint_xz;
      if((tmin)<0) tmin = tmin + 180.;
      if((tmax)>180) tmax = tmax - 180.;
      rmin_xz = z->at(i)*TMath::Cos(tmin*PI/180.)+x->at(i)*TMath::Sin(tmin*PI/180.);
      rmax_xz = z->at(i)*TMath::Cos(tmax*PI/180.)+x->at(i)*TMath::Sin(tmax*PI/180.);

      rinf = min( rmean_xz - binr_xz, rmean_xz + binr_xz);
      rsup = max( rmean_xz - binr_xz, rmean_xz + binr_xz);

      if((r0_xz>=rinf || rmin_xz>=rinf || rmax_xz>=rinf) && (r0_xz<=rsup || rmin_xz<=rsup || rmax_xz<=rsup))
      {
        r0_yz = z->at(i)*TMath::Cos(thetamean_yz*PI/180.)+y->at(i)*TMath::Sin(thetamean_yz*PI/180.);
        tmin = thetamean_yz-bint_yz;
        tmax = thetamean_yz+bint_yz;
        if((tmin)<0) tmin = tmin + 180.;
        if((tmax)>180) tmax = tmax - 180.;
        rmin_yz = z->at(i)*TMath::Cos(tmin*PI/180.)+y->at(i)*TMath::Sin(tmin*PI/180.);
        rmax_yz = z->at(i)*TMath::Cos(tmax*PI/180.)+y->at(i)*TMath::Sin(tmax*PI/180.);

        rinf = min( rmean_yz - binr_yz, rmean_yz + binr_yz);
        rsup = max( rmean_yz - binr_yz, rmean_yz + binr_yz);

        if((r0_yz>=rinf || rmin_yz>=rinf || rmax_yz>=rinf) && (r0_yz<=rsup || rmin_yz<=rsup || rmax_yz<=rsup))
        {
          x_out->push_back(x->at(i));
          y_out->push_back(y->at(i));
          z_out->push_back(z->at(i));
          q_out->push_back(q->at(i));
        }
      }
    }
  }
}

// Calculation of the minimal distance between 2 lines in 3D space & calculation of mid-point=>vertex of interaction
void NPL::Tracking::vertex(double *p, double *pp, double &xv,double &yv,double &zv,double &min_dist, double &Theta_tr1, double &Theta_tr2, double &Phi1, double &Phi2, double *VectorTrack1, double *VectorTrack2) {
  double a1 = p[0];
  double a2 = p[2];
  double b1 = p[1];
  double b2 = p[3];
  double ap1 = pp[0];
  double ap2 = pp[2];
  double bp1 = pp[1];
 double bp2 = pp[3];

  // calcul the closest pointis between track1 && track2
  double alpha, beta, A, B, C;

  alpha = (bp1*(a1-ap1)+bp2*(a2-ap2))/(bp1*bp1 + bp2*bp2 + 1);
  beta = (bp1*b1+bp2*b2+1)/(bp1*bp1 + bp2*bp2 + 1);

  A = beta*(bp1*bp1 + bp2*bp2 + 1) - (bp1*b1 + bp2*b2 + 1);
  B = (b1*b1 + b2*b2 + 1) - beta*(bp1*b1+bp2*b2+1);
  C = beta*(bp1*(ap1-a1) + bp2*(ap2-a2)) - (b1*(ap1-a1) + b2*(ap2-a2));

  double sol1, solf1;
  double x,y,z,xp,yp,zp;

  sol1 = -(A*alpha + C)/(A*beta + B);
  solf1 = alpha + beta* sol1;

  // point from track1
  x = a1 + b1*sol1;
  y = a2 + b2*sol1;
  z = sol1;
  // point from track2
  xp = ap1 + bp1*solf1;
  yp = ap2 + bp2*solf1;
  zp = solf1;
  // vertex (mid-ditance between the 2 points) 
  xv = (x+xp)/2.;
  yv = (y+yp)/2.;
  zv = (z+zp)/2.;

  // calulate Theta and Phi 
  double xa,ya,za,zap,xap,yap; 
  double zpoint = 3000;
  xa = a1 + b1*zpoint;
  ya = a2 + b2*zpoint;

  xap = ap1 + bp1*zpoint;
  yap = ap2 + bp2*zpoint;

  //3D unit vectors of tracks 
  
  VectorTrack1[0] = xa-x;
  VectorTrack1[1] = ya-y;
  VectorTrack1[2] = zpoint-z;
  
  VectorTrack2[0] = xap-xp;
  VectorTrack2[1] = yap-yp;
  VectorTrack2[2] = zpoint-zp;
  
  // Aldric Revel version :
  Theta_tr1 = TMath::ATan(TMath::Sqrt((x-a1)*(x-a1)+(y-a2)*(y-a2))/TMath::Abs(sol1));
  Theta_tr2 = TMath::ATan(TMath::Sqrt((xp-ap1)*(xp-ap1)+(yp-ap2)*(yp-ap2))/TMath::Abs(solf1));
  Phi1 = TMath::ATan2(b2,b1);
  Phi2 = TMath::ATan2(bp2,bp1);
  
  Phi1 = 180*Phi1/TMath::Pi();
  Phi2 = 180*Phi2/TMath::Pi();
  Theta_tr1 = Theta_tr1*180./TMath::Pi();
  Theta_tr2 = Theta_tr2*180./TMath::Pi();
  min_dist = sqrt(pow((x-xp),2) + pow((y-yp),2) + pow((z-zp),2));

}

void NPL::Tracking::ParFor_Vertex(double *a, double *b, double *parFit) {
  // input 2 points in 3D, output parFit format for the use in vertex function
  double APX, APY, APZ, X0, Y0, APX0, APY0;

  APX=a[0]-b[0]; APY=a[1]-b[1]; APZ=a[2]-b[2];
  APX0=APX/APZ; APY0=APY/APZ;
  X0=a[0]-a[2]*APX0; Y0=a[1]-a[2]*APY0; //Z0=0

  parFit[0]=X0;
  parFit[1]=APX0;
  parFit[2]=Y0;
  parFit[3]=APY0;
}
