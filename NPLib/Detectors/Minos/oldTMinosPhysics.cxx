/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Developed by C. Santamaria / CEA Saclay                  *
 *                                                                           *
 * Creation Date  : 2014/11/24                                               *
 * Last update    : 2019/09 implemeted in NPTool by Cyril Lenain             *
 *                  lenain@lpccaen.in2p3.fr                                  *     
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Minos Treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TMinosPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "TMinosClust.h"
#include "TMinosResult.h"

//   ROOT
#include "TChain.h"
#include "TROOT.h"
#include "TMath.h"
/* #include "TMinuit2.h" */
#include "Math/Vector3D.h"
#include "TStyle.h"
#include "TClonesArray.h" 
#include "Math/Minimizer.h"
#include "Math/Factory.h"

TMinosPhysics* current_phy = 0;

///////////////////////////////////////////////////////////////////////////
TMinosPhysics::TMinosPhysics()
  : m_EventData(new TMinosData),
  m_PreTreatedData(new TMinosData),
  m_EventPhysics(this),
  Tracking_functions(new Tracking),
  m_E_RAW_Threshold(0), // adc channels
  m_E_Threshold(0),     // MeV
  hfit(new TH1F("hfit","hfit",512,0,512)),
  grxztmp(new TGraph()),
  gryztmp(new TGraph()),
  npoint_temp(0),
  allevt_2pfiltered(0),
  m_NumberOfDetectors(0) {
    current_phy = this;
    data_result.SetClass("TMinosResult");
    fitdata.SetClass("TMinosClust");
  }

int NclusterFit;

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TMinosPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::AddDetector(double R, double Theta, TVector3 POS, double Phi, double TargetLenght){
  // Compute the TVector3 corresponding
  // Call the cartesian method
} 

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();
}

// definition of the fit function for the signal Q(t)
double TMinosPhysics::conv_fit(double *x, double *p){
  static double p0, p1, p2, x0;
  p0=p[0]; p1=p[1]; p2=p[2];
  x0 = x[0];
  if(x0 > p1 && x0<512.) 
    return (p0 * exp(-3.*(x0-p1)/p2)  * sin((x0-p1)/p2) * ((x0-p1)/p2)*((x0-p1)/p2)*((x0-p1)/p2) + 250);
  else return (250);
}

///////////////////////////////////////////////////////////////////////////
void TMinosPhysics::PreTreat() {

  // This method applies thresholds and calibrations 
  // Isolate 2D tracks with the modified HoughTransformation
  // Calculate z for each pad 

  data_result.Clear();
  fitdata.Clear();

  for(int i = 0; i <2 ; i++){
    parFit1.push_back(-1000);
    parFit2.push_back(-1000);
    parFit3.push_back(-1000);
    parFit4.push_back(-1000);
  }

  Xvertex = -1000;
  Yvertex = -1000;
  Zvertex = -1000;
  Theta_12= -1000;

  Theta_1= -1000;
  Theta_2= -1000;

  Phi1 = -1000;
  Phi2 = -1000;

  Dmin = -1000;

  TMinosClust *minosfitdata;
  ////////////////////fit variables//////////////////////  

  /* fit_function = new TF1("fit_function",conv_fit, 0, 511, 3); */
  fit_function->Clear(); 
  int fit2DStatus = 0.;
  double fit_function_max = 0., fit_function_Tpad = 0.;
  double Chi2 = 0.;
  double Tshaping = 333.9;
  const double TimeBinElec = 30.; // in ns
  if(SimulationBool){ // Simulated data
    UperFitLim = 1000000000.; //   
    MINOSthresh = 899.;
    VDrift = 0.03475;  
    Z_Shift = 0; 
  }
  else{ // Experiment data
    MINOSthresh = 100;
    UperFitLim = 100000.; 
    coef_Vdrift = 1.125; 
    ZRot_Minos = 35; // Rotation of Minos along Z axis s034
    //to fixe VDrift
    Z_Shift = -18.95 + 4;
  }
  
  ///////////////////////////////////////////////////////

  ClearPreTreatedData();

  /// Read Q(t) signal to fit and extract DriftTime

  unsigned int NbOfPad = m_EventData->GetPadNumberMult();
  filled =0;trackNbr=0;
  trackNbr_FINAL=0;
  ringsum=0;
  zmax=0.;Iteration=0;
  for (unsigned int e = 0; e < NbOfPad; e++){

    x_mm = 0.; y_mm = 0.; z_mm = 0.;maxCharge = 0.; t_pad =0; q_pad=0;
    x_mm = m_EventData->GetPadX(e);
    y_mm = m_EventData->GetPadY(e);

    if ( !(abs(x_mm)<0.0001 && abs(y_mm)<0.0001) ) { 
      for (UShort_t o = 0; o < m_EventData->GetTime(e).size() ; ++o){
        if(m_EventData->GetCharge(e)[o]> maxCharge){
          maxCharge = m_EventData->GetCharge(e)[o];
        }
      }
      if(maxCharge >= MINOSthresh ) { // to remove rings 1,2,3,4 and 16,17,18      
        /* if((maxCharge >= MINOSthresh && round((sqrt(x_mm*x_mm+y_mm*y_mm)-44.15)/2.1) < 16 && round((sqrt(x_mm*x_mm+y_mm*y_mm)-44.15)/2.1) >4 )) { // to remove rings 1,2,3,4 and 16,17,18 */      
        Xpad.push_back(x_mm);
        Ypad.push_back(y_mm);
        Qpad.push_back(maxCharge);
        filled++;  // 
      }
    }
  } //end of NbOfPad 

  if(filled>0){ // Nbr of activated pads
    while(Xpad.size()>=7 && Iteration< 20 ) {   //  !!!!!!!!  
      /* while(Xpad.size()>=10 && Iteration< 20 ) {   //  !!!!!!!! */  
      filter_result = 0; // Nbr of pads in the Track
      Iteration++;
      
      static vector<double> XpadTemp, YpadTemp, QpadTemp;
      static vector <int> clusterringboolTemp;
      XpadTemp.clear();                
      YpadTemp.clear();                
      QpadTemp.clear();                
      clusterringboolTemp.clear();      

      filter_result = Tracking_functions->Hough_modified(&Xpad, &Ypad, &Qpad, &XpadTemp, &YpadTemp, &QpadTemp, &clusterringboolTemp);

      if(filter_result<0) break;
      if(filter_result>7 ){ 
        trackNbr++;
        for(int ik=0; ik<filter_result; ik++) { // selected pads
          XpadNew.push_back(XpadTemp[ik]);                
          YpadNew.push_back(YpadTemp[ik]);                
          ZpadNew.push_back(-10000);                 
          QpadNew.push_back(QpadTemp[ik]);                 
          clusterringbool.push_back(clusterringboolTemp[ik]);                
          clusternbr.push_back(Iteration);
          clusterpads.push_back(filter_result);
        }
      }
    }
  }
  for(unsigned int il=0; il<XpadNew.size(); il++){
    minosfitdata = (TMinosClust*)fitdata.ConstructedAt(il);
    minosfitdata->Set(XpadNew[il], YpadNew[il],
        -10000, -10000, QpadNew[il], 
        clusternbr[il],clusterpads[il], 0.);
    /* ZpadNew.push_back(-10000.); */
  }  

  if(trackNbr == 2 || trackNbr == 1){  ///////////

    /* if(filled==0) cerr << "Error !!!" << endl; */
    //-------------------------------------------------------
    //  STEP 4.1: Fitting taken pads for Qmax and Ttrig info
    //-------------------------------------------------------
    for(unsigned int i=0 ; i< NbOfPad; i++ ) {
      hfit->Reset();
      bool fitbool = false;
      x_mm = m_EventData->GetPadX(i);
      y_mm = m_EventData->GetPadY(i);

      for(unsigned int jj=0; jj<XpadNew.size(); jj++) {
        if( abs(XpadNew[jj]-x_mm)<0.0001 && abs(YpadNew[jj]-y_mm)<0.0001) {
          fitbool = true;
          indexfill=jj;
          break;
        }
      }
      if( fitbool==true ) {

        /* TGraph gfit(  m_EventData->GetTime(i).size() ,&(m_EventData->GetTime(i)[0]),&(m_EventData->GetCharge(i)[0])); */
        /* for (UShort_t o = 0; o < m_EventData->GetTime(i).size() ; ++o){ */
        /*   gfit.SetPoint(o,m_EventData->GetTime(i)[o],m_EventData->GetCharge(i)[o]+250); */
        /*   if(m_EventData->GetCharge(i)[o]+250> hfit_max){ */
        /*     hfit_max_T = m_EventData->GetTime(i)[o]; */
        /*     hfit_max=m_EventData->GetCharge(i)[o]+250; */
        /*   } */
        /* } */  
        
        for(Int_t j=0; j< m_EventData->GetTime(i).size(); j++) {
          if(m_EventData->GetCharge(i)[j]>=0){
            hfit->SetBinContent(hfit->FindBin(m_EventData->GetTime(i)[j]), m_EventData->GetCharge(i)[j]+250);
          }
        }

        // Fitting the hfit histogram of last channel if not empty
        if(hfit->GetSumOfWeights()>3000) {
          hfit->GetXaxis()->SetRange(0,510);
          hfit_max = hfit->GetMaximum();
          hfit_max_T = hfit->GetMaximumBin();
          T_min=-1;
          T_max=-1;

          // Find the T_min & T_max limits of the signal non zero
          for(int h=hfit_max_T;h>0;h--) {
            if(T_min == -1 && (hfit->GetBinContent(h))<=250 ) {
              T_min = h;
              break;
            }
          }
          for(int h=hfit_max_T;h<510;h++) {
            if(T_max == -1 && (hfit->GetBinContent(h)) < 1  ) {
              T_max = h;
              break;
            }
          }

          if((hfit_max_T-3.5*(Tshaping/TimeBinElec)) > T_min) T_min = hfit_max_T-2*Tshaping/TimeBinElec;
          if((hfit_max_T+10) < T_max || T_max==-1) T_max = hfit_max_T+10.;
          T_min = max(T_min,0.);
          if(T_max>510) T_max = 510;

          // Set fit parameters
          fit_function->SetParameter(0, hfit_max-250.);
          fit_function->SetParameter(1,hfit_max_T - Tshaping/TimeBinElec);
          fit_function->SetParameter(2, Tshaping/TimeBinElec);
          fit_function->SetParLimits(0,0,UperFitLim);
          fit_function->SetParLimits(1,-20,512);  
          fit_function->SetParLimits(2,0,512);

          fit2DStatus = hfit->Fit(fit_function,"QN","",T_min,T_max);
          /* fit2DStatus = gfit.Fit(fit_function,"QN","",hfit_max-2*Tshaping,hfit_max+2*Tshaping); */

          double fit_function_max = 0., fit_function_Tpad = 0.;
          if(fit2DStatus==0) {
            Chi2 = fit_function->GetChisquare();
            fit_function_max = fit_function->GetMaximum();
            fit_function_Tpad = fit_function->GetParameter(1);
          }
          //attribute q_pad and z_mm value
          if(fit2DStatus!=0 || fit_function_max<=20. || 
              fit_function_max >= UperFitLim || fit_function_Tpad<=0.15 || 
              fit_function_Tpad>=513. || fit_function->GetParameter(2)<=0.15 || 
              fit_function->GetParameter(2)>=513.) {
            q_pad = hfit_max-250.; 
            z_mm = -10000;
          }
          else {
            t_pad = fit_function_Tpad;
            if(SimulationBool)z_mm = t_pad*TimeBinElec*VDrift;  // for simu
            else{
              int ring = round((sqrt(x_mm*x_mm + y_mm*y_mm)-44.15)/2.1);
              z_mm =(t_pad*TimeBinElec+DelayTrig[ring])*(VdriftperRing[ring]*coef_Vdrift); 
            }
            Q_Pad.push_back(fit_function_max-250);
            T_Pad.push_back(t_pad*TimeBinElec);
            X_Pad.push_back(x_mm);
            Y_Pad.push_back(y_mm);
            Z_Pad.push_back(z_mm);
            q_pad = fit_function_max-250.;

          }
          ZpadNew[indexfill] = z_mm + Z_Shift; 
          QpadNew[indexfill] = q_pad;
          minosfitdata = (TMinosClust*)fitdata.ConstructedAt(indexfill);
          minosfitdata->Set(XpadNew[indexfill], YpadNew[indexfill], t_pad*TimeBinElec, z_mm, q_pad, clusternbr[indexfill], clusterpads[indexfill], Chi2);
        } // end if SumOfWeights > 0
      } // end if fitbool==True
      else continue;
    } // end of NbOfPad

    //-------------------------------------------------------
    //  STEP 3.2:  Filtering the tracks off possible noise 
    //             with Hough3D (3*2D planes)
    //-------------------------------------------------------

    cluster_temp = 0;
    int ringtouch[19]={0};      
    static vector<double> xin, yin, zin, qin;
    static vector<double> xout, yout, zout, qout;
    static vector<double> xoutprime, youtprime, zoutprime, qoutprime;

    for(unsigned int i=0;i<(XpadNew.size());i++){
      if(xin.size()>0 && ((cluster_temp != int(clusternbr[i]) && i!=0) || i==(XpadNew.size()-1))){ // We fill xin until the next cluster 
        Tracking_functions->Hough_3D(&xin, &yin, &zin, &qin, &xout, &yout, &zout, &qout);
        for(unsigned int ij=0; ij<xout.size();ij++){
          if(zout[ij]>zmax) {zmax = zout[ij];}
          ringtouch[int(round((sqrt(xout[ij]*xout[ij]+yout[ij]*yout[ij])-44.15)/2.1))]++; // Corr by Cyril
        }
        for(int ko=0; ko<19; ko++){
          if(ringtouch[ko]>0) ringsum++;
        }
        // Tracks of interest: >10 pads and >=12 rings hit
        if(xout.size()>10 && ringsum>=8){
          npoint=0;
          trackNbr_FINAL++;
          double charge_temp=0.;
          double lenght_temp=0;
          double zintarget=300;
          double zouttarget=0;
          int indexin=0; int indexout=0;
          for(unsigned int ij=0; ij<xout.size(); ij++){
            point.SetXYZ(xout[ij],yout[ij],zout[ij]);
            if(!SimulationBool)point.RotateZ(ZRot_Minos*TMath::DegToRad());//15.6*TMath::DegToRad() CHANGE####
            xoutprime.push_back(point.X());youtprime.push_back(point.Y());zoutprime.push_back(point.Z());
            TOTxoutprime.push_back(point.X());
            // Added line to see events used for fit of lines in tree
            TOTyoutprime.push_back(point.Y());
            TOTzoutprime.push_back(point.Z());
            TOTqout.push_back(qout[ij]);
            // save cluster_temp
            // save numbtrack
            trackclusternbr.push_back(clusternbr[i]);
            tracknbr.push_back( trackNbr_FINAL ) ;

            grxztmp->SetPoint(npoint,zoutprime[ij],xoutprime[ij]);
            gryztmp->SetPoint(npoint,zoutprime[ij],youtprime[ij]);
            charge_temp += qout[ij];
            minosdata_result = (TMinosResult*)data_result.ConstructedAt(array_final);
            minosdata_result->Set(xoutprime[ij], youtprime[ij], zoutprime[ij], qout[ij], trackNbr_FINAL, xout.size(), zmax);
            array_final++;
            npoint++;

            if(zoutprime[ij]<zintarget) {zintarget=zoutprime[ij];indexin=ij;}
            if(zoutprime[ij]>zouttarget) {zouttarget=zoutprime[ij];indexout=ij;}
          }

          if(xout.size()>0)lenght_temp=sqrt((zoutprime[indexin]-zoutprime[indexout])*(zoutprime[indexin]-zoutprime[indexout]) +(youtprime[indexin]-youtprime[indexout])*(youtprime[indexin]-youtprime[indexout]) +(xoutprime[indexin]-xoutprime[indexout])*(xoutprime[indexin]-xoutprime[indexout]));

          grxz.push_back(*grxztmp);
          gryz.push_back(*gryztmp);
          grxztmp->Set(0);
          /* if(maxCharge >= MINOSthresh ) { */      
          gryztmp->Set(0);

          chargeTot.push_back(charge_temp);
          lenght.push_back(lenght_temp);                        
        } // end of if(xout.size()>10 &&)
        
        xin.clear(); yin.clear(); zin.clear(); qin.clear(); 
        xout.clear(); yout.clear(); zout.clear(); qout.clear();
        xoutprime.clear(); youtprime.clear(); zoutprime.clear();
        
        npoint_temp=0; ringsum=0; zmax=0.;
        for(int ko=0; ko<18; ko++) 
          ringtouch[ko] = 0;
        
        } // if(xin.size()>0 && (( cluster_temp .... )

        cluster_temp = clusternbr[i]; // Number of the track/cluster
        if(!(clusterpads[i]>=10 && clusterringbool[i]==1 && ZpadNew[i]>-100 && ZpadNew[i]<=310)) continue;
        else
        {
          xin.push_back(XpadNew[i]);
          yin.push_back(YpadNew[i]);
          zin.push_back(ZpadNew[i]);
          qin.push_back(QpadNew[i]);
          npoint_temp++;
        }

      }//end of PadNews

      /* //------------------------------------------------------- */
      /* //  STEP 3.2:  Fitting the filtered tracks in 3D */ 
      /* //             (weight by charge, TMinuit) */
      /* //------------------------------------------------------- */


      if(trackNbr_FINAL== 2 || trackNbr_FINAL == 1){

        //////////Minimization in 2D to reconstruct track lines
        allevt_2pfiltered++;
        for(int itr= 0 ; itr < trackNbr_FINAL; itr++) {
          pStart[0]=0; pStart[2]=0; pStart[1]=1; pStart[3]=3;

          min = new TMinuit(4);
          min->SetPrintLevel(-1);
          arglist[0] = 3;

          Tracking_functions->FindStart(pStart,chi,fitStatus, &grxz.at(itr), &gryz.at(itr));

          NclusterFit = itr+1;
          current_phy=this;

          min->SetFCN(SumDistance);

          // Set starting values and step sizes for parameters
          min->mnparm(0,"x0",pStart[0],0.1,-500,500,iflag);
          min->mnparm(1,"Ax",pStart[1],0.1,-10,10,iflag);
          min->mnparm(2,"y0",pStart[2],0.1,-500,500,iflag);
          min->mnparm(3,"Ay",pStart[3],0.1,-10,10,iflag);

          arglist[0] = 200; // number of function calls
          arglist[1] = 0.000001; // tolerance

          min->mnexcm("MIGRAD",arglist,2,iflag); // minimization with MIGRAD
          min->mnstat(amin,edm,errdef,nvpar,nparx,iflag);  //returns current status of the minimization

          // get fit parameters
          for(int i = 0; i <4; i++) min->GetParameter(i,parFit_temp[i],err_temp[i]);

          /* if( (parFit_temp[0] >-499 && parFit_temp[0]<499) && (parFit_temp[2] >-499 && parFit_temp[2]<499)) { */
          /* parFit1[itr] = push_back(parFit_temp[0]); */
          /* parFit2.push_back(parFit_temp[1]); */
          /* parFit3.push_back(parFit_temp[2]); */
          /* parFit4.push_back(parFit_temp[3]); */
          /* } */
          parFit1[itr]=parFit_temp[0];
          parFit2[itr]=parFit_temp[1];
          parFit3[itr]=parFit_temp[2];
          parFit4[itr]=parFit_temp[3];

          delete min;
        }

        static double ParTrack1[4], ParTrack2[4];
        static double VectorTrack11[3], VectorTrack22[3];

        ParTrack1[0] = parFit1[0];
        ParTrack1[1] = parFit2[0];
        ParTrack1[2] = parFit3[0];
        ParTrack1[3] = parFit4[0];

        if(trackNbr_FINAL==2){
          ParTrack2[0] = parFit1[1];
          ParTrack2[1] = parFit2[1];
          ParTrack2[2] = parFit3[1];
          ParTrack2[3] = parFit4[1];
        }
        else if(trackNbr_FINAL==1){// The vertex is still calculated with Zaxis
          ParTrack2[0] = 0;
          ParTrack2[1] = 0;
          ParTrack2[2] = 0;
          ParTrack2[3] = 0;
        }

        Dmin=-100, Theta_1 = -1, Theta_2 = -1;

        Tracking_functions->vertex(ParTrack1, ParTrack2, Xvertex, Yvertex, Zvertex, Dmin, Theta_1, Theta_2, Phi1, Phi2, VectorTrack11, VectorTrack22);

        VectorTrack1.SetXYZ(VectorTrack11[0],VectorTrack11[1],VectorTrack11[2]);
        VectorTrack2.SetXYZ(VectorTrack22[0],VectorTrack22[1],VectorTrack22[2]);
        VectorTrack1 = VectorTrack1.Unit();
        VectorTrack2 = VectorTrack2.Unit();
        Theta_12 = VectorTrack1.Angle(VectorTrack2)*180/TMath::Pi();
      }// end if trackNbr_FINAL>=1

    } // end loop if 0<trackNbr < 5
    // instantiate CalibrationManager
    static CalibrationManager* Cal = CalibrationManager::getInstance();
  }

  ///////////////////////////////////////////////////////////////////////////
  void TMinosPhysics::ReadAnalysisConfig() {
    bool ReadingStatus = false;
    // path to file
    string FileName = "./configs/ConfigMinos.dat";

    // open analysis config file
    ifstream AnalysisConfigFile;
    AnalysisConfigFile.open(FileName.c_str());

    if (!AnalysisConfigFile.is_open()) {
      cout << " No ConfigMinos.dat found: Default parameter loaded for Analayis " << FileName << endl;
      return;
    }
    cout << " Loading user parameter for Analysis from ConfigMinos.dat " << endl;

    // Save it in a TAsciiFile
    TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
    asciiConfig->AppendLine("%%% ConfigMinos.dat %%%");
    asciiConfig->Append(FileName.c_str());
    asciiConfig->AppendLine("");
    // read analysis config file
    string LineBuffer,DataBuffer,whatToDo;
    while (!AnalysisConfigFile.eof()) {
      // Pick-up next line
      getline(AnalysisConfigFile, LineBuffer);

      // search for "header"
      string name = "ConfigMinos";
      if (LineBuffer.compare(0, name.length(), name) == 0) 
        ReadingStatus = true;

      // loop on tokens and data
      while (ReadingStatus ) {
        whatToDo="";
        AnalysisConfigFile >> whatToDo;

        // Search for comment symbol (%)
        if (whatToDo.compare(0, 1, "%") == 0) {
          AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
        }

        else if (whatToDo=="E_RAW_THRESHOLD") {
          AnalysisConfigFile >> DataBuffer;
          m_E_RAW_Threshold = atof(DataBuffer.c_str());
          cout << whatToDo << " " << m_E_RAW_Threshold << endl;
        }

        else if (whatToDo=="E_THRESHOLD") {
          AnalysisConfigFile >> DataBuffer;
          m_E_Threshold = atof(DataBuffer.c_str());
          cout << whatToDo << " " << m_E_Threshold << endl;
        }

        else {
          ReadingStatus = false;
        }
      }
    }

  }

  double TMinosPhysics::distance2(double x,double y,double z, double *p) {
    // distance line point is D= | (xp-x0) cross  ux |
    // where ux is direction of line and x0 is a point in the line (like t = 0)
    ROOT::Math::XYZVector xp(x,y,z); //point of the track
    ROOT::Math:: XYZVector x0(p[0], p[2], 0. );
    ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); //line
    ROOT::Math::XYZVector u = (x1-x0).Unit();
    double d2 = ((xp-x0).Cross(u)) .Mag2();
    return d2;
  }

  void TMinosPhysics::SumDistance(int &, double *, double & sum, double * par,  int) {
    TMinosPhysics* phy = current_phy;
    int nused=0;
    double qtot=0;
    sum = 0;

    for(int i=0; i < phy->data_result.GetEntriesFast(); i++)
    {
      phy->minosdata_result = (TMinosResult*)phy->data_result.At(i);
      if(phy->minosdata_result->n_Cluster==NclusterFit)
      {
        float x=phy->minosdata_result->x_mm;
        float y=phy->minosdata_result->y_mm;
        float z=phy->minosdata_result->z_mm;
        float q=phy->minosdata_result->Chargemax;
        //if(nused<2)cout<<minosdata_result->n_Cluster<<" "<<x<<" "<<y<<" "<<z<<" "<<q<<endl;
        double d = TMinosPhysics::distance2(x, y, z, par);
        sum += d*q;       
        nused++;
        qtot+=q;
      }
    }
    //sum/=nused;
    sum/=qtot;
    return;
  }

  ///////////////////////////////////////////////////////////////////////////
  void TMinosPhysics::Clear() {

    Xpad.clear();
    Ypad.clear();
    Qpad.clear();
    XpadNew.clear();
    YpadNew.clear();
    ZpadNew.clear();
    QpadNew.clear();

    clusterringboolTemp.clear();
    clusterringbool.clear();
    clusternbr.clear();
    clusterpads.clear();  
    hfit->Reset();

    /* xoutprime.clear(); */
    /* youtprime.clear(); */
    /* zoutprime.clear(); */

    trackclusternbr.clear();
    tracknbr.clear();
    TOTxoutprime.clear();
    TOTyoutprime.clear();
    TOTzoutprime.clear();
    TOTqout.clear(); 

    lenght.clear();
    chargeTot.clear();
    parFit1.clear();
    parFit2.clear();
    parFit3.clear();
    parFit4.clear();

    grxz.clear();
    gryz.clear();

    hfit->Reset();

    filled=0;
    indexfill = 0;
    ChargeBin = 0.;
    maxCharge = 0.;
    Iteration=0;
    filter_result=0;
    fit2DStatus=0;
    trackNbr=0;
    trackNbr_FINAL=0;
    x_mm = 0.; y_mm = 0.; z_mm = 0.; q_pad = 0.; t_pad = 0.;
    array_final=0;
    ringsum=0;
    zmax=0.;

  }

  ///////////////////////////////////////////////////////////////////////////
  void TMinosPhysics::ReadConfiguration(NPL::InputParser parser) {
  }


  ///////////////////////////////////////////////////////////////////////////
  void TMinosPhysics::AddParameterToCalibrationManager() {
    CalibrationManager* Cal = CalibrationManager::getInstance();
    for (int i = 0; i < m_NumberOfDetectors; ++i) {
      Cal->AddParameter("Minos", "D"+ NPL::itoa(i+1)+"_ENERGY","Minos_D"+ NPL::itoa(i+1)+"_ENERGY");
      Cal->AddParameter("Minos", "D"+ NPL::itoa(i+1)+"_TIME","Minos_D"+ NPL::itoa(i+1)+"_TIME");
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  void TMinosPhysics::InitializeRootInputRaw() {

    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("Minos",  true );
    inputChain->SetBranchAddress("Minos", &m_EventData );

    fit_function = new TF1("fit_function",conv_fit, 0, 511, 3);

    if(NPOptionManager::getInstance()->HasDefinition("simulation")){
      cout << "Considering input data as simulation"<< endl;
      SimulationBool = true;
    }
    else{
      cout << "Considering input data as real" << endl;

      SimulationBool = false;

      ifstream calibFile2("Vdrift.txt");
      string buffer2;
      getline(calibFile2, buffer2);
      double vdriftR;
      int i = 0;
      while(calibFile2 >> vdriftR){
        VdriftperRing[i] = vdriftR; // ns, s034 par.
        i++;
      }

      ifstream calibFile("Time_Offset.txt");
      string buffer;
      getline(calibFile, buffer);
      double offset;
      i = 0;
      while(calibFile >> offset){
        DelayTrig[i] = offset; // ns, s034 par.
        i++;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  void TMinosPhysics::InitializeRootInputPhysics() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchAddress("Minos", &m_EventPhysics);
  }

  ///////////////////////////////////////////////////////////////////////////
  void TMinosPhysics::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    outputTree->Branch("Minos", "TMinosPhysics", &m_EventPhysics);
  }

  ////////////////////////////////////////////////////////////////////////////////
  //            Construct Method to be pass to the DetectorFactory              //
  ////////////////////////////////////////////////////////////////////////////////
  NPL::VDetector* TMinosPhysics::Construct() {
    return (NPL::VDetector*) new TMinosPhysics();
  }

  ////////////////////////////////////////////////////////////////////////////////
  //            Registering the construct method to the factory                 //
  ////////////////////////////////////////////////////////////////////////////////
  extern "C"{
    class proxy_Minos{
      public:
        proxy_Minos(){
          NPL::DetectorFactory::getInstance()->AddToken("Minos","Minos");
          NPL::DetectorFactory::getInstance()->AddDetector("Minos",TMinosPhysics::Construct);
        }
    };

  proxy_Minos p_Minos;
  }

