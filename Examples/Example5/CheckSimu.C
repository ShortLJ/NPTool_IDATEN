#include <iostream>
#include <ctime>
#include <cstdlib>
using namespace std;

// ROOT headers
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

// nptool headers
#include "TReactionConditions.h"
#include "TInteractionCoordinates.h"
#include "NPReaction.h"
#include "NPQFS.h"

using namespace NPL;

TCanvas* canvas2;

void CheckSimu(const char * fname = "Example5"){
  // for the style 
  gStyle->SetOptStat(0);
  gROOT->SetStyle("nptool");     
  gROOT->ForceStyle(true);  
  gStyle->SetPalette(1);

  // Open output ROOT file from NPTool simulation run
  TString path = gSystem->Getenv("NPTOOL");
  path += "/Outputs/Simulation/";
  TString inFileName = fname;
  if (!inFileName.Contains("root")) inFileName += ".root";
  TFile *inFile = new TFile(path + inFileName);
  if (!inFile->IsOpen()) exit(1);
  TTree *tree   = (TTree*) inFile->Get("SimulatedTree");

  // Connect the branches of the TTree and activate then if necessary
  // TReactionConditions branch
  TReactionConditions *reacCond = 0;
  tree->SetBranchAddress("ReactionConditions", &reacCond);
  tree->SetBranchStatus("ReactionConditions", 1);
 
  // TInteractionCoordinates branch
  TInteractionCoordinates *interCoord = 0;
  tree->SetBranchAddress("InteractionCoordinates", &interCoord);
  tree->SetBranchStatus("InteractionCoordinates", 1);

  // Declare histograms
  // for emitted particle
  TH1F *hEmThetaCM = new TH1F("hEmThetaCM", "Light Ejectile Theta CM",180, 0, 180);
  TH1F *hEmInternalMomX = new TH1F("hEmInternalMomX", "Internal Momentum (X) of removed cluster",500, -500, 500);
  TH1F *hEmInternalMomY = new TH1F("hEmInternalMomY", "Internal Momentum (Y) of removed cluster",500, -500, 500);
  TH1F *hEmInternalMomZ = new TH1F("hEmInternalMomZ", "Internal Momentum (Z) of removed cluster",500, -500, 500);
  TH1F *hEmTheta1IF = new TH1F("hEmTheta1IF", "Ejectile1 Theta (reaction frame)", 180, 0, 180);
  TH1F *hEmPhi1IF   = new TH1F("hEmPhi1IF",   "Ejectile1 Phi (reaction frame)",   360, 0, 360);
  TH1F *hEmTheta1WF = new TH1F("hEmTheta1WF", "Ejectile1 Theta (world frame)", 180, 0, 180);
  TH1F *hEmPhi1WF   = new TH1F("hEmPhi1WF",   "Ejectile1 Phi (world frame)",   360, 0, 360);
  TH1F *hEmTheta2IF = new TH1F("hEmTheta2IF", "Ejectile2 Theta (reaction frame)", 180, 0, 180);
  TH1F *hEmPhi2IF   = new TH1F("hEmPhi2IF",   "Ejectile2 Phi (reaction frame)",   360, 0, 360);
  TH1F *hEmTheta2WF = new TH1F("hEmTheta2WF", "Ejectile2 Theta (world frame)", 180, 0, 180);
  TH1F *hEmPhi2WF   = new TH1F("hEmPhi2WF",   "Ejectile2 Phi (world frame)",   360, 0, 360);

  TH1F *hEmOpAngle     = new TH1F("hEmOpAngle",  "Opening angle betw. (1) and (2)", 100, 0, 100);
  TH1F *hdPhi     = new TH1F("hdPhi",  "Phi angle difference between (1) and (2)", 200, 0, 200);
  TH2F *hEmE1Theta1  = new TH2F("hEmE1Theta1",  "Kinematics (1)", 900, 0, 90, 1000, 0, 500);
  TH2F *hEmE2Theta2  = new TH2F("hEmE1Theta1",  "Kinematics (2)", 900, 0, 90, 1000, 0, 500);

  TH2F *hEmE1VsE2 = new TH2F("hEmE1VsE2", " E1 VS E2 (reaction frame)", 300, 0, 300,300,0,300);
  TH2F *hEmTheta1VsTheta2 = new TH2F("hEmTheta1VsTheta2", " Theta1 VS Theta2 (reaction frame)", 360, 0, 90,360,0,90);
  TH2F *hEmPhi1VsPhi2 = new TH2F("hEmPhi1VsPhi2", " Phi1 VS Phi2 (reaction frame)", 360, -180, 180,360,-180,180);

  // Read the TTree
  Int_t nentries = tree->GetEntries();
  cout << endl << " TTree contains " << nentries << " events" << endl;

  for (Int_t i = 0; i < nentries; i++) {
    if (i%10000 == 0 && i!=0)  {
      cout.precision(5);
      Double_t percent = (Double_t)i/nentries ;
      cout  << "\r Progression:  " << percent*100 << " %" << flush;
    }
    else if (i==nentries-1)  cout << "\r Progression:" << " 100%" << endl;

    // Get entry
    tree->GetEntry(i);

    // Fill histos
    // ejected particles
    hEmThetaCM  -> Fill(reacCond->GetThetaCM());
    hEmInternalMomX  -> Fill(reacCond->GetInternalMomentum().X());
    hEmInternalMomY  -> Fill(reacCond->GetInternalMomentum().Y());
    hEmInternalMomZ  -> Fill(reacCond->GetInternalMomentum().Z());

    hEmTheta1IF -> Fill(reacCond->GetThetaLab_BeamFrame(0));
    hEmTheta1WF  -> Fill(reacCond->GetThetaLab_WorldFrame(0));
    hEmE1Theta1  -> Fill(reacCond->GetThetaLab_BeamFrame(0), reacCond->GetKineticEnergy(0));
    hEmTheta2IF  -> Fill(reacCond->GetThetaLab_BeamFrame(1));
    hEmTheta2WF  -> Fill(reacCond->GetThetaLab_WorldFrame(1));
    hEmE2Theta2  -> Fill(reacCond->GetThetaLab_BeamFrame(1), reacCond->GetKineticEnergy(1));

    hEmTheta1VsTheta2   -> Fill(reacCond->GetTheta(1), reacCond->GetTheta(0));
    hEmPhi1VsPhi2   -> Fill(reacCond->GetPhi(1), reacCond->GetPhi(0));
    hEmE1VsE2   -> Fill(reacCond->GetKineticEnergy(1), reacCond->GetKineticEnergy(0));

    double theta1 = reacCond->GetThetaLab_BeamFrame(0)*TMath::Pi()/180.;
    double theta2 = reacCond->GetThetaLab_BeamFrame(1)*TMath::Pi()/180.;
    double phi1 = reacCond->GetPhiLab_BeamFrame(0)*TMath::Pi()/180.;
    double phi2 = reacCond->GetPhiLab_BeamFrame(1)*TMath::Pi()/180.;
    double Opang = acos( sin(theta1)*sin(theta2)*cos(phi1-phi2) +
                         cos(theta1)* cos(theta2) );                      
    hEmOpAngle->Fill(Opang*180./TMath::Pi());

    double df = fabs(phi1-phi2)*180./TMath::Pi();

    if(df>0 && df <= 180) hdPhi->Fill(df);
    else hdPhi->Fill(360 - df);
  }


  // Display emmitted paricles histograms
  canvas2 = new TCanvas("canvas2", "Emmited particles properties in reaction frame",1000,1000);
  canvas2->Divide(3,3);

  canvas2->cd(1);
  hEmThetaCM->SetXTitle("#theta_{c.m.}");
  hEmThetaCM->SetYTitle("counts / 1^{#circ}");
  hEmThetaCM->GetYaxis()->SetTitleOffset(1.18);
  hEmThetaCM->GetYaxis()->SetRangeUser(0,400);
  hEmThetaCM->Draw();
  
  canvas2->cd(2);
  hEmE1Theta1->SetXTitle("#theta_{1}");
  hEmE1Theta1->SetYTitle("E_{1} (MeV)");
  hEmE1Theta1->Draw("colz");

  canvas2->cd(3);
  hdPhi->Draw();
  hdPhi->SetXTitle("#phi_{1} - #phi_{2} (deg)");
  hdPhi->SetYTitle("Counts");

  canvas2->cd(4);
  hEmTheta1VsTheta2->Draw("colz");
  hEmTheta1VsTheta2->SetXTitle("#theta_{1} (deg)");
  hEmTheta1VsTheta2->SetYTitle("#theta_{2} (deg)");
  NPL::QFS qfs;
  qfs.ReadConfigurationFile("Example5.reaction");
  qfs.SetMomentumSigma(0.);
  TGraph* Kine = qfs.GetTheta2VsTheta1(1);
  Kine->SetLineWidth(2);
  Kine->SetLineColor(kRed);
  Kine->SetLineStyle(2);
  Kine->Draw("csame");
  TGraph* Kine2 = qfs.GetTheta2VsTheta1(10);
  Kine2->SetLineColor(kRed);
  Kine2->SetMarkerColor(kRed);
  Kine2->SetMarkerStyle(20);
  Kine2->SetMarkerSize(1.3);
  Kine2->Draw("Psame");




  canvas2->cd(5);
  hEmPhi1VsPhi2->Draw("colz");
  hEmPhi1VsPhi2->SetXTitle("#phi_{1} (deg)");
  hEmPhi1VsPhi2->SetYTitle("#phi_{2} (deg)");
  TGraph* KinePhi = qfs.GetPhi2VsPhi1(1);
  KinePhi->SetMarkerColor(kRed);
  KinePhi->SetMarkerSize(0.4);
  KinePhi->Draw("Psame");


  canvas2->cd(6);
  hEmOpAngle->Draw();
  hEmOpAngle->SetXTitle("Opening angle (1-2)  (deg)");
  hEmOpAngle->SetYTitle("Counts");

  canvas2->cd(7);
  hEmInternalMomX->Draw();
  hEmInternalMomX->SetXTitle("P_{x} (MeV/c)");
  hEmInternalMomX->SetYTitle("Counts");

  canvas2->cd(8);
  hEmInternalMomY->Draw();
  hEmInternalMomY->SetXTitle("P_{y} (MeV/c)");
  hEmInternalMomY->SetYTitle("Counts");

  canvas2->cd(9);
  hEmInternalMomZ->Draw();
  hEmInternalMomZ->SetXTitle("P_{z} (MeV/c)");
  hEmInternalMomZ->SetYTitle("Counts");



}


