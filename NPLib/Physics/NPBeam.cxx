/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien MATTA contact address: matta@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date   : January 2013                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with Beam:                                               *
 * User can enter various parameter, such as emittance or use ASCII or root  *
 * TH1F distribution                                                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <vector>

// NPL header
#include "NPBeam.h"
#include "NPFunction.h"
#include "NPOptionManager.h"

// ROOT Header
#include "TDirectory.h"

// Use CLHEP System of unit and Physical Constant
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"

#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

using namespace NPL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
Beam::Beam(){
  fExcitationEnergy = 0; 
  fEnergy = 0;
  fSigmaEnergy = -1 ;
  fMeanX = 0 ;
  fMeanY = 0 ;
  fSigmaX = -1;
  fSigmaY = -0;
  fMeanThetaX = 0 ;
  fMeanPhiY = 0 ;
  fSigmaThetaX = -1 ;
  fSigmaPhiY = -1 ;
  fTargetSize = 0 ;
  fEffectiveTargetSize = 0 ;
  fTargetThickness = 0 ;
  fEffectiveTargetThickness = 0 ;
  fTargetAngle = 0 ;
  fTargetZ = 0 ;
  fZEmission=-1*NPUNITS::m;
  fZProfile=0;
  fVerboseLevel = NPOptionManager::getInstance()->GetVerboseLevel();

  // case of user given distribution
  // do that to avoid warning from multiple Hist with same name...
  int offset = 0;
  while(gDirectory->FindObjectAny(Form("EnergyHist_%i",offset))!=0)
    ++offset;

  fEnergyHist  = new TH1F(Form("EnergyHist_%i",offset),"EnergyHist",1,0,1);
  fXThetaXHist = new TH2F(Form("XThetaXHis_%i",offset),"XThetaXHis",1,0,1,1,0,1);
  fYPhiYHist   = new TH2F(Form("YPhiYHist_%i",offset),"YPhiYHist",1,0,1,1,0,1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Beam::Beam(string name){
  SetUp(name);
  fEnergy = 0;
  fExcitationEnergy = 0; 
  fSigmaEnergy = -1 ;
  fMeanX = 0 ;
  fMeanY = 0 ;
  fSigmaX = -1;
  fSigmaY = -0;
  fMeanThetaX = 0 ;
  fMeanPhiY = 0 ;
  fSigmaThetaX = -1 ;
  fSigmaPhiY = -1 ;
  fTargetSize = 0 ;
  fEffectiveTargetSize = 0 ;
  fTargetThickness = 0 ;
  fEffectiveTargetThickness = 0 ;
  fTargetAngle = 0 ;
  fTargetZ = 0 ;
  fZEmission=-1*NPUNITS::m;
  fZProfile=0;
  fVerboseLevel = NPOptionManager::getInstance()->GetVerboseLevel();

  // case of user given distribution
  // do that to avoid warning from multiple Hist with same name...
  int offset = 0;
  while(gDirectory->FindObjectAny(Form("EnergyHist_%i",offset))!=0)
    ++offset;

  fEnergyHist  = new TH1F(Form("EnergyHist_%i",offset),"EnergyHist",1,0,1);
  fXThetaXHist = new TH2F(Form("XThetaXHis_%i",offset),"XThetaXHis",1,0,1,1,0,1);
  fYPhiYHist   = new TH2F(Form("YPhiYHist_%i",offset),"YPhiYHist",1,0,1,1,0,1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Beam::~Beam(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Beam::ReadConfigurationFile(string Path){
  NPL::InputParser parser(Path);
  ReadConfigurationFile(parser);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Beam::ReadConfigurationFile(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Beam");
  if(blocks.size()>0 && NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// Beam found " << endl; 

  vector<string> token   = {"Particle"};
  vector<string> energyA = {"Energy","SigmaEnergy"};
  vector<string> energyP = {"EnergyProfilePath"};
  vector<string> emmitA  = {"SigmaThetaX","SigmaPhiY","SigmaX","SigmaY","MeanThetaX","MeanPhiY","MeanX","MeanY"};
  vector<string> emmitP  = {"XThetaXProfilePath","YPhiYProfilePath"};


  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      SetUp(blocks[i]->GetString("Particle"));

      if(blocks[i]->HasToken("ExcitationEnergy"))
        fExcitationEnergy = blocks[i]->GetDouble("ExcitationEnergy","MeV");

      if(blocks[i]->HasToken("ZEmission"))
        fZEmission = blocks[i]->GetDouble("ZEmission","mm");

      if(blocks[i]->HasToken("ZProfile"))
        fZProfile = blocks[i]->GetDouble("ZProfile","mm");

      // Energy analytic
      if(blocks[i]->HasTokenList(energyA)){
        fEnergy = blocks[i]->GetDouble("Energy","MeV");
        fSigmaEnergy= blocks[i]->GetDouble("SigmaEnergy","MeV");
      }

      // Energy profile
      else if(blocks[i]->HasTokenList(energyP)){
        vector<string> FileName = blocks[i]->GetVectorString("EnergyProfilePath");
        SetEnergyHist( Read1DProfile(FileName[0], FileName[1]));
      }

      else{
        cout << "ERROR: check your input file formatting " << endl; 
        exit(1);
      }

      // Emmitance analytic
      if(blocks[i]->HasTokenList(emmitA)){
        fSigmaThetaX = blocks[i]->GetDouble("SigmaThetaX","deg");
        fSigmaPhiY = blocks[i]->GetDouble("SigmaPhiY","deg");
        fSigmaX = blocks[i]->GetDouble("SigmaX","mm");
        fSigmaY = blocks[i]->GetDouble("SigmaY","mm");
        fMeanThetaX = blocks[i]->GetDouble("MeanThetaX","deg");
        fMeanPhiY = blocks[i]->GetDouble("MeanPhiY","deg");
        fMeanX = blocks[i]->GetDouble("MeanX","mm");
        fMeanY = blocks[i]->GetDouble("MeanY","mm");
      }
      // Emmitance profile
      else if(blocks[i]->HasTokenList(emmitP)){
        vector<string> XThetaX= blocks[i]->GetVectorString("XThetaXProfilePath");
        SetXThetaXHist( Read2DProfile(XThetaX[0], XThetaX[1]));
        vector<string> YPhiY= blocks[i]->GetVectorString("YPhiYProfilePath");
        SetYPhiYHist( Read2DProfile(YPhiY[0], YPhiY[1]));
        
      }

      else{
        cout << "ERROR: check your input file formatting \033[0m" << endl; 
        exit(1);
      }
    }

    else{
      cout << "ERROR: check your input file formatting \033[0m" << endl;
      exit(1);
    }
  }
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "\033[0m" ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Beam::GenerateRandomEvent(double& E, double& X, double& Y, double& Z, double& ThetaX, double& PhiY ){
  double X0,Y0;
  X0 = Y0 = 1*km;

  // ENERGY //
  // Gaussian energy distribution 
  if(fSigmaEnergy!=-1){
    E=-1;
    while(E<0)
      E = gRandom->Gaus(fEnergy,fSigmaEnergy);
  }
  // User Profile
  else
    E = fEnergyHist->GetRandom();
  
  // POSITION/DIRECTION AT Z PROFILE//
  // Gaussian Distribution
  if(fSigmaX!=-1){
      NPL::RandomGaussian2D(fMeanX, fMeanThetaX, fSigmaX, fSigmaThetaX, X0, ThetaX);
      NPL::RandomGaussian2D(fMeanY, fMeanPhiY, fSigmaY, fSigmaPhiY, Y0, PhiY);
  }

  // Profile
  else{
      fXThetaXHist->GetRandom2(X,ThetaX);
      fYPhiYHist->GetRandom2(Y,PhiY);
  }
  // Direction
  double Xdir = sin(ThetaX); // cos(90-x) = sin(x)
  double Ydir = sin(PhiY); 
  double Zdir = sqrt(1-Xdir*Xdir-Ydir*Ydir); // alpha^2 + beta^2 + gamma^2 = 1
  // Stretch factor to extend unitary vector from ZEmission to ZProfile
  double S = fZProfile-fZEmission ;
  TVector3 BeamDir(Xdir*S/Zdir,Ydir*S/Zdir,S);
 
  Xdir = BeamDir.X();
  Ydir = BeamDir.Y();
  // POSITION/DIRECTION AT Z EMISSION// 
  X = X0 - Xdir;
  Y = Y0 - Ydir;
  Z = fZEmission;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Beam::Print() const {


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Beam::SetTargetSize(double TargetSize){
  fTargetSize = TargetSize;
  fEffectiveTargetSize = fTargetSize*cos(fTargetAngle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Beam::SetTargetThickness(double TargetThickness){
  fTargetThickness = TargetThickness;
  fEffectiveTargetThickness = fTargetThickness/cos(fTargetAngle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Beam::SetTargetAngle(double TargetAngle){
  fTargetAngle = TargetAngle;
  fEffectiveTargetSize = fTargetSize*cos(fTargetAngle);
  fEffectiveTargetThickness = fTargetThickness/cos(fTargetAngle);
}
