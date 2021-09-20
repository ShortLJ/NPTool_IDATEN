/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  F. Flavigny contact address: flavigny@lpccaen.in2p3.fr *
 *                                                                           *
 * Creation Date   : April 2019                                              *
 * Last update     : Nov 2019                                                *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with Quasi Free Scattering Reaction in which a cluster   *
 *  or a nucleon is removed from a projectile  by interaction with a target  *
 *  nucleon (proton target in general)                                       *
 *                                                                           *
 *  First step (dissociation):  A -> B + a                                   *
 *  Second step (scattering) :  a + T -> 1 + 2                               *
 *  Labeling is:                                                             *
 *                                                                           *
 *              A --> T  ==> B + (a -> T) =>  B + 1 + 2                      *
 *                                                                           *
 *  where:                                                                   *
 *    +  A is the beam nucleus                                               *
 *    +  T is the target nucleon (proton)                                    *
 *                                                                           *
 *    +  B is the residual fragment (beam-like)                              *
 *    +  1 is the scattered target nucleon  (former T)                       *
 *    +  2 is the knocked-out cluster/nucleon (noted a) in the intermediate  *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    +  Adapted from original event generator from V. Panin (R3B collab)    *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <vector>

#include "NPQFS.h"
#include "NPCore.h"
#include "NPOptionManager.h"
#include "NPFunction.h"

// Use CLHEP System of unit and Physical Constant
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

// ROOT
#include"TF1.h"

ClassImp(QFS)

////////////////////////////////////////////////////////////////////////////////

QFS::QFS(){

    //------------- Default Constructor -------------
    fVerboseLevel         = NPOptionManager::getInstance()->GetVerboseLevel();
    fBeamEnergy           = 0;
    fThetaCM              = 0;
    fPhiCM                = 0;
 
    fExcitationA          = 0;
    fExcitationB          = 0;
    fMomentumSigma        = 0;
    fInternalMomentum     = {0., 0.,0. };
    fshootB=false;
    fshoot1=true;
    fshoot2=true;
    fisotropic = true;

    fUseExInGeant4=true;

    fTheta2VsTheta1 = 0;
    fPhi2VsPhi1 = 0;

    fPerpMomentumHist = NULL;
    fParMomentumHist = NULL;
    fDeexcitation = false;
}

////////////////////////////////////////////////////////////////////////////////

QFS::~QFS(){

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void QFS::ReadConfigurationFile(string Path){
  ifstream ReactionFile;
  string GlobalPath = getenv("NPTOOL");
  string StandardPath = GlobalPath + "/Inputs/EventGenerator/" + Path;
  ReactionFile.open(Path.c_str());
  if (!ReactionFile.is_open()) {
    ReactionFile.open(StandardPath.c_str());
    if(ReactionFile.is_open()) {
      Path = StandardPath;
    }
    else {cout << "QFS File " << Path << " not found" << endl;exit(1);}
  }
  NPL::InputParser parser(Path);
  ReadConfigurationFile(parser);
}

////////////////////////////////////////////////////////////////////////////////

void QFS::ReadConfigurationFile(NPL::InputParser parser){

  cout << " In QFS ReadConfiguration " << endl;
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("QFSReaction");
  if(blocks.size()>0 && NPOptionManager::getInstance()->GetVerboseLevel())
      cout << endl << "\033[1;35m//// QFS reaction found " << endl;

  vector<string> token1 = {"Beam","Target","Scattered","KnockedOut","Heavy"};
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
      if(blocks[i]->HasTokenList(token1)){
          int v = NPOptionManager::getInstance()->GetVerboseLevel();
          NPOptionManager::getInstance()->SetVerboseLevel(0);
          fParticleA.ReadConfigurationFile(parser);
          NPOptionManager::getInstance()->SetVerboseLevel(v);

          fBeamEnergy= fParticleA.GetEnergy();
          GetParticle(blocks[i]->GetString("Beam"),parser);
          fParticleT = GetParticle(blocks[i]->GetString("Target"),parser);
          fParticleB = GetParticle(blocks[i]->GetString("Heavy"),parser);
          fParticle1 = GetParticle(blocks[i]->GetString("Scattered"),parser);
          fParticle2 = GetParticle(blocks[i]->GetString("KnockedOut"),parser);
      }
      else{
          cout << "ERROR: check your input file formatting \033[0m" << endl;
          exit(1);
      }
      if(blocks[i]->HasToken("ExcitationEnergyBeam")){
          fExcitationA = blocks[i]->GetDouble("ExcitationEnergyBeam","MeV");
      }
      if(blocks[i]->HasToken("ExcitationEnergyHeavy")){
          fExcitationB = blocks[i]->GetDouble("ExcitationEnergyHeavy","MeV");
      }
      if(blocks[i]->HasToken("MomentumSigma")){
          fMomentumSigma = blocks[i]->GetDouble("MomentumSigma","MeV");
      }
      if(blocks[i]->HasToken("ShootHeavy")){
          fshootB = blocks[i]->GetInt("ShootHeavy");
      }
      if(blocks[i]->HasToken("ShootLight")){
          fshoot1 = blocks[i]->GetInt("ShootLight");
          fshoot2 = blocks[i]->GetInt("ShootLight");
      }
      if(blocks[i]->HasToken("ShootLight1")){
          fshoot1 = blocks[i]->GetInt("ShootLight1");
      }
      if(blocks[i]->HasToken("ShootLight2")){
          fshoot2 = blocks[i]->GetInt("ShootLight2");
      }
      if(blocks[i]->HasToken("PerpMomentumPath")){
          vector<string> file_perp = blocks[i]->GetVectorString("PerpMomentumPath");
          TH1F* Perptemp = Read1DProfile(file_perp[0], file_perp[1]);
          SetPerpMomentumHist(Perptemp);
      }
      if(blocks[i]->HasToken("ParMomentumPath")){
          vector<string> file_par = blocks[i]->GetVectorString("ParMomentumPath");
          TH1F* Partemp = Read1DProfile(file_par[0], file_par[1]);
          SetParMomentumHist(Partemp);
      }
      if(blocks[i]->HasToken("Deexcitation")){
        fDeexcitation = blocks[i]->GetInt("Deexcitation");
      }
  }

  cout << "\033[0m" ;
}

////////////////////////////////////////////////////////////////////////////////
Particle QFS::GetParticle(string name, NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("DefineParticle",name);
  unsigned int size = blocks.size();
  if(size==0)
    return NPL::Particle(name);
  else if(size==1){
    cout << " -- User defined nucleus " << name << " -- " << endl;
    vector<string> token = {"SubPart","BindingEnergy"};
    if(blocks[0]->HasTokenList(token)){
      NPL::Particle N(name,blocks[0]->GetVectorString("SubPart"),blocks[0]->GetDouble("BindingEnergy","MeV"));
      if(blocks[0]->HasToken("ExcitationEnergy"))
        N.SetExcitationEnergy(blocks[0]->GetDouble("ExcitationEnergy","MeV"));
      if(blocks[0]->HasToken("SpinParity"))
        N.SetSpinParity(blocks[0]->GetString("SpinParity").c_str());
      if(blocks[0]->HasToken("Spin"))
        N.SetSpin(blocks[0]->GetDouble("Spin",""));
      if(blocks[0]->HasToken("Parity"))
        N.SetParity(blocks[0]->GetString("Parity").c_str());
      if(blocks[0]->HasToken("LifeTime"))
        N.SetLifeTime(blocks[0]->GetDouble("LifeTime","s"));

    cout << " -- -- -- -- -- -- -- -- -- -- --" << endl;
      return N;
    }
  }
  else{
    NPL::SendErrorAndExit("NPL::QFS","Too many nuclei define with the same name");
    return NPL::Particle();
  }
  return NPL::Particle();
}


////////////////////////////////////////////////////////////////////////////////////////////
void QFS::CalculateVariables(){

  if(fBeamEnergy < 0)
    fBeamEnergy = 0 ; 

    //cout<<"---- COMPUTE ------"<<endl;
   // cout<<"--CM--"<<endl; 

    mA =  fParticleA.Mass();           // Beam mass in MeV
    mT =  fParticleT.Mass();           // Target mass in MeV 
    fParticleB.SetExcitationEnergy(fExcitationB);
    mB =  fParticleB.Mass();           // Heavy residual mass in MeV 
    m1 =  mT;                        // scattered target nucleon (same mass);
    m2 =  fParticle2.Mass();           // knocked out cluster mass in MeV 
    ma =  m2;                        // intermediate cluster mass in MeV (same);
 
    double TA = fBeamEnergy;                 // Beam kinetic energy
    double PA = sqrt(TA*(TA+2*mA));          // Beam momentum (norm)
    double EA = sqrt(mA*mA + PA*PA);         // Beam total energy
    fEnergyImpulsionLab_A = TLorentzVector(0.,0.,PA,EA);
    
    // Internal momentum of removed cluster/nucleon (Pa) and recoil (PB)
    // here fInternalMomentum contains PB (recoil fragment momentum)
    // readout from the input file (theoretical)
    PB.SetX(fInternalMomentum.X());
    PB.SetY(fInternalMomentum.Y());
    PB.SetZ(fInternalMomentum.Z());
    Pa.SetXYZ( (-PB.X()) , (-PB.Y()) , (-PB.Z()) );

    // Off-shell mass of the bound nucleon from E conservation
    // in virtual dissociation of A -> B + a
    double buffer = mA*mA + mB*mB - 2*mA*sqrt(mB*mB+Pa.Mag2()) ; 
    if(buffer<=0) { cout<<"ERROR off shell mass ma_off=\t"<<buffer<<endl; return;}
    ma_off = sqrt(buffer);

    //deduced total energies of "a" and "B" in restframe of A
    double Ea = sqrt(ma_off*ma_off + Pa.Mag2());
    double EB = sqrt(mB*mB + PB.Mag2());

    fEnergyImpulsionCM_a = TLorentzVector(Pa,Ea);
    fEnergyImpulsionCM_B = TLorentzVector(PB,EB);

    fEnergyImpulsionLab_a = TLorentzVector(Pa,Ea);
    fEnergyImpulsionLab_B = TLorentzVector(PB,EB);
    fEnergyImpulsionLab_a.Boost(0,0,fEnergyImpulsionLab_A.Beta());
    fEnergyImpulsionLab_B.Boost(0,0,fEnergyImpulsionLab_A.Beta());
    Ea_lab = fEnergyImpulsionLab_a.E();
    EB_lab = fEnergyImpulsionLab_B.E();
    Pa_lab = fEnergyImpulsionLab_a.Vect();
    PB_lab = fEnergyImpulsionLab_B.Vect();
   
    // Scattering part (2-body kinematics)
    // virtual cluster of mass "ma_off" scattering on target T
    // to give scattered  cluster with real mass (ma=m2)
    // and scattered target (mT=m1)

    fQValue =ma_off+mT-m1-m2;

    s = ma_off*ma_off + mT*mT + 2*mT*Ea_lab ; 
    fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));
    fEcm = sqrt(s) - m1 -m2;
    if(fEcm<=0) { cout<<"ERROR Ecm negative =\t"<<fEcm<<endl;Dump(); return;}

    ECM_a = (s + ma_off*ma_off - mT*mT)/(2*sqrt(s));
    ECM_T = (s + mT*mT - ma_off*ma_off)/(2*sqrt(s));
    ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
    ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));

    pCM_a = sqrt(ECM_a*ECM_a - ma_off*ma_off);
    pCM_T = sqrt(ECM_T*ECM_T - mT*mT);
    pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
    pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);


}

////////////////////////////////////////////////////////////////////////////////////////////

void QFS::KineRelativistic(double& ThetaLab1, double& PhiLab1, double& KineticEnergyLab1, double& ThetaLab2, double& PhiLab2, double& KineticEnergyLab2){

    CalculateVariables();

    double thetaCM_1 = fThetaCM;
    double thetaCM_2 =  M_PI - thetaCM_1;
    double phiCM_2 = fPhiCM;

    TVector3 z_axis(0.,0.,1.);

/*  // OTHER WAY of doing
    TVector3 pCM_2_temp(0.,0.,1.);
    pCM_2_temp.SetMag(pCM_2);
    pCM_2_temp.SetTheta(thetaCM_2);
    pCM_2_temp.SetPhi(phiCM_2);
    fEnergyImpulsionCM_2	= TLorentzVector(pCM_2_temp,ECM_2);
    fEnergyImpulsionCM_1	= TLorentzVector(-1*pCM_2_temp,ECM_1);
*/
    fEnergyImpulsionCM_2	= TLorentzVector(
                                        pCM_2*sin(thetaCM_2)*cos(phiCM_2),
                                        pCM_2*sin(thetaCM_2)*sin(phiCM_2),
                                        pCM_2*cos(thetaCM_2),
                                        ECM_2);
  
    fEnergyImpulsionCM_1	= fTotalEnergyImpulsionCM - fEnergyImpulsionCM_2;

    //-- Boost in the direction of the moving cluster "a" --//
    BetaCM = Pa_lab.Mag() / (Ea_lab + mT);
    fEnergyImpulsionLab_1 = fEnergyImpulsionCM_1;
    fEnergyImpulsionLab_1.Boost(0,0,BetaCM);
    fEnergyImpulsionLab_2 = fEnergyImpulsionCM_2;
    fEnergyImpulsionLab_2.Boost(0,0,BetaCM);

    //-- Rotation to go from cluster frame to beam frame --//
    TVector3 direction = Pa_lab.Unit();
    fEnergyImpulsionLab_1.RotateUz(direction);
    fEnergyImpulsionLab_2.RotateUz(direction);
/*
    // Angle in the Lab frame
    ThetaLab1 = fEnergyImpulsionLab_1.Angle(fEnergyImpulsionLab_A.Vect());
    //ThetaLab1 = fEnergyImpulsionLab_1.Angle(z_axis);
    if (ThetaLab1 < 0) ThetaLab1 += M_PI;
    ThetaLab2 = fEnergyImpulsionLab_2.Angle(fEnergyImpulsionLab_A.Vect());
    //ThetaLab2 = fEnergyImpulsionLab_2.Angle(z_axis);
    if (fabs(ThetaLab1) < 1e-6) ThetaLab1 = 0;
    ThetaLab2 = fabs(ThetaLab2);
    if (fabs(ThetaLab2) < 1e-6) ThetaLab2 = 0;

    PhiLab1 = M_PI + fEnergyImpulsionLab_1.Vect().Phi(); 
    if (fabs(PhiLab1) < 1e-6) PhiLab1 = 0;
    PhiLab2 = M_PI + fEnergyImpulsionLab_2.Vect().Phi(); 
    if (fabs(PhiLab2) < 1e-6) PhiLab2 = 0;
*/

    ThetaLab1 = fEnergyImpulsionLab_1.Angle(fEnergyImpulsionLab_A.Vect());
    ThetaLab2 = fEnergyImpulsionLab_2.Angle(fEnergyImpulsionLab_A.Vect());
    PhiLab1 = fEnergyImpulsionLab_1.Vect().Phi(); 
    PhiLab2 = fEnergyImpulsionLab_2.Vect().Phi(); 

    // Kinetic Energy in the lab frame
    KineticEnergyLab1 = fEnergyImpulsionLab_1.E() - m1;
    KineticEnergyLab2 = fEnergyImpulsionLab_2.E() - m2;
    // test for total energy conversion
    //if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_1.E()+fEnergyImpulsionLab_2.E())) > 1e-6)
    //    cout << "Problem for energy conservation" << endl;
    
    //Dump();

}

////////////////////////////////////////////////////////////////////////////////////////////
void QFS::Dump(){

    cout<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<"------------ DUMP QFS --------------"<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<endl;
    cout<<"Cluster/recoil momentum (in the beam nucleus frame):"<<endl; 
    cout<<"Pa=\t("<<Pa.Px()<<","<<Pa.Py()<<","<<Pa.Pz()<<") MeV/c"<<endl;
    cout<<"PB=\t("<<PB.Px()<<","<<PB.Py()<<","<<PB.Pz()<<") MeV/c"<<endl;
    cout<<endl;
    cout<<"Off-shell mass of the bound nucleon from E conservation "<<endl;
    cout<<" in virtual dissociation of A -> B + a"<<endl;
    cout<<"ma=\t"<<ma<<endl;
    cout<<"ma_off=\t"<<ma_off<<endl;
    cout<<"mB=\t"<<mB<<endl;
    cout<<"Deduced total energies of a and B in restframe of A"<<endl;
    cout<<"Ea=\t"<<fEnergyImpulsionCM_a.E()<<"\tMeV"<<endl;
    cout<<"EB=\t"<<fEnergyImpulsionCM_B.E()<<"\tMeV"<<endl;
    cout<<endl;
    cout<<"-- Boosted in lab frame with beam on Z axis --"<<endl; 
    cout<<"Beta_z=\t"<<fEnergyImpulsionLab_A.Beta()<<endl;
    cout<<"Pa_lab=\t("<<Pa_lab.Px()<<","<<Pa_lab.Py()<<","<<Pa_lab.Pz()<<") MeV/c"<<endl;
    cout<<"PB_lab=\t("<<PB_lab.Px()<<","<<PB_lab.Py()<<","<<PB_lab.Pz()<<") MeV/c"<<endl;
    cout<<"Ea_lab=\t"<<Ea_lab<<"\tMeV"<<endl;
    cout<<"EB_lab=\t"<<EB_lab<<"\tMeV"<<endl;
    cout<<endl; 
    cout<<"-- Scattering off virtual cluster a of virtual mass --"<<endl; 
    cout<<"-- ma_off and energy Ea_lab on target T at rest ------"<<endl;
    cout<<"fQValue=\t"<<fQValue<<endl;
    cout<<"s=\t"<<s<<endl;
    cout<<"Ecm=\t"<<fEcm<<endl;
    cout<<"ea*=\t"<<ECM_a<<endl;
    cout<<"pa*=\t"<<pCM_a<<endl;
    cout<<"eT*=\t"<<ECM_T<<endl;
    cout<<"pT*=\t"<<pCM_T<<endl;
    cout<<"e1*=\t"<<ECM_1<<endl;
    cout<<"p1*=\t"<<pCM_1<<endl;
    cout<<"e2*=\t"<<ECM_2<<endl;
    cout<<"p2*=\t"<<pCM_2<<endl;
    cout<<"beta_cm=\t"<<BetaCM<<endl;
    cout<<endl;
    cout<<"-- Emitted Particles --"<<endl;
    cout<<"Theta_cm:"<<fThetaCM*180./TMath::Pi()<<endl;
    cout<<"Phi_cm:"<<fPhiCM*180./TMath::Pi()<<endl;
    cout<<"P1_CM=\t("<<fEnergyImpulsionCM_1.Px()<<","<<fEnergyImpulsionCM_1.Py()<<","<<fEnergyImpulsionCM_1.Pz()<<")"<<endl;
    cout<<"P2_CM=\t("<<fEnergyImpulsionCM_2.Px()<<","<<fEnergyImpulsionCM_2.Py()<<","<<fEnergyImpulsionCM_2.Pz()<<")"<<endl;
    cout<<"E1_lab=\t"<<fEnergyImpulsionLab_1.E()<<endl;
    cout<<"E2_lab=\t"<<fEnergyImpulsionLab_2.E()<<endl;
    cout<<"P1_lab=\t("<<fEnergyImpulsionLab_1.Px()<<","<<fEnergyImpulsionLab_1.Py()<<","<<fEnergyImpulsionLab_1.Pz()<<")"<<endl;
    cout<<"P2_lab=\t("<<fEnergyImpulsionLab_2.Px()<<","<<fEnergyImpulsionLab_2.Py()<<","<<fEnergyImpulsionLab_2.Pz()<<")"<<endl;
    cout<<"Theta1:\t"<<fEnergyImpulsionLab_1.Vect().Theta()*180./TMath::Pi()<<endl;
    cout<<"Theta2:\t"<<fEnergyImpulsionLab_2.Vect().Theta()*180./TMath::Pi()<<endl;
    cout<<"Phi1:\t"<<fEnergyImpulsionLab_1.Vect().Phi()*180./TMath::Pi()<<endl;
    cout<<"Phi2:\t"<<fEnergyImpulsionLab_2.Vect().Phi()*180./TMath::Pi()<<endl;

}
////////////////////////////////////////////////////////////////////////////////////////////

TVector3 QFS::ShootInternalMomentum(){

  //Shoot a momentum (vector) for the internal cluster in the beam-at-rest frame
  // (1) if only a width is provided: shoot in 3 independant Gaussian
  // (2) if input histos are provided: use them instead of option (1)  
  // Remark : if both width and input histos are provided only histos are considered 

  TVector3 momentum = {0,0,0};
  double  PerpMomentum =0;
  double  ParMomentum =0;
  double  angle_tmp =0;

  momentum.SetX(gRandom->Gaus(0.,fMomentumSigma));
  momentum.SetY(gRandom->Gaus(0.,fMomentumSigma));
  momentum.SetZ(gRandom->Gaus(0.,fMomentumSigma));

  if(fPerpMomentumHist){
      PerpMomentum=fPerpMomentumHist->GetRandom();
      angle_tmp = gRandom->Rndm()*2*M_PI;
      momentum.SetX(PerpMomentum * TMath::Cos(angle_tmp));
      momentum.SetY(PerpMomentum * TMath::Sin(angle_tmp));
  }
  if(fParMomentumHist){
      ParMomentum=fParMomentumHist->GetRandom();
      momentum.SetZ(ParMomentum);
  }

  //cout << " Shooting Random Momentum: "  << endl;
  //cout<<"Px:"<<momentum.X() << endl;
  //cout<<"Py:"<<momentum.Y() << endl;
  //cout<<"Pz:"<<momentum.Z() << endl;
  SetInternalMomentum(momentum);
  return momentum;
}

////////////////////////////////////////////////////////////////////////////////////////////

TGraph* QFS::GetTheta2VsTheta1(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta1,phi1,E1,theta2,phi2,E2;
  SetPhiCM(0.*TMath::Pi()/180.);

  for (double angle=0 ; angle <= 180 ; angle+=AngleStep_CM){
    SetThetaCM(angle*TMath::Pi()/180.);
    KineRelativistic(theta1, phi1, E1, theta2, phi2, E2);
    vx.push_back(theta1*180./M_PI);
    vy.push_back(theta2*180./M_PI);
  }
  fTheta2VsTheta1 = new TGraph(vx.size(),&vx[0],&vy[0]);

  return(fTheta2VsTheta1);
}

////////////////////////////////////////////////////////////////////////////////////////////

TGraph* QFS::GetPhi2VsPhi1(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta1,phi1,E1,theta2,phi2,E2;
  SetThetaCM(0.*TMath::Pi()/180.);

  for (double angle=-180 ; angle <= 180 ; angle+=AngleStep_CM){
      SetPhiCM(angle*TMath::Pi()/180.);
      KineRelativistic(theta1, phi1, E1, theta2, phi2, E2);
      vx.push_back(phi1*180./M_PI);
      vy.push_back(phi2*180./M_PI);
  }
  fPhi2VsPhi1 = new TGraph(vx.size(),&vx[0],&vy[0]);

  return(fPhi2VsPhi1);
}

///////////////////////////////////////////////////////////////////////////////
// Check whenever the reaction is allowed at a given energy
bool QFS::IsAllowed(){//double Energy){
  //double AvailableEnergy = Energy + fParticle1.Mass() + fParticle2.Mass();
  //double RequiredEnergy  = fParticle3.Mass() + fParticle4.Mass();

  //if(AvailableEnergy>RequiredEnergy)
    return true;
  //else
  //  return false;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
///////////////// Old R3B method not using TLorentz Vector  ////////////////////////////////
/////////////////// (used as a reference)///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void QFS::TestR3B()
{
CalculateVariablesOld();
}
void QFS::CalculateVariablesOld(){

  if(fBeamEnergy < 0)
    fBeamEnergy = 0 ; 

    //cout<<"---- COMPUTE ------"<<endl;
   // cout<<"--CM--"<<endl; 

    mA =  fParticleA.Mass();            // Beam mass in MeV
    mT =  fParticleT.Mass();           // Target mass in MeV 
    fParticleB.SetExcitationEnergy(fExcitationB);
    mB =  fParticleB.Mass();           // Heavy residual mass in MeV 
    m1 =  mT;                        // scattered target nucleon (same mass);
    m2 =  fParticle2.Mass();           // knocked out cluster mass in MeV 
    ma =  m2;                        // intermediate cluster mass in MeV (same);
 
    double TA = fBeamEnergy;                 // Beam kinetic energy
    double PA = sqrt(TA*(TA+2*mA));          // Beam momentum (norm)
    double EA = sqrt(mA*mA + PA*PA);         // Beam total energy
    fEnergyImpulsionLab_A = TLorentzVector(0.,0.,PA,EA);
    
    //Internal momentum of removed cluster/nucleon
    //gRandom->SetSeed(0);
    //double mom_sigma = 0; // MeV/c
    //Pa.SetX(gRandom->Gaus(0.,fMomentumSigma));
    //Pa.SetY(gRandom->Gaus(0.,fMomentumSigma));
    //Pa.SetZ(gRandom->Gaus(0.,fMomentumSigma));
    Pa.SetX(50);
    Pa.SetY(50);
    Pa.SetZ(50);



    //Internal momentum of heavy recoil after removal
    PB.SetXYZ( (-Pa.X()) , (-Pa.Y()) , (-Pa.Z()) );

    // Off-shell mass of the bound nucleon from E conservation
    // in virtual dissociation of A -> B + a
    double buffer = mA*mA + mB*mB - 2*mA*sqrt(mB*mB+Pa.Mag2()) ; 
    if(buffer<=0) { cout<<"ERROR off shell mass ma_off=\t"<<buffer<<endl; return;}
    ma_off = sqrt(buffer);

    //deduced total energies of "a" and "B" in restframe of A
    double Ea = sqrt(ma_off*ma_off + Pa.Mag2());
    double EB = sqrt(mB*mB + PB.Mag2());

    fEnergyImpulsionCM_a = TLorentzVector(Pa,Ea);
    fEnergyImpulsionCM_B = TLorentzVector(PB,EB);

    fEnergyImpulsionLab_a = TLorentzVector(Pa,Ea);
    fEnergyImpulsionLab_B = TLorentzVector(PB,EB);
    fEnergyImpulsionLab_a.Boost(0,0,fEnergyImpulsionLab_A.Beta());
    fEnergyImpulsionLab_B.Boost(0,0,fEnergyImpulsionLab_A.Beta());
    Ea_lab = fEnergyImpulsionLab_a.E();
    EB_lab = fEnergyImpulsionLab_B.E();
    Pa_lab = fEnergyImpulsionLab_a.Vect();
    PB_lab = fEnergyImpulsionLab_B.Vect();

   
    // Scattering part (2-body kinematics)
    // virtual cluster of mass "ma_off" scattering on target T
    // to give scattered  cluster with real mass (ma=m2)
    // and scattered target (mT=m1)

    fQValue =ma_off+mT-m1-m2;

    s = ma_off*ma_off + mT*mT + 2*mT*Ea_lab ; 
    fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));
    fEcm = sqrt(s) - m1 -m2;
    if(fEcm<=0) { cout<<"ERROR Ecm negative =\t"<<fEcm<<endl;Dump(); return;}

    vector<double> theta1;
    vector<double> theta2;
    vector<double> phi1;
    vector<double> phi2;

    //for(int i=0; i<=180; i++){
    int i = 30;
        KineR3B(s, ma_off, mT, ma, (double)i);
        if(!good) { cout<<"ERROR CM calculations!!!"<<endl; return;}

        cout<<endl;
        cout<<"------------------------------------"<<endl;
        cout<<"------------ DUMP R3B --------------"<<endl;
        cout<<"------------------------------------"<<endl;
        cout<<endl;
        cout<<"Cluster/recoil momentum (in the beam nucleus frame):"<<endl; 
        cout<<"Pa=\t("<<Pa.Px()<<","<<Pa.Py()<<","<<Pa.Pz()<<") MeV/c"<<endl;
        cout<<"PB=\t("<<PB.Px()<<","<<PB.Py()<<","<<PB.Pz()<<") MeV/c"<<endl;
        cout<<endl;
        cout<<"Off-shell mass of the bound nucleon from E conservation "<<endl;
        cout<<" in virtual dissociation of A -> B + a"<<endl;
        cout<<"ma=\t"<<ma<<endl;
        cout<<"ma_off=\t"<<ma_off<<endl;
        cout<<"mB=\t"<<mB<<endl;
        cout<<"Deduced total energies of a and B in restframe of A"<<endl;
        cout<<"Ea=\t"<<fEnergyImpulsionCM_a.E()<<"\tMeV"<<endl;
        cout<<"EB=\t"<<fEnergyImpulsionCM_B.E()<<"\tMeV"<<endl;
        cout<<endl;
        cout<<"-- Boosted in lab frame with beam on Z axis --"<<endl; 
        cout<<"Beta_z=\t"<<fEnergyImpulsionLab_A.Beta()<<endl;
        cout<<"Pa_lab=\t("<<Pa_lab.Px()<<","<<Pa_lab.Py()<<","<<Pa_lab.Pz()<<") MeV/c"<<endl;
        cout<<"PB_lab=\t("<<PB_lab.Px()<<","<<PB_lab.Py()<<","<<PB_lab.Pz()<<") MeV/c"<<endl;
        cout<<"Ea_lab=\t"<<Ea_lab<<"\tMeV"<<endl;
        cout<<"EB_lab=\t"<<EB_lab<<"\tMeV"<<endl;
        cout<<endl; 
        cout<<"-- Scattering off virtual cluster a of virtual mass --"<<endl; 
        cout<<"-- ma_off and energy Ea_lab on target T at rest ------"<<endl;
        cout<<"fQValue=\t"<<fQValue<<endl;
        cout<<"s=\t"<<s<<endl;
        cout<<"Ecm=\t"<<fEcm<<endl;
        cout<<endl;


        TVector3 P1_cm(0.,0.,1.), P2_cm(0.,0.,1.);
        P2_cm.SetMag(p_clust);
        P2_cm.SetTheta(theta_clust);
        //TRandom3 ra;
        //ra.SetSeed(0);
        //P2_cm.SetPhi(ra.Uniform(-1*TMath::Pi(),+1*TMath::Pi()));
        P2_cm.SetPhi(20.*TMath::Pi()/180.);
        P1_cm.SetX(-P2_cm.X());
        P1_cm.SetY(-P2_cm.Y());
        P1_cm.SetZ(-P2_cm.Z());

        cout<<"P1_CM=\t("<<P1_cm.X()<<","<<P1_cm.Y()<<","<<P1_cm.Z()<<")"<<endl;
        cout<<"P2_CM=\t("<<P2_cm.X()<<","<<P2_cm.Y()<<","<<P2_cm.Z()<<")"<<endl;
 
        // Calculate relative to direction of quasi-particle (cluster)
        
        double beta_cm = -Pa_lab.Mag() / (Ea_lab + mT);
        double gamma_cm = 1/sqrt(1-beta_cm*beta_cm);

        pair<double,double> lor_a1 = Lorentz(gamma_cm,beta_cm,e_scat,P1_cm.Z());
        pair<double,double> lor_a2 = Lorentz(gamma_cm,beta_cm,e_clust,P2_cm.Z());

        P1_cm.SetZ(lor_a1.second);
        P2_cm.SetZ(lor_a2.second);

        //Rotating back to beam direction
        TVector3 P1_L = Rotations(P1_cm, Pa_lab);
        TVector3 P2_L = Rotations(P2_cm, Pa_lab);
        
        //TVector3 P1_L = P1_cm;
        //TVector3 P2_L = P2_cm;
        //TVector3 direction = Pa.Unit();
        //P1_L.RotateUz(direction);
        //P1_L.RotateUz(direction);

        cout<<"----Calculate variables output------"<<endl;
       cout<<"--CM--"<<endl;
        cout<<"theta1*=\t"<<theta_scat*180/TMath::Pi()<<endl;
        cout<<"theta2*=\t"<<theta_clust*180/TMath::Pi()<<endl;
        cout<<"e1*=\t"<<e_scat<<endl;
        cout<<"p1*=\t"<<p_scat<<endl;
        cout<<"e2*=\t"<<e_clust<<endl;
        cout<<"p2*=\t"<<p_clust<<endl;
        cout<<"T=\t"<<T<<endl;
        cout<<"beta_cm=\t"<<beta_cm<<endl;
        cout<<"gamma_cm=\t"<<gamma_cm<<endl;

        cout<<"--LAB (cluster dir)--"<<endl;
        cout<<"P1_lab=\t("<<P1_cm.X()<<","<<P1_cm.Y()<<","<<P1_cm.Z()<<")"<<endl;
        cout<<"P2_lab=\t("<<P2_cm.X()<<","<<P2_cm.Y()<<","<<P2_cm.Z()<<")"<<endl;
        cout<<"Theta1:\t"<<P1_cm.Theta()*180./TMath::Pi()<<endl;
        cout<<"Theta2:\t"<<P2_cm.Theta()*180./TMath::Pi()<<endl;

        cout<<"--LAB--"<<endl;
        cout<<"Pa_lab=\t("<<Pa_lab.X()<<","<<Pa_lab.Y()<<","<<Pa_lab.Z()<<")"<<endl;
        cout<<"P1_L=\t("<<P1_L.X()<<","<<P1_L.Y()<<","<<P1_L.Z()<<")"<<endl;
        cout<<"P2_L=\t("<<P2_L.X()<<","<<P2_L.Y()<<","<<P2_L.Z()<<")"<<endl;
        cout<<"Theta1L:\t"<<P1_L.Theta()*180./TMath::Pi()<<endl;
        cout<<"Theta2L:\t"<<P2_L.Theta()*180./TMath::Pi()<<endl;
        cout<<"Phi1L:\t"<<P1_L.Phi()*180./TMath::Pi()<<endl;
        cout<<"Phi2L:\t"<<P2_L.Phi()*180./TMath::Pi()<<endl;


        //cout<<P1_cm.Theta()*180./TMath::Pi()<<"\t"<<P2_cm.Theta()*180./TMath::Pi()<<endl;
        //cout<<P1_L.Phi()*180./TMath::Pi()<<"\t"<<P2_L.Phi()*180./TMath::Pi()<<endl;
        
       theta1.push_back(P1_L.Theta()*180./TMath::Pi());
       theta2.push_back(P2_L.Theta()*180./TMath::Pi());
      
       double temp_phi1 = P1_L.Phi(); 
       double temp_phi2 = P2_L.Phi(); 
       phi1.push_back(180. + P1_L.Phi()*180./TMath::Pi());
       phi2.push_back(180. + P2_L.Phi()*180./TMath::Pi());

    //}
    TGraph* fTheta2VsTheta1 = new TGraph(theta1.size(),&theta1[0],&theta2[0]);
    TGraph* fPhi2VsPhi1 = new TGraph(phi1.size(),&phi1[0],&phi2[0]);
    fTheta2VsTheta1->SetName("Theta2VsTheta1");
    fPhi2VsPhi1->SetName("Phi2VsPhi1");
    TFile* f = new TFile("graphs.root","RECREATE");
    fTheta2VsTheta1->Write();
    fPhi2VsPhi1->Write();
    f->Close();


}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
void QFS::KineR3B(double s,double m2off,double m1,double m2,double thetacm)
{
     if(thetacm>180 || thetacm<0){
        cout << "\nERROR! ThetaCM (in deg) should be between 0 and 180"<<endl;
        return;
    }

	e_clust = 0;
	p_clust = 0;
	theta_clust = 0;
	e_scat = 0;
	p_scat = 0;
	theta_scat = 0;
	T = 0;
	good = false;


	double X = s;
	double Y = m2off*m2off;
	double Z = m1*m1;
	double sqrt_s = sqrt(s);

	// Kinematics before the scattering process
	// (with one off-shell mass)
	double p2_off = sqrt(function(X,Y,Z))/2/sqrt_s;
	double p1_off = p2_off;
	// CM energies
	double e1_off = (s+Z-Y)/2/sqrt_s;
	double e2_off = (s+Y-Z)/2/sqrt_s;

	// Now take the real masses (after scattering)
	Y = m2*m2;  Z = m1*m1;
	//And check whether the kinematical function is ok
	//for this specific kinematical case
	double ERROR_CI = function(X,Y,Z);
	if(ERROR_CI <= 0.){
		cout << "\nERROR!!! Kinematical function is negative!";
		return;
	}

	// Kinematics after the scattering process
	// (with all real masses)
	double p2 = sqrt(function(X,Y,Z))/2/sqrt_s;
	double p1 = p2;
	double e1 = (s+Z-Y)/2/sqrt_s;
	double e2 = (s+Y-Z)/2/sqrt_s;

	// Let's consider momentum transfer <t> from the
	// target particle 1 to the cluster 2
    double t = 2*(m1*m1 - e1_off*e1 + p1_off*p1*cos(thetacm*TMath::Pi()/180.));

    //CM scattering angles
    double theta1 = thetacm*TMath::Pi()/180.;
    double theta2 = TMath::Pi() - theta1;

    e_clust = e2;
    p_clust = p2;
    theta_clust = theta2;

    e_scat = e1;
    p_scat = p1;
    theta_scat = theta1;

    T = t;
    good = true;


}



double QFS::function(double x,double y,double z)
{	
	double lambda = x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
	return lambda;
}

//---- Two consecutive rotations 
//first around Z on <phi>, then around new X' on <theta> (1=Pcm, 2=Pa in lab)
TVector3 QFS::Rotations(TVector3 v1,TVector3 v2) 
{
	double CT = v2.Z()/v2.Mag(); // cos(theta) of v2 wrt. Z-axis
	double ST = sqrt(1-CT*CT);   // sin(theta)
	double CF = v2.X()/v2.Mag()/ST;
	double SF = v2.Y()/v2.Mag()/ST;

	TVector3 v3;
	double _v3x =  v1.X()*CT*CF - v1.Y()*SF + v1.Z()*ST*CF;
	double _v3y =  v1.X()*CT*SF + v1.Y()*CF + v1.Z()*ST*SF;
	double _v3z = -v1.X()*ST   +  v1.Z()*CT;
	v3.SetXYZ(_v3x,_v3y,_v3z);

    cout<<"--- ROTATION---"<<endl;
    cout<<"CT=\t"<<CT<<endl;
    cout<<"ST=\t"<<ST<<endl;
    cout<<"CF=\t"<<CF<<endl;
    cout<<"SF=\t"<<SF<<endl;
    cout<<"v3x=\t"<<_v3x<<endl;
    cout<<"v3y=\t"<<_v3y<<endl;
    cout<<"v3z=\t"<<_v3z<<endl;
	return v3;
}



pair<double, double> QFS::Lorentz(double gamma,double beta,double e,double p)
{
	double eL = gamma*e - gamma*beta*p;
	double pL = gamma*p - gamma*beta*e;
	return make_pair(eL, pL);
}

