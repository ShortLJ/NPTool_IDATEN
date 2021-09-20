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
 * Creation Date   : March 2009                                              *
 * Last update     : January 2011                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with Two Body transfert Reaction                         *
 *  Physical parameter (Particle mass) are loaded from the nubtab03.asc file   *
 *  (2003 nuclear table of isotopes mass).                                   *
 *                                                                           *
 *  KineRelativistic: Used in NPSimulation                                   *
 *  A relativistic calculation is made to compute Light and Heavy nuclei     *
 *  angle given the Theta CM angle.                                          *
 *                                                                           *
 *  ReconstructRelativistic: Used in NPAnalysis                              *
 *  A relativistic calculation is made to compute Excitation energy given the*
 *  light angle and energy in Lab frame.                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + Redesign using LorentzVector by Pierre Morfouace                     *
 *    + 20/01/2011: Add support for excitation energy for light ejectile     *
 *                  (N. de Sereville)                                        *
 *    + Based on previous work by N.de Sereville                             *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <vector>

#include "NPReaction.h"
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

ClassImp(Reaction)

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  Reaction::Reaction(){
    //------------- Default Constructor -------------

    // Need to be done before initializePrecomputeVariable
    fKineLine3 = 0 ;
    fKineLine4 = 0 ;
    fLineBrho3 = 0 ;
    fTheta3VsTheta4 = 0;
    fAngleLine = 0;

    //
    fBeamEnergy           = 0;
    fThetaCM              = 0;
    fExcitation1          = 0;
    fExcitation3          = 0;
    fExcitation4          = 0;
    fQValue               = 0;
    fVerboseLevel         = NPOptionManager::getInstance()->GetVerboseLevel();
    initializePrecomputeVariable();

    fCrossSectionHist = NULL;
    fExcitationEnergyHist = NULL;
    fDoubleDifferentialCrossSectionHist = NULL ;

    fshoot3=true;
    fshoot4=true;
    fUseExInGeant4=true;

    fLabCrossSection=false; // flag if the provided cross-section is in the lab or not

  }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// This constructor aim to provide a fast way to instantiate a reaction without input file
// The string should be of the form "A(b,c)D@E" with E the ernegy of the beam in MeV
Reaction::Reaction(string reaction){
  // Instantiate the parameter to default
  // Analyse the reaction and extract revelant information
  string A,b,c,D,E;
  unsigned int i=0;
  for(; i < reaction.length() ; i++){
    if(reaction.compare(i,1,"(")!=0) A.push_back(reaction[i]);
    else break;
  }

  i++;
  for(; i < reaction.length() ; i++){
    if(reaction.compare(i,1,",")!=0) b.push_back(reaction[i]);
    else break;
  }

  i++;
  for(; i < reaction.length() ; i++){
    if(reaction.compare(i,1,")")!=0) c.push_back(reaction[i]);
    else break;
  }

  i++;
  for(; i < reaction.length() ; i++){
    if(reaction.compare(i,1,"@")!=0) D.push_back(reaction[i]);
    else break;
  }

  i++;
  for(; i < reaction.length() ; i++){
    E.push_back(reaction[i]);
  }

  fKineLine3 = 0 ;
  fKineLine4 = 0 ;
  fLineBrho3 = 0;
  fTheta3VsTheta4 = 0;
  fAngleLine = 0;
  fParticle1 = Beam(A);
  fParticle2 = Particle(b);
  fParticle3 = Particle(c);
  fParticle4 = Particle(D);
  fBeamEnergy = atof(E.c_str());
  fThetaCM              = 0;
  fExcitation1          = 0;
  fExcitation3          = 0;
  fExcitation4          = 0;
  fQValue               = 0;
  fVerboseLevel         = NPOptionManager::getInstance()->GetVerboseLevel();
  initializePrecomputeVariable();

  // do that to avoid warning from multiple Hist with same name...  int offset = 0;
  int offset = 0;
  while(gDirectory->FindObjectAny(Form("EnergyHist_%i",offset))!=0)
    ++offset;

  fCrossSectionHist = new TH1F(Form("EnergyHist_%i",offset),"Reaction_CS",1,0,180);
  fDoubleDifferentialCrossSectionHist = NULL ;

  fshoot3=true;
  fshoot4=true;

  fLabCrossSection=false;

  initializePrecomputeVariable();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Reaction::~Reaction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool Reaction::CheckKinematic(){
  double theta = fThetaCM;
  if (m1 > m2) theta = M_PI - fThetaCM;
  fEnergyImpulsionCM_3	= TLorentzVector(pCM_3*sin(theta),0,pCM_3*cos(theta),ECM_3);
  fEnergyImpulsionCM_4	= fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  if ( fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E()))> 1e-6){
    cout << "Problem with energy conservation" << endl;
    return false;
  }
  else{
    //cout << "Kinematic OK" << endl;
    return true;
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double Reaction::ShootRandomThetaCM(){
  double theta; // CM
  if(fDoubleDifferentialCrossSectionHist){
    // Take a slice in energy
    TAxis* Y = fDoubleDifferentialCrossSectionHist->GetYaxis();
    int binY;

    // Those test are there for the tail event of the energy distribution
    // In case the energy is outside the range of the 2D histo we take the
    // closest available CS
    if(Y->FindBin(fBeamEnergy) > Y->GetLast())
      binY = Y->GetLast()-1;

    else if(Y->FindBin(fBeamEnergy) < Y->GetFirst())
      binY = Y->GetFirst()+1;

    else
      binY = Y->FindBin(fBeamEnergy);

    TH1D* Proj = fDoubleDifferentialCrossSectionHist->ProjectionX("proj",binY,binY);
    SetThetaCM( theta=Proj->GetRandom()*deg );
  }
  else if (fLabCrossSection){
    double thetalab=-1;
    double energylab=-1;
    while(energylab<0){
      thetalab=fCrossSectionHist->GetRandom()*deg; //shoot in lab
      energylab=EnergyLabFromThetaLab(thetalab);   //get corresponding energy
    }
    theta = EnergyLabToThetaCM(energylab, thetalab); //transform to theta CM
    SetThetaCM( theta );
  }
  else{
    // When root perform a Spline interpolation to shoot random number out of
    // the distribution, it can over shoot and output a number larger that 180
    // this lead to an additional signal at 0-4 deg Lab, especially when using a
    // flat distribution.
    // This fix it.
    theta=181;
    if(theta/deg>180)
      theta=fCrossSectionHist->GetRandom();
    //cout << " Shooting Random ThetaCM "  << theta << endl;
    SetThetaCM( theta*deg );
  }


  return theta;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Reaction::ShootRandomExcitationEnergy(){
  if(fExcitationEnergyHist){
    SetExcitation4(fExcitationEnergyHist->GetRandom());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Reaction::KineRelativistic(double& ThetaLab3, double& KineticEnergyLab3,
    double& ThetaLab4, double& KineticEnergyLab4){
  // 2-body relativistic kinematics: direct + inverse
  // EnergieLab3,4 : lab energy in MeV of the 2 ejectiles
  // ThetaLab3,4   : angles in rad
  // case of inverse kinematics

  double theta = fThetaCM;
  if (m1 > m2) theta = M_PI - fThetaCM;

  fEnergyImpulsionCM_3	= TLorentzVector(pCM_3*sin(theta),0,pCM_3*cos(theta),ECM_3);
  fEnergyImpulsionCM_4	= fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;

  fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
  fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
  fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
  fEnergyImpulsionLab_4.Boost(0,0,BetaCM);

  // Angle in the lab frame
  ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
  if (ThetaLab3 < 0) ThetaLab3 += M_PI;

  ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
  if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
  ThetaLab4 = fabs(ThetaLab4);
  if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;

  // Kinetic Energy in the lab frame
  KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - m3;
  KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - m4;

  // test for total energy conversion
  if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
    cout << "Problem for energy conservation" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double Reaction::ReconstructRelativistic(double EnergyLab, double ThetaLab){
  // EnergyLab in MeV
  // ThetaLab in rad
  double E3 = m3 + EnergyLab;
  double p_Lab_3 = sqrt(E3*E3 - m3*m3);
  fEnergyImpulsionLab_3 = TLorentzVector(p_Lab_3*sin(ThetaLab),0,p_Lab_3*cos(ThetaLab),E3);
  fEnergyImpulsionLab_4 = fTotalEnergyImpulsionLab - fEnergyImpulsionLab_3;


  double Eex = fEnergyImpulsionLab_4.Mag() - fParticle4.Mass();

  return Eex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Return ThetaCM
double  Reaction::EnergyLabToThetaCM(double EnergyLab, double ThetaLab){
  double E3 = m3 + EnergyLab;
  double p_Lab_3 = sqrt(E3*E3 - m3*m3);

  fEnergyImpulsionLab_3 = TLorentzVector(p_Lab_3*sin(ThetaLab),0,p_Lab_3*cos(ThetaLab),E3);
  fEnergyImpulsionCM_3 = fEnergyImpulsionLab_3;
  fEnergyImpulsionCM_3.Boost(0,0,-BetaCM);

  double ThetaCM = M_PI - fEnergyImpulsionCM_1.Angle(fEnergyImpulsionCM_3.Vect());
  
  return(ThetaCM);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Return EnergyLab
double  Reaction::EnergyLabFromThetaLab(double ThetaLab){
  //ThetaLab in rad
 
  // IMPORTANT NOTICE: This function is not suitable for reaction A(c,d)B
  // where M(c) < M(d), i.e. (p,d) or (d,3He) since in this case the same 
  // angle observed in the lab corresponds to two different energies.
  //
  // If this situation is encountered here:
  // 1- the final energy will be determined randomly with a 50%-50% 
  //     split for either energies.
  // 2- Any angle outside of the allowed range (in this case) will return -1.
  // 3- The user lab distribution will be used between angles {0,max} and the
  //    spatial distribution of the high energy and low energy branch
  //    will be the same. 
  // 4- Practically, in an experiment using a downstream spectrometer 
  //    only one of the energy branches is considered. 
  //
  // !! If both of the branches are needed the user SHOULD use the center-of-mass distribution.

  //
  // This calculation uses the formalism from J.B Marion and F.C Young 
  // (Book: Nucler Reaction Analysis, Graphs and Tables)

  //Treat Exeptions
  if(fBeamEnergy==0) return m4*fQValue/(m3+m4);

  
  double A,B,C,D, ThetaLabMax=181*deg;
  double Q = fQValue ;
  double T1 = fBeamEnergy;
  double TT = fBeamEnergy+Q;
  double sign=+1;
  
  A = m1*m4*(T1/TT) /(m1+m2)/(m3+m4);
  B = m1*m3*(T1/TT) /(m1+m2)/(m3+m4);
  C = m2*m3*(1+ (m1*Q)/(m2*TT)) /(m1+m2)/(m3+m4);
  D = m2*m4*(1+ (m1*Q)/(m2*TT)) /(m1+m2)/(m3+m4);

  if ( fabs(A+B+C+D-1)>1e-6  or  fabs(A*C-B*D)>1e-6 ) {
     cout << " Reaction balance is wrong in NPReaction object." << endl;
     exit(-1);
  }
 
  if(B>D) {
    ThetaLabMax = asin(sqrt(D/B));
    if(gRandom->Rndm()<0.5) sign=-1;
  }
  if(ThetaLab>ThetaLabMax) return -1;

  double cosine = cos(ThetaLab);
  double sine2 = pow(sin(ThetaLab),2); 
  double factor = sqrt(D/B-sine2);
  double EnergyLab = TT * B * pow( cosine + sign*factor , 2 );

  //cout << " Angle/energy: " << ThetaLab/deg << " " << EnergyLab << endl ; 
  //cin.get();

  return EnergyLab ;  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Reaction::Print() const{
  // Print informations concerning the reaction

  cout << "Reaction : " << fParticle2.GetName() << "(" << fParticle1.GetName()
    << "," << fParticle3.GetName() << ")" << fParticle4.GetName() << "  @  "
    << fBeamEnergy << " MeV"
    << endl   ;

  cout << "Exc Particle 1 = " << fExcitation1 << " MeV" << endl;
  cout << "Exc Particle 3 = " << fExcitation3 << " MeV" << endl;
  cout << "Exc Particle 4 = " << fExcitation4 << " MeV" << endl;
  cout << "Qgg = " << fQValue << " MeV" << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Reaction::ReadConfigurationFile(string Path){
  ifstream ReactionFile;
  string GlobalPath = getenv("NPTOOL");
  string StandardPath = GlobalPath + "/Inputs/EventGenerator/" + Path;
  ReactionFile.open(Path.c_str());
  if (!ReactionFile.is_open()) {
    ReactionFile.open(StandardPath.c_str());
    if(ReactionFile.is_open()) {
      Path = StandardPath;
    }
    else {cout << "Reaction File " << Path << " not found" << endl;exit(1);}
  }
  NPL::InputParser parser(Path);
  ReadConfigurationFile(parser);
}
////////////////////////////////////////////////////////////////////////////////
Particle Reaction::GetParticle(string name, NPL::InputParser parser){
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
    NPL::SendErrorAndExit("NPL::Reaction","Too many nuclei define with the same name");
  }

  return (NPL::Particle());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void Reaction::ReadConfigurationFile(NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("TwoBodyReaction");
  if(blocks.size()>0 && NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// Two body reaction found " << endl;

  vector<string> token1 = {"Beam","Target","Light","Heavy"};
  vector<string> token2 = {"Beam","Target","Particle3","Particle4"};
  double CSHalfOpenAngleMin= 0*deg ;
  double CSHalfOpenAngleMax= 180*deg;
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token1)){
      int v = NPOptionManager::getInstance()->GetVerboseLevel();
      NPOptionManager::getInstance()->SetVerboseLevel(0);
      fParticle1.ReadConfigurationFile(parser);
      NPOptionManager::getInstance()->SetVerboseLevel(v);

      fBeamEnergy= fParticle1.GetEnergy();
      fParticle2 = GetParticle(blocks[i]->GetString("Target"),parser);
      fParticle3 = GetParticle(blocks[i]->GetString("Light"),parser);
      fParticle4 = GetParticle(blocks[i]->GetString("Heavy"),parser);
    }
    else if(blocks[i]->HasTokenList(token2)){
      fParticle1.SetVerboseLevel(0);
      fParticle1.ReadConfigurationFile(parser);
      fBeamEnergy= fParticle1.GetEnergy();

      fParticle2 = GetParticle(blocks[i]->GetString("Target"),parser);
      fParticle3 = GetParticle(blocks[i]->GetString("Particle3"),parser);
      fParticle4 = GetParticle(blocks[i]->GetString("Particle4"),parser);
    }
    else{
      cout << "ERROR: check your input file formatting \033[0m" << endl;
      exit(1);
    }

    if(blocks[i]->HasToken("ExcitationEnergyBeam")){
      fExcitation1 = blocks[i]->GetDouble("ExcitationEnergyBeam","MeV");
    }
    else if(blocks[i]->HasToken("ExcitationEnergy1")){
      fExcitation1 = blocks[i]->GetDouble("ExcitationEnergy1","MeV");
    }

    if(blocks[i]->HasToken("ExcitationEnergyLight"))
      fExcitation3 = blocks[i]->GetDouble("ExcitationEnergyLight","MeV");
    else if(blocks[i]->HasToken("ExcitationEnergy3"))
      fExcitation3 = blocks[i]->GetDouble("ExcitationEnergy3","MeV");

    if(blocks[i]->HasToken("ExcitationEnergyHeavy"))
      fExcitation4 = blocks[i]->GetDouble("ExcitationEnergyHeavy","MeV");
    else if(blocks[i]->HasToken("ExcitationEnergy4"))
      fExcitation4 = blocks[i]->GetDouble("ExcitationEnergy4","MeV");

    if(blocks[i]->HasToken("ExcitationEnergyDistribution")){
      vector<string> file = blocks[i]->GetVectorString("ExcitationEnergyDistribution");
      fExcitationEnergyHist = Read1DProfile(file[0],file[1]);
      fExcitation4 = 0 ;
    }

    if(blocks[i]->HasToken("CrossSectionPath")){
      vector<string> file = blocks[i]->GetVectorString("CrossSectionPath");
      TH1F* CStemp = Read1DProfile(file[0], file[1]);

      // multiply CStemp by sin(theta)
      TF1* fsin = new TF1("sin",Form("1/(sin(x*%f/180.))",M_PI),0,180);
      CStemp->Divide(fsin,1);
      SetCrossSectionHist(CStemp);
      delete fsin;
    }

    if(blocks[i]->HasToken("LabCrossSectionPath")){
      fLabCrossSection=true;

      vector<string> file = blocks[i]->GetVectorString("LabCrossSectionPath");
      TH1F* CStemp = Read1DProfile(file[0], file[1]);

      // multiply CStemp by sin(theta)
      TF1* fsin = new TF1("sin",Form("1/(sin(x*%f/180.))",M_PI),0,180);
      CStemp->Divide(fsin,1);
      SetCrossSectionHist(CStemp);
      delete fsin;
    }


    if(blocks[i]->HasToken("DoubleDifferentialCrossSectionPath")){
      vector<string> file = blocks[i]->GetVectorString("DoubleDifferentialCrossSectionPath");
      TH2F* CStemp = Read2DProfile(file[0],file[1]);

      // multiply CStemp by sin(theta)
      // X axis is theta CM
      // Y axis is beam energy
      // Division affect only X axis
      TF1* fsin = new TF1("sin",Form("1/(sin(x*%f/180.))",M_PI),0,180);
      CStemp->Divide(fsin,1);

      SetDoubleDifferentialCrossSectionHist(CStemp);
      delete fsin;
    }

    if(blocks[i]->HasToken("HalfOpenAngleMin")){
      CSHalfOpenAngleMin = blocks[i]->GetDouble("HalfOpenAngleMin","deg");
    }
    if(blocks[i]->HasToken("HalfOpenAngleMax")){
      CSHalfOpenAngleMax = blocks[i]->GetDouble("HalfOpenAngleMax","deg");
    }
    if(blocks[i]->HasToken("Shoot3")){
      fshoot3 = blocks[i]->GetInt("Shoot3");
    }
    if(blocks[i]->HasToken("Shoot4")){
      fshoot4 = blocks[i]->GetInt("Shoot4");
    }
    if(blocks[i]->HasToken("ShootHeavy")){
      fshoot4 = blocks[i]->GetInt("ShootHeavy");
    }
    if(blocks[i]->HasToken("ShootLight")){
      fshoot3 = blocks[i]->GetInt("ShootLight");
    }
    if(blocks[i]->HasToken("UseExInGeant4")){
      // This option will not change the Ex of the produced ion in G4 Tracking
      // This is to be set to true when using a Ex distribution without decay
      // Otherwise the Ion Table size grew four ech event slowing down the
      // simulation
      fUseExInGeant4 = blocks[i]->GetInt("UseExInGeant4");
    }

  }
  SetCSAngle(CSHalfOpenAngleMin/deg,CSHalfOpenAngleMax/deg);
  initializePrecomputeVariable();
  cout << "\033[0m" ;
}

////////////////////////////////////////////////////////////////////////////////////////////
void Reaction::initializePrecomputeVariable(){

  if(fBeamEnergy < 0)
    fBeamEnergy = 0 ;

  if(fExcitation1>=0) fParticle1.SetExcitationEnergy(fExcitation1); // Write over the beam excitation energy

  //fParticle1.GetExcitationEnergy() is a copy of fExcitation1
  m1 = fParticle1.Mass() + fParticle1.GetExcitationEnergy();// in case of isomeric state, 
  m2 = fParticle2.Mass(); // Target
  m3 = fParticle3.Mass() + fExcitation3;
  m4 = fParticle4.Mass() + fExcitation4;
  fQValue =m1+m2-m3-m4;

  s = m1*m1 + m2*m2 + 2*m2*(fBeamEnergy + m1);
  fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));
  fEcm = sqrt(s) - m1 -m2;

  ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
  ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));
  ECM_3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
  ECM_4 = (s + m4*m4 - m3*m3)/(2*sqrt(s));

  pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
  pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);
  pCM_3 = sqrt(ECM_3*ECM_3 - m3*m3);
  pCM_4 = sqrt(ECM_4*ECM_4 - m4*m4);

  fImpulsionLab_1 = TVector3(0,0,sqrt(fBeamEnergy*fBeamEnergy + 2*fBeamEnergy*m1));
  fImpulsionLab_2 = TVector3(0,0,0);

  fEnergyImpulsionLab_1= TLorentzVector(fImpulsionLab_1,m1+fBeamEnergy);
  fEnergyImpulsionLab_2= TLorentzVector(fImpulsionLab_2,m2);

  fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;

  BetaCM = fTotalEnergyImpulsionLab.Beta();

  fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
  fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);

  fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
  fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);

}

////////////////////////////////////////////////////////////////////////////////////////////
void Reaction::SetParticle3(double EnergyLab, double ThetaLab){
  double p3 = sqrt(pow(EnergyLab,2) + 2*m3*EnergyLab);

  fEnergyImpulsionLab_3 = TLorentzVector(p3*sin(ThetaLab),0,p3*cos(ThetaLab),EnergyLab+m3);
  fEnergyImpulsionLab_4 = fTotalEnergyImpulsionLab - fEnergyImpulsionLab_3;

  fParticle3.SetEnergyImpulsion(fEnergyImpulsionLab_3);
  fParticle4.SetEnergyImpulsion(fEnergyImpulsionLab_4);

  fThetaCM = EnergyLabToThetaCM(EnergyLab, ThetaLab);
  fExcitation4 = ReconstructRelativistic(EnergyLab, ThetaLab);
}

////////////////////////////////////////////////////////////////////////////////////////////
TGraph* Reaction::GetKinematicLine3(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta3,E3,theta4,E4;

  for (double angle=0 ; angle < 360 ; angle+=AngleStep_CM){
    SetThetaCM(angle*deg);
    KineRelativistic(theta3, E3, theta4, E4);
    fParticle3.SetKineticEnergy(E3);

    if(E3>0){
      vx.push_back(theta3/deg);
      vy.push_back(E3);
    }
  }
  fKineLine3 = new TGraph(vx.size(),&vx[0],&vy[0]);

  return(fKineLine3);
}

////////////////////////////////////////////////////////////////////////////////////////////
TGraph* Reaction::GetKinematicLine4(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta3,E3,theta4,E4;

  for (double angle=0 ; angle < 360 ; angle+=AngleStep_CM){
    SetThetaCM(angle*deg);
    KineRelativistic(theta3, E3, theta4, E4);
    fParticle4.SetKineticEnergy(E4);
    if(E4>0){
      vx.push_back(theta4/deg);
      vy.push_back(E4);
    }
  }
  fKineLine4= new TGraph(vx.size(),&vx[0],&vy[0]);

  return(fKineLine4);
}

////////////////////////////////////////////////////////////////////////////////////////////
TGraph* Reaction::GetTheta3VsTheta4(double AngleStep_CM)
{

  vector<double> vx;
  vector<double> vy;
  double theta3,E3,theta4,E4;

  for (double angle=0 ; angle < 360 ; angle+=AngleStep_CM){
    SetThetaCM(angle*deg);
    KineRelativistic(theta3, E3, theta4, E4);

    vx.push_back(theta3/deg);
    vy.push_back(theta4/deg);
  }
  fTheta3VsTheta4= new TGraph(vx.size(),&vx[0],&vy[0]);
  return(fTheta3VsTheta4);
}

////////////////////////////////////////////////////////////////////////////////////////////
TGraph* Reaction::GetBrhoLine3(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta3,E3,theta4,E4;
  double Brho;

  for (double angle=0 ; angle < 360 ; angle+=AngleStep_CM){
    SetThetaCM(angle*deg);
    KineRelativistic(theta3, E3, theta4, E4);
    fParticle3.SetKineticEnergy(E3);
    Brho = fParticle3.GetBrho();

    vx.push_back(theta3/deg);
    vy.push_back(Brho);
  }
  fLineBrho3= new TGraph(vx.size(),&vx[0],&vy[0]);
  return(fLineBrho3);
}

////////////////////////////////////////////////////////////////////////////////////////////
TGraph* Reaction::GetThetaLabVersusThetaCM(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta3,E3,theta4,E4;

  for (double angle=0 ; angle < 360 ; angle+=AngleStep_CM){
    SetThetaCM(angle*deg);
    KineRelativistic(theta3, E3, theta4, E4);

    vx.push_back(fThetaCM/deg);
    vy.push_back(theta3/deg);
  }

  fAngleLine= new TGraph(vx.size(),&vx[0],&vy[0]);
  return(fAngleLine);
}
////////////////////////////////////////////////////////////////////////////////////////////
TGraph* Reaction::GetELabVersusThetaCM(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta3,E3,theta4,E4;

  for (double angle=0 ; angle < 360 ; angle+=AngleStep_CM){
    SetThetaCM(angle*deg);
    KineRelativistic(theta3, E3, theta4, E4);

    vx.push_back(E3);
    vy.push_back(fThetaCM/deg);
  }

  fAngleLine= new TGraph(vx.size(),&vx[0],&vy[0]);
  return(fAngleLine);
}


////////////////////////////////////////////////////////////////////////////////////////////
Double_t Reaction::GetTotalCrossSection() const {
  Double_t stot = fCrossSectionHist->Integral("width"); // take bin width into account (in deg!)
  stot *= M_PI/180; // correct so that bin width is in rad
  stot *= 2*M_PI;   // integration over phi

  return stot;
}



////////////////////////////////////////////////////////////////////////////////////////////
void Reaction::PrintKinematic(){
  int size = 360;
  double theta3,E3,theta4,E4,Brho3,Brho4;

  cout << endl;
  cout << "*********************** Print Kinematic ***********************" << endl;
  cout << "ThetaCM" << "	" << "ThetaLab" << " " << "EnergyLab3" << "	" << "Brho3" << "	" << "EnergyLab4" << "	" << "Brho4" << endl;
  for (int i = 0; i < size; ++i){
    SetThetaCM(((double)i)/2*deg);
    KineRelativistic(theta3, E3, theta4, E4);

    fParticle3.SetKineticEnergy(E3);
    Brho3 = fParticle3.GetBrho();

    fParticle4.SetKineticEnergy(E4);
    Brho4 = fParticle4.GetBrho();

    cout << (double)i/2 << "	" << theta3/deg << "	" << E3 << "	" << Brho3 << "		" << E4 << "	" << Brho4 << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////
void Reaction::SetCSAngle(double CSHalfOpenAngleMin,double CSHalfOpenAngleMax){
  if(fCrossSectionHist){
    for (int i = 0 ; i< fCrossSectionHist->GetNbinsX(); i++){
      if( fCrossSectionHist->GetBinCenter(i) > CSHalfOpenAngleMax || fCrossSectionHist->GetBinCenter(i) < CSHalfOpenAngleMin){
        fCrossSectionHist->SetBinContent(i,0);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Check whenever the reaction is allowed at the given energy
bool Reaction::IsAllowed(double Energy){
  double AvailableEnergy = Energy + fParticle1.Mass() + fParticle2.Mass();
  double RequiredEnergy  = fParticle3.Mass() + fParticle4.Mass();

  if(AvailableEnergy>RequiredEnergy)
    return true;
  else
    return false;
}
