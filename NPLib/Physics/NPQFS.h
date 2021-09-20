#ifndef NPQFS_h
#define NPQFS_h
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
 *  First step (dissociation):  A -> B + c                                   *
 *  Second step (scattering) :  c + T -> 1 + 2                               *
 *  Labeling is:                                                             *
 *                                                                           *
 *              A --> T  ==> B + (c -> T) =>  B + 1 + 2                      *
 *                                                                           *
 *  where:                                                                   *
 *    +  A is the beam nucleus                                               *
 *    +  T is the target nucleon (proton)                                    *
 *                                                                           *
 *    +  B is the residual fragment (beam-like)                              *
 *    +  1 is the scattered target nucleon  (former T)                       *
 *    +  2 is the knocked-out cluster/nucleon (noted c) in the intermediate  *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    +  Adapted from original event generator from V. Panin (R3B collab)    *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>

// NPL
#include "NPParticle.h"
#include "NPBeam.h"
#include "NPInputParser.h"
using namespace NPL;

// ROOT header
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;

namespace NPL{

  class QFS{

    public:  // Constructors and Destructors
      QFS();
      ~QFS();

    private:
    int fVerboseLevel;
    Beam     fParticleA;                 // Beam (A)
    Particle  fParticleT;                 // Target (T)
    Particle  fParticleB;                 // Beam-like ejectile (B)
    Particle  fParticle1;                 // Target-like ejectile (1)
    Particle  fParticle2;                 // Knocked-out nucleon/cluster (2)
    double   fQValue;                  // Q-value in MeV
    double   fEcm;                     // Ecm in MeV
    double   fThetaCM;                 // Center-of-mass theta angle in radian
    double   fPhiCM;                   // Center-of-mass Phi angle in radian
    double   fBeamEnergy;              // Beam energy in MeV
    double   fExcitationA;             // Excitation energy in MeV of the beam, useful for isomers 
    double   fExcitationB;             // Excitation energy in MeV of beam-like ejectile
    double   fMomentumSigma;           // Width of the momentum distribution  (sigma)             
    TVector3 fInternalMomentum;        // Internal momentum of the removed cluster            
    TH1F*    fPerpMomentumHist;        // Perpendicular momentum distribution in beam at rest frame
    TH1F*    fParMomentumHist;         // Parallel momentum distribution in beam at rest frame
    double   fisotropic;               

    TGraph* fTheta2VsTheta1;
    TGraph* fPhi2VsPhi1;

    // used for MC simulations
    bool fshootB; // shoot beam-like ejectile
    bool fshoot1; // shoot light ejectile &
    bool fshoot2; // shoot light ejectile 2
    bool fUseExInGeant4; 
    bool fDeexcitation;
    
    public:
    Particle GetParticle(string name, NPL::InputParser parser);
    void ReadConfigurationFile(string Path);
    void ReadConfigurationFile(NPL::InputParser);
    void CalculateVariables();
    void TestR3B();
    void KineRelativistic(double &ThetaLab1, double &PhiLab1, double &KineticEnergyLab1, double &ThetaLab2, double &PhiLab2, double &KineticEnergyLab2);
    TVector3 ShootInternalMomentum();
    bool IsAllowed();
    void Dump();

    private: // intern precompute variable
    double mA;
    double mB;
    double ma_off;
    double ma;
    double mT;
    double m1;
    double m2;

    TVector3 Pa, PB, Pa_lab, PB_lab;
    double Ea_lab,EB_lab;

    // Lorentz Vector
    TLorentzVector fEnergyImpulsionLab_A;
    TLorentzVector fEnergyImpulsionLab_B;
    TLorentzVector fEnergyImpulsionLab_a;
    TLorentzVector fEnergyImpulsionLab_T;
    TLorentzVector fEnergyImpulsionLab_1;
    TLorentzVector fEnergyImpulsionLab_2;
    TLorentzVector fTotalEnergyImpulsionLab;

    TLorentzVector fEnergyImpulsionCM_A;
    TLorentzVector fEnergyImpulsionCM_B;
    TLorentzVector fEnergyImpulsionCM_a;
    TLorentzVector fEnergyImpulsionCM_T;
    TLorentzVector fEnergyImpulsionCM_1;
    TLorentzVector fEnergyImpulsionCM_2;
    TLorentzVector fTotalEnergyImpulsionCM;

    // Impulsion Vector3
    TVector3 fImpulsionLab_a;
    TVector3 fImpulsionLab_T;
    TVector3 fImpulsionLab_1;
    TVector3 fImpulsionLab_2;

    TVector3 fImpulsionCM_a;
    TVector3 fImpulsionCM_T;
    TVector3 fImpulsionCM_1;
    TVector3 fImpulsionCM_2;

    // CM Energy composante & CM impulsion norme
    Double_t ECM_a;
    Double_t ECM_T;
    Double_t ECM_1;
    Double_t ECM_2;
    Double_t pCM_a;
    Double_t pCM_T;
    Double_t pCM_1;
    Double_t pCM_2;

    // Mandelstam variable
    Double_t s;

    // Center of Mass Kinematic
    Double_t BetaCM;

    public:
    //SETTERS
    void SetBeamEnergy(const double& eBeam) {fBeamEnergy = eBeam;}
    void SetThetaCM(const double& angle) {fThetaCM = angle;}
    void SetPhiCM(const double& angle) {fPhiCM = angle;}
    void SetInternalMomentum(const TVector3& mom) {fInternalMomentum = mom;}
    void SetMomentumSigma(const double& sigma) {fMomentumSigma = sigma;}
    void SetPerpMomentumHist  (TH1F*  PerpMomentumHist)
        {delete fPerpMomentumHist; fPerpMomentumHist   = PerpMomentumHist;}
    void SetParMomentumHist  (TH1F*  ParMomentumHist)
        {delete fParMomentumHist; fParMomentumHist   = ParMomentumHist;}

    //GETTERS
    Particle*  GetParticleA()               {return &fParticleA;}
    Particle*  GetParticleT()               {return &fParticleT;}
    Particle*  GetParticleB()               {return &fParticleB;}
    Particle*  GetParticle1()               {return &fParticle1;}
    Particle*  GetParticle2()               {return &fParticle2;}
    bool     GetShoot1()         const        {return fshoot1;}
    bool     GetShoot2()         const        {return fshoot2;}
    bool     GetShootB()         const        {return fshootB;}
    bool     GetDeexcitation()   const        {return fDeexcitation;}
    double   GetThetaCM()        const        {return fThetaCM;}
    double   GetPhiCM()          const        {return fPhiCM;}
    double   GetMomentumSigma()  const        {return fMomentumSigma;}
    bool      GetUseExInGeant4() const { return fUseExInGeant4; }
    double   GetExcitationA() const           {return fExcitationA;}
    double   GetExcitationB() const           {return fExcitationB;}
    TVector3 GetInternalMomentum() const   {return fInternalMomentum;}
 
    TLorentzVector*  GetEnergyImpulsionLab_A() {return &fEnergyImpulsionLab_A;}
    TLorentzVector*  GetEnergyImpulsionLab_T() {return &fEnergyImpulsionLab_T;}
    TLorentzVector*  GetEnergyImpulsionLab_a() {return &fEnergyImpulsionLab_a;}
    TLorentzVector*  GetEnergyImpulsionLab_1() {return &fEnergyImpulsionLab_1;}
    TLorentzVector*  GetEnergyImpulsionLab_2() {return &fEnergyImpulsionLab_2;}
    TLorentzVector*  GetEnergyImpulsionLab_B() {return &fEnergyImpulsionLab_B;}

    TLorentzVector*  GetEnergyImpulsionCM_A() {return &fEnergyImpulsionCM_A;}
    TLorentzVector*  GetEnergyImpulsionCM_T() {return &fEnergyImpulsionCM_T;}
    TLorentzVector*  GetEnergyImpulsionCM_a() {return &fEnergyImpulsionCM_a;}
    TLorentzVector*  GetEnergyImpulsionCM_1() {return &fEnergyImpulsionCM_1;}
    TLorentzVector*  GetEnergyImpulsionCM_2() {return &fEnergyImpulsionCM_2;}
    TLorentzVector*  GetEnergyImpulsionCM_B() {return &fEnergyImpulsionCM_B;}
 
    TGraph* GetTheta2VsTheta1(double AngleStep_CM=1);
    TGraph* GetPhi2VsPhi1(double AngleStep_CM=1);



    
    //TO REMOVE AT SOME POINT WHEN CLASS IS ROBUSTLY TESTED
    private:
    // R3B Methods and Variables used as a starting point for this class (useful for checks)
    void CalculateVariablesOld();
    pair<double,double> Lorentz(double, double, double, double);
    double function(double, double, double);
    void KineR3B(double, double, double, double, double);
    TVector3 Rotations(TVector3,TVector3);
    double e_clust;
    double p_clust;
    double theta_clust;
    double e_scat;
    double p_scat;
    double theta_scat;
    bool good;
    double T;
    

    ClassDef(QFS,0)
  };
}
#endif
