#ifndef GEF_CLASS
#define GEF_CLASS
/***********************************************************************
*                                                                      *
* Original Author    : Diego Ramos      -> diego.ramos@ganil.fr        *
* Adapted for NPTool : Pierre Morfouace -> pierre.morfouace2@cea.fr    *
* Creation Date      : 25/09/2020                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Description:                                                         *
* This class generates fission fragments based on a simplified         *
* model of the GEF code.                                               *
* This is not the full GEF code. User should use it with caution.      *
*                                                                      *
*----------------------------------------------------------------------*
* Comment:                                                             *
*                                                                      *
*                                                                      *
************************************************************************/
// C++ header
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <exception>
#include <vector>
#include <list>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <string>
#include <ctime>
#include <fcntl.h>
#include <unistd.h> 
#include <sys/stat.h>
#include <time.h>
#include <typeinfo>

// ROOT header
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TKey.h>
#include <TRandom3.h>

// NPTool header
#include "NPParticle.h"

using namespace std;

class GEF
{
public:
  GEF(NPL::Particle FissNucl);
  ~GEF(void);
  
  inline double GetBrhoffl(void) {return Brhoffl;}
  inline double GetBrhoffh(void) {return Brhoffh;}
  inline double GetThffl(void)   {return Thlab_light;}
  inline double GetThffh(void)   {return Thlab_heavy;}
  inline double GetPhffl(void)   {return Phlab_light;}
  inline double GetPhffh(void)   {return Phlab_heavy;}
  inline float  GetKEffl(void)   {return KElab_light;}
  inline float  GetKEffh(void)   {return KElab_heavy;}
  inline float  GetVffl(void)    {return vlab_light;}
  inline float  GetVffh(void)    {return vlab_heavy;}
  inline int    GetZffl(void)    {return I_Z_light_sci;}
  inline int    GetZffh(void)    {return I_Z_heavy_sci;}
  inline int    GetAffl(void)    {return I_A_light_sci;}
  inline int    GetAffh(void)    {return I_A_heavy_sci;}

  float* GetNeutronEnergyFrag1() {return Array_E_n1_frag1;}
  float* GetNeutronEnergyFrag2() {return Array_E_n2_frag2;}
  float* GetGammaEnergyFrag1()   {return Array_Eg0_light;} 
  float* GetGammaEnergyFrag2()   {return Array_Eg0_heavy;} 
  float GetTKE() {return TKE;}
  float GetKE1() {return Ekinlight_sci;}
  float GetKE2() {return Ekinheavy_sci;}

  void Treat(void);
  void InitCompound(double vEx, double vEfis, double vLfis=0, double vThfis=0, double vPhfis=0);
  bool IsValid(int Z, int A);

private:
  NPL::Particle FIS;
  Int_t Afis;
  Int_t Zfis;
  Float_t Ex;
  Float_t Lfis; //hbar
  Float_t Efis;
  Float_t Thfis;
  Float_t Phfis;

  Float_t BEldmTF[204][137];  //[N][Z]
  void ReadBEldmTF(void);
  Float_t BEexp[204][137];    //[N][Z]
  void ReadBEexp(void);
  Float_t DEFOtab[237][137];  //[N][Z]
  void ReadDEFOtab(void);
  Float_t ShellMO[204][137];  //[N][Z]
  void ReadShellMO(void);
  Float_t R_SPI[112][284];  //[Z][A] //Spin ground-states
  void ReadNucProp(void);
  
  void ReadGEFParameters(void);
  void Init(void);
  
  float Max(float a, float b);
  float Min(float a, float b);
  float Erf(float x);
  float Erfc(float x);
  float Tanh(float x);
  float Coth(float x);
  float Log10(float R);
  float Bell(float xpos, float xleft, float xright);
  float Csng(int a);
  int CInt(float a);
  float Mod(float a, float b);
  float FLOOR(float a);
  float CEIL(float a);
  int Sgn(float a);

  void FissionGENERAL(void);
  void CentralZ_FM(void);
  void MeanDefScission(void);
  void MeanZ_asA(void);
  void RelationsZ_A_FM(void);
  void PotCurv_FM(void);
  void EnergyTrans(void);
  void FissionBarriers(void);
  void BarriersEx_FM(void);
  void TCollective(void);
  void PolarizationStiffness(void);
  void MeanValues_FM(void);
  void EnergyDependence(void);
  void Yields_FM(void);
  void MassWidths_FM(void);
  void ShellEff_FM(void);
  void IntrinsicT(void);
  void IntrinsicEx(void);
  void FF_angularMomentum(void);
  
  void FissionModeChoise(void);
  void ZAChoise(void);
  void N_Saddle_Scission(void);
  void FragmentsEnergy(void);
  void FromCMtoLAB(void);
  void FissionProbability(void);
  
  
  float Getyield(float E_rel,float E_ref, float T_low, float T_high);
  float F1(float Z_S_A);
  float F2(float Z_S_A);
  float Masscurv(float Z, float A, float RL, float kappa);
  float Masscurv1(float Z,float A, float RL, float kappa);
  float De_Saddle_Scission(float Z_square_over_Athird,float ESHIFTSASCI);
  float TEgidy(float A,float DU, float Fred);
  float TRusanov(float E, float A);
  float LyMass(float Z, float A, float beta);
  float LyPair(int Z, int A);
  float TFPair(int Z, int A);
  float Pmass(float Z, float A, float beta);
  float FEDEFOLys(float Z, float A, float beta);
  float FEDEFOP(float Z, float A, float beta);
  float LDMass(float Z, float A, float beta);
  float AME2012(int IZ, int IA);
  float U_SHELL(int Z, int A);
  float U_SHELL_exp(int IZ, int IA);   
  float U_SHELL_EO_exp(int IZ, int IA);
  float U_MASS(float Z, float A);
  float ECOUL(float Z1, float A1, float beta1, float Z2, float A2, float beta2, float d);
  float beta_light(int Z, float betaL0, float betaL1);
  float beta_heavy(int Z, float betaH0, float betaH1);
  float Z_equi(int ZCN, int A1, int A2, float beta1, float beta2, float d, int Imode);
  void Beta_opt_light(float A1, float A2, float Z1, float Z2, float d, float beta2_imposed, float &beta1_opt);	     
  void Beta_Equi(float A1, float A2, float Z1, float Z2, float d, float beta1prev, float beta2prev, float &beta1opt, float &beta2opt);
  void Eva(int Ilh, float Z_CN, float A_CN, float E_INIT, float T, float J_frag, float &Z_RES, float &A_RES, float &E_FINAL, float *Array_En, float *Array_Tn, float *Array_Eg0); // it is not well translated from Visual BASIC  need to be checked
   float u_accel(float A1, float Z1, float A2, float Z2, float TKE , float E0, float Tn); 
  float P_Egamma_high(float Zi, float Ai, float Ei);
  float P_Egamma_low(float Zi, float Ai, float Ei);     
  float U_Ired(float Z, float A);
  float U_IredFF(float Z, float A);
  float U_I_Shell(float Z, float A);
  float U_alev_ld(float Z, float A);
  float U_Temp(float Z, float A, float E, int Ishell, int Ipair, float Tscale, float Econd);
  double U_levdens_Egidy(int Z, int A, float E, int Ishell, int Ipair, float Tscale, float Econd, float af_an);
  double U_levdens_FG(int Z, int A, float E, int Ishell, int Ipair, float Tscale, float Econd, float af_an);
  double U_levdens(int Z, int A, float E, int Ishell, float Ipair, float Tscale, float Econd, float af_an);
  float U_Temp2(float Z, float A, float E, float Rshell, float  Rpair, float Tscale, float Econd);
  float E0_GDR(float Z, float A);
  float Width_GDR(float E0);
   //float Sigma_GDR(float Z, float A, float E, float E0, float WidthK);//function not defined in GEF
   float Efac_def_GDR(float Beta, float Gamma, float K);  
  // float TK_GDR(float Z, float A, float E);   //function not defined in GEF
   float GgGtot(float Z, float A, float E, float Egamma);
   float E_next(float T1, float T2, float E1, float E2, float A1, float A2);
   int EVEN_ODD(float R_ORIGIN,float R_EVEN_ODD);
   float U_Even_Odd(int I_Channel, float PEO);
   float BFTF(float RZ, float RA, int I_Switch);
   float BFTFA(float RZ, float RA, int I_Switch);
   float BFTFB(float RZ, float RA, int I_Switch);
   float Gaussintegral(float R_x, float R_sigma);

  float U_Box(float x, float sigma, float width);
  float U_Box2(float x, float sigma1, float sigma2, float width);
  float U_Gauss(float x, float sigma);
  float U_Gauss_abs(float x, float sigma);
  float U_Gauss_mod(float x,float sigma);
  float U_LinGauss(float x, float R_Sigma);
  float Round(float R, int N);        
  float Pexplim(float R_lambda, float xmin, float xmax);
  float PBox(float Mean, float Sigma, float Bottom); 
  float PBox2(float Mean, float Sigma1, float Sigma2,  float Bottom);
  float PPower(int Order, float Rmin, float Rmax);    
  float PPower_Griffin_v(int Order, float Rmin, float Rmax);    
  float PPower_Griffin_E(int Order, float Rmin,float  Rmax);    
  float PGauss(float Mean, float Sigma);
  float PLinGauss(float R_Sigma);
  float PExp(float R_Tau);
  float PMaxwell(float R_T);
  float PMaxwellMod(float R_T, float R_A);
  float PMaxwellv(float R_T);
  long int Modulo(unsigned long int I, unsigned long int J);
  //int PLoss(unsigned long int IL); //I don't know yet how to translate this function to C++
  bool U_Valid(int I_Z, int I_A);
  float U_Delta_S0(int I_Z, int I_A); 


  float pi;
  TRandom3 *rn;
  
  //Internal variables
  
  int I_N_CN; // Neutron number of fissioning nucleus  
  float T_Coll_Mode_1,T_Coll_Mode_2,T_Coll_Mode_3,T_Coll_Mode_4;
  float T_asym_Mode_1,T_asym_Mode_2,T_asym_Mode_3,T_asym_Mode_4,T_asym_Mode_0;
  float Sigpol_Mode_1,Sigpol_Mode_2,Sigpol_Mode_3,Sigpol_Mode_4;
  float R_Z_Curv_S0,R_Z_Curv1_S0,R_A_Curv1_S0;
  float ZC_Mode_0,ZC_Mode_1,ZC_Mode_2,ZC_Mode_3,ZC_Mode_4;
  float ZC_Mode_3_shift;
  float SigZ_Mode_0,SigZ_Mode_1,SigZ_Mode_2,SigZ_Mode_3,SigZ_Mode_4;
  float SigZ_SL4;
  float SN,Sprot;
  float E_exc_S0_prov,E_exc_S1_prov,E_exc_S2_prov,E_exc_S3_prov,E_exc_S4_prov;
  float E_exc_S11_prov,E_exc_S22_prov;
  float E_exc_Barr;
  float E_LD_S1,E_LD_S2,E_LD_S3,E_LD_S4;
  float R_Shell_S1_eff,R_Shell_S2_eff,R_Shell_S3_eff,R_Shell_S4_eff;
  float Yield_Norm;
  float R_E_exc_eff;
  float R_Z_Heavy,R_Z_Light;
  int I_Mode;
  float T_Pol_Mode_0,T_Pol_Mode_1,T_Pol_Mode_2,T_Pol_Mode_3,T_Pol_Mode_4;
  float E_Min_Barr;
  float RI;
  float rbeta, beta1, beta2;
  float rbeta_ld, rbeta_shell;
  float ZUCD;
  float Z;
  float E_tunn;
  float beta1_opt,beta2_opt,beta1_prev,beta2_prev;
  float Z1,Z2;
  int IZ1,IN1,IZ2,IN2;
  float A1,A2;
  int IA1,IA2;
  float E_defo;
  float R_Pol_Curv_S0, R_Pol_Curv_S1, R_Pol_Curv_S2,R_Pol_Curv_S3,R_Pol_Curv_S4;
  float RA,RZ;
  float SigA_Mode_0, SigA_Mode_1, SigA_Mode_2,SigA_Mode_3,SigA_Mode_4;
  float AC_Mode_0, AC_Mode_1, AC_Mode_2, AC_Mode_3, AC_Mode_4;
  float R_A_heavy, R_A_light;
  float RZpol;
      
  float Eexc_light,Eexc_heavy;
 
  float T_intr_Mode_0,T_intr_Mode_1_heavy,T_intr_Mode_1_light;
  float T_intr_Mode_2_heavy,T_intr_Mode_2_light;
  float T_intr_Mode_3_heavy,T_intr_Mode_3_light;
  float T_intr_Mode_4_heavy,T_intr_Mode_4_light;
  float T;
  float DU0,DU1,DU2,DU3,DU4;
    
  float E_intr,E_intr_light,E_intr_heavy;
  float E_intr_light_S0_mac,E_intr_heavy_S0_mac,RW_mac;
  float E_intr_light_S0_mic,E_intr_heavy_S0_mic;
  float E_intr_light_mean,E_intr_heavy_mean;
  float Eexc_heavy_mean, Eexc_light_mean;
  float Ecoll,Ecoll_mean,Ecoll_heavy,Ecoll_light;
  const static float E_EXC_TRUE,E_EXC_ISO;  // Excitation energy of isomeric state
  float Qvalue,TKEmin;
  float Rtest;
  float Theavy,Tlight;
    
  float T_low_S1_used;
  float SigA_Mode_11,SigA_Mode_22;
  int Ngtot;
  int Nglight;
  int Ngheavy;
  float Egtot1000;
  float S1_enhance, S2_enhance, S1_enhance_S2;
  float DZ_S2_lowE;    
  int I_A_CN,I_Z_CN;
  float R_E_exc_used;
  
  float TKE, TKE_post, Ekinlight, Ekinheavy, Ekinlight_post, Ekinheavy_post;
  float Ekinlight_sci, Ekinheavy_sci;
  float vkinlight,vkinheavy;

  //  float Beta(-1 To 6,1 To 2,150);
      // -1: microscopic; 0: macroscopic for S0 fission channel
  float Beta[8][3][151]; //******remember to add Beta[a+1][b][c] in the translation***********
  float Edefo[6][3][151];//******remember to add Edefo[a+1][b][c] in the translation***********
  float Zmean[5][3][351];
  float Zshift[5][3][351];
  float ZshiftOriginal[5][3][351];
  float Temp[5][3][351];
  float TempFF[5][3][351];
  float EShell[5][3][351];
  float PEOZ[7][3][351];
  float PEON[7][3][351];   //pre-neutron evaporation
  float EPART[7][3][351];
  float SpinRMSNZ[7][3][201][151];
  
  // Input parameters: 

  //Dim As String kin   // Key input
    //Dim As String Cyesno
  const static int I_thread=0;  // Thread number of this process   
      // Dim Shared As UByte B_Error_On = 0        // Error analysis required  
      //Dim As UByte B_Error_Analysis
  int N_Error_Max;    // Number of different random parameter sets      
  int I_Error;         // Counts the parameter sets  
  // I_Error runs from 0 to N_Error_Max - 1!
   //Dim As UByte B_Random_on = 0       // Write ENDF random files  
  int I_DelGam;
  int P_Z_CN;           //`Z of fissioning nucleus  
  int P_A_CN;          // A of fissioning nucleus  
  float P_E_exc;                   // Energy above lowest outer barrier EB  
    
  //Dim As UByte B_Double_Covar = 0    // Option: covariances for yields of 2 fragments  
  int I_Double_Covar;      // Sequence of calculations for 2 fragments  
  int N_Double_Covar;
  int P_Z_CN_Double; 
  int P_A_CN_Double;
  float P_E_exc_Double;
    
  float P_I_rms_CN;                  // rms initial angular momentum  


 //Model parameters of GEF

  float  Emax_valid;      // Maximum allowed excitation energy  
  float Eexc_min_multi;            // Threshold for calc. of multi-chance fission  
  float _Delta_S0;         // Shell effect for SL, for individual systems  
  float EOscale;  // Scaling factor for even-odd structure in yields  
  int Emode ;      // 0: E over BF_B; 1: E over gs; 2: E_neutron; 12: E_proton  
  float D_Par_Fac;          // Scales the variation of perturbed parameters   
  float _P_DZ_Mean_S1;
  float _P_DZ_Mean_S2;
  float _P_DZ_Mean_S3;     // Shift of mean Z of Mode 3  
  float _P_DZ_Mean_S4;  // Shell for structure at A around 190  
  float  ZC_Mode_4L;  // enhances S1  
  float _P_Z_Curv_S1;
  float P_Z_Curvmod_S1;    // Scales energy-dependent shift   
  float _P_Z_Curv_S2;      
  float _S2leftmod;     // Asymmetry in diffuseness of S2 mass peak   
  float P_Z_Curvmod_S2;    // Scales energy-dependent shift  
  float _P_A_Width_S2;   // A width of Mode 2 (box)  
  float P_Cut_S2;        // Divide S2 into two modes, S2a and S2b  
  float _P_Z_Curv_S3; 
  float P_Z_Curvmod_S3;    // Scales energy-dependent shift  
  float P_Z_Curv_SL4;
  float P_Z_Sigma_SL4; 
  float _P_Z_Curv_S4;
  float P_Z_Curvmod_S4;   // Scales energy-dependent shift  
  float _P_Shell_S1;     // Shell effect for Mode 1 (S1)  
  float _P_Shell_S2;     // Shell effect for Mode 2 (S2)  
  float _P_Shell_S3;     // Shell effect for Mode 3 (SA)  
  float P_Shell_SL4;    // Shell enhancing S1  
  float _P_Shell_S4;    // Shell effect for Mode 4  
  float P_S4_NZmod;     // Variation of S4 shell with N_CN (reference: 180Hg)  
  float PZ_S3_olap_pos;     // Pos. of S3 shell in light fragment (in Z)  
  float PZ_S3_olap_curv; 
  float ETHRESHSUPPS1;    
  float ESIGSUPPS1;      
  float Level_S11;         // Level for mode S11  
  float Shell_fading;    // fading of shell effect with E*  
  float _T_low_S1;  
  float _T_low_S2;       // Slope parameter for tunneling  
  float _T_low_S3;       // Slope parameter for tunneling  
  float _T_low_S4;       // Slope parameter for tunneling  
  float _T_low_SL;       // Slope parameter for tunneling  
  float T_low_S11;      // Slope parameter for tunneling  
  float _P_att_pol;      // Attenuation length of 132Sn shell  
  float P_att_pol2;  
  float P_att_pol3; 
  float _P_att_rel;     // Relative portion of attenuation  
  float dE_Defo_S1;     // Deformation energy expense for Mode 1  
  float dE_Defo_S2;     // Deformation energy expense for Mode 2  
  float dE_Defo_S3;     // Deformation energy expense for Mode 3  
  float dE_Defo_S4;     // Deformation energy expense for Mode 4  
  float betaL0; 
  float betaL1; 
  float betaH0;        // Offset for deformation of heavy fragment  
  float betaH1; 
  float kappa;         // N/Z dedendence of A-asym. potential  
  float TCOLLFRAC;     // Tcoll per energy gain from saddle to scission  
  float ECOLLFRAC; 
  float TFCOLL;   
  float TCOLLMIN; 
  float ESHIFTSASCI_intr;   // Shift of saddle-scission energy   
  float ESHIFTSASCI_coll;    // Shift of saddle-scission energy  
  float EDISSFRAC; 
  float Epot_shift; 
  float SIGDEFO;   
  float SIGDEFO_0; 
  float SIGDEFO_slope; 
  float SIGENECK;        // Width of TXE by fluctuation of neck length  
  float EexcSIGrel; 
  float DNECK;             // Tip distance at scission / fm  
  float FTRUNC50;          // Truncation near Z = 50  
  float ZTRUNC50;          // Z value for truncation  
  float FTRUNC28;          // Truncation near Z = 28  
  float ZTRUNC28;          // Z value for truncation  
  float ZMAX_S2;           // Maximum Z of S2 channel in light fragment  
  float NTRANSFEREO;       // Steps for E sorting for even-odd effect  
  float NTRANSFERE;        // Steps for E sorting for energy division  
  float Csort;             // Smoothing of energy sorting  
  float PZ_EO_symm;        // Even-odd effect in Z at symmetry  
  float PN_EO_Symm;        // Even-odd effect in N at symmetry  
  float R_EO_THRESH;       // Threshold for asymmetry-driven even-odd effect 
  float R_EO_SIGMA; 
  float R_EO_MAX;          // Maximum even-odd effect  
  float _POLARadd;         // Offset for enhanced polarization  
  float POLARfac;          // Enhancement of polarization of ligu. drop  
  float T_POL_RED;         // Reduction of temperature for sigma(Z)  
  float _HOMPOL;           // hbar omega of polarization oscillation  
  float ZPOL1;             // Extra charge polarization of S1  
  float P_n_x;             // Enhanced inverse neutron x section  
  float Tscale; 
  float Econd;    
  float T_orbital;         // From orbital ang. momentum  
  float _Jscaling;         // General scaling of fragment angular momenta  
  float Spin_odd;          // RMS Spin enhancement for odd Z      
  float Esort_extend;      // Extension of energy range for E-sorting  
  float Esort_slope;       // Onset of E-sorting around symmetry       
  float Esort_slope_S0;    // Onset of E-sorting around symmetry for S0 channel       

  // Control parameters: 
  float B_F;                  // Fission barrier  
  float B_F_ld;           // Fission barrier, liquid drop  
  float E_B;              // Outer fission barrier  
  float E_B_ld;           // Outer fission barrier, liquid drop  
  float R_E_exc_Eb;       // Energy above outer barrier  
  float R_E_exc_GS;       // Energy above ground state  
  float P_Z_Mean_S0;      // Mean Z of Mode 1  
  float P_Z_Mean_S1;        // Mean Z of Mode 1  
  float P_Z_Mean_S2;        // Mean Z of Mode 2  
  float P_Z_Mean_S3;       // Mean Z of Mode 3  
  float P_Z_Mean_S4;       // Mean Z of Mode 4  
  float NC_Mode_0;        // Mean N of symm. Mode  
  float NC_Mode_1;        // Mean N of Mode 1  
  float NC_Mode_2;        // Mean N of Mode 2  
  float NC_Mode_3;        // Mean N of Mode 3  
  float NC_Mode_4;
  float B_S1;             // Barrier S1, relative to SL  
  float B_S2;             // Barrier S2, relative to SL  
  float B_S3;             // Barrier S3, relative to SL  
  float B_S4;
  float B_S11;            // Barrier S11, relative to SL  
  float B_S22;            // Barrier S22, relative to SL  
  float DES11ZPM;         // Mod. of eff. barrier due to ZPM in overlap  
  float Delta_NZ_Pol;      // Polarization for 132Sn  
  float Yield_Mode_0;     // Relative yield of SL  
  float Yield_Mode_1;     // Relative yield of S1  
  float Yield_Mode_2;     // Relative yield of S2  
  float Yield_Mode_3;     // Relative yield of S3  
  float Yield_Mode_4;     // Relative yield of S4  
  float Yield_Mode_11;    // Relative yield of S11  
  float Yield_Mode_22;    // Relative yield of S22  
  float P_POL_CURV_S0;    // Stiffnes in N/Z  
  float T_Coll_Mode_0;    // Effective collective temperature  
  float E_Exc_S0;         // Energy over barrier of symmetric channel  
  float E_Exc_S1;         // Energy over barrier of S1 channel  
  float E_Exc_S2;         // Energy over barrier of S2 channel  
  float E_Exc_S3;         // Energy over barrier of S3 channel  
  float E_Exc_S4;         // Energy over barrier of S4 channel  
  float E_Exc_S11;        // Energy over barrier of S11 channel  
  float E_Exc_S22;        // Energy over barrier of S22 channel  
  float E_POT_SCISSION;   // Potential-energy gain saddle-scission  
  float E_diss_Scission;  // Dissipated energy between saddle and scission  
  float EINTR_SCISSION;   // Intrinsic excitation energy at scission  
  float EeffS1;           // Governs S1 reduction  
  float Sigpol_Mode_0;    // Width of isobaric Z distribution 


  //more parameters
  float P_DZ_Mean_S1;
  float P_DZ_Mean_S2;
  float P_DZ_Mean_S3;
  float P_DZ_Mean_S4;
  float P_Z_Curv_S1;
  float P_Z_Curv_S2;
  float S2leftmod;
  float P_A_Width_S2;
  float P_Z_Curv_S3;
  float P_Z_Curv_S4;
  float Delta_S0;
  float P_Shell_S1;
  float P_Shell_S2;
  float P_Shell_S3;
  float P_Shell_S4;
  float T_low_S1;
  float T_low_S2;
  float T_low_S3;
  float T_low_S4;
  float T_low_SL;
  float P_att_pol;
  float P_att_rel;
  float HOMPOL;
  float POLARadd;  
  float Jscaling;
  
  float Escission_lim;
      
  float Spin_CN;  
  float Spin_pre_fission;
  float Spin_gs_light;
  float Spin_gs_heavy;

  float R_E_intr_S1, R_E_intr_S2, R_E_intr_S3;   // intrinsic exc. energies at barrier
  float R_E_intr_S4;
  float R_Att[7];                              // attenuation of shell
  float R_Att_Sad[7];     

  float Etot,E1FG,E1ES;
  float Rincr1P,Rincr1N,Rincr2,Rincr2P,Rincr2N;
  float T1,T2,E1,E2;
  float E_coll_saddle[7];
  float Ediff;
  float DT;

  float R_A_help,RN;
  int Iguess;
	
  float R_N_heavy;
  int I_Z_sad, I_N_sad, I_A_sad;
  int I_Z_heavy_sad, I_A_heavy_sad;
  int I_Z_light_sad, I_A_light_sad;
  int I_N_heavy_sad, I_N_light_sad;
  int I_Z_sci, I_N_sci, I_A_sci;    //sci: values at scission, before prompt-neutron emission
  int I_Z_light_sci, I_N_light_sci, I_A_light_sci;
  int I_Z_heavy_sci, I_N_heavy_sci, I_A_heavy_sci;
  int I_Z_post, I_N_post, I_A_post; 
  int I_Z_light_post, I_N_light_post, I_A_light_post;
  int I_Z_heavy_post, I_N_heavy_post, I_A_heavy_post;
  float ESIGDEFOlight,ESIGDEFOheavy;
  float RS;

  float I_nu_ss;
  //Input/Return parameters of subroutine Eva:
  float Array_E_n1_frag1[51];  // neutron energy array in fragment frame
  float Array_v_n1_frag1[51];
  float Array_vlong_n1_frag1[51];  // neutron velocity array longitudinal
  float Array_vperp_n1_frag1[51];  // neutron velocity array perpendicular
  float Array_E_n2_frag2[51];  // neutron energy array in fragment frame
  float Array_v_n2_frag2[51];
  float Array_vlong_n2_frag2[51];  // neutron velocity array longitudinal
  float Array_vperp_n2_frag2[51];  // neutron velocity array perpendicular
  float Array_Tn[51];  // neutron decay times after scission
  float Array_Eg0_light[101];  // statistical gamma energy array
  float Array_Eg0_heavy[101];  // statistical gamma energy array

  float TXE,Erotlight,Erotheavy,TXElight,TXEheavy;
  float IfragEff_light,IfragEff_heavy;

  float Thcm_light,Thcm_heavy;
  float Phcm_light,Phcm_heavy;
  float KElab_light, KElab_heavy;
  float Glab_light, Glab_heavy;
  float Blab_light, Blab_heavy;
  float vlab_light,vlab_heavy;
  int Qffl,Qffh;
  float Thlab_light,Thlab_heavy;
  float Phlab_light,Phlab_heavy;
  float Brhoffl,Brhoffh;
  //float Thlablab_light,Thlablab_heavy;
  //float Phlablab_light,Phlablab_heavy;
  float PF;
  
};
#endif
