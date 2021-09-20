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

#include "GEF.h"

GEF::GEF(NPL::Particle FissNucl)
{
  FIS = FissNucl;
  Afis = FIS.GetA();
  Zfis = FIS.GetZ();

  ReadBEldmTF();
  ReadBEexp();
  ReadDEFOtab();
  ReadShellMO();
  ReadNucProp();
  ReadGEFParameters();

  pi = TMath::Pi();
  rn = new TRandom3();
  Init();
  FissionGENERAL();
}

GEF::~GEF(void){}

void GEF::Init(void)
{
  Ngtot = 0;
  Nglight = 0;
  Ngheavy = 0;
  Egtot1000 = 0;
  DZ_S2_lowE = 0;    

  N_Error_Max = 10;       
  I_Error = 0;       
  I_DelGam = 0;
  P_Z_CN = Zfis;       
  P_A_CN = Afis;       
  P_E_exc = 0;               
  P_I_rms_CN = 0;
  I_A_CN = P_A_CN;
  I_Z_CN = P_Z_CN;

  // Use nominal parameter values:
  P_DZ_Mean_S1 = _P_DZ_Mean_S1;
  P_DZ_Mean_S2 = _P_DZ_Mean_S2;
  P_DZ_Mean_S3 = _P_DZ_Mean_S3;
  P_DZ_Mean_S4 = _P_DZ_Mean_S4;
  P_Z_Curv_S1 = _P_Z_Curv_S1;
  P_Z_Curv_S2 = _P_Z_Curv_S2;
  S2leftmod = _S2leftmod;
  P_A_Width_S2 = _P_A_Width_S2;
  P_Z_Curv_S3 = _P_Z_Curv_S3;
  P_Z_Curv_S4 = _P_Z_Curv_S4;
  Delta_S0 = _Delta_S0;
  P_Shell_S1 = _P_Shell_S1;
  P_Shell_S2 = _P_Shell_S2;
  P_Shell_S3 = _P_Shell_S3;
  P_Shell_S4 = _P_Shell_S4;
  T_low_S1 = _T_low_S1;
  T_low_S2 = _T_low_S2;
  T_low_S3 = _T_low_S3;
  T_low_S4 = _T_low_S4;
  T_low_SL = _T_low_SL;
  P_att_pol = _P_att_pol;
  P_att_rel = _P_att_rel;
  HOMPOL = _HOMPOL;
  POLARadd = _POLARadd;
  Jscaling = _Jscaling;


  // Control parameters: 
  B_F = 0;              // Fission barrier  
  B_F_ld = 0;           // Fission barrier, liquid drop  
  E_B = 0;              // Outer fission barrier  
  E_B_ld = 0;           // Outer fission barrier, liquid drop  
  R_E_exc_Eb = 0;       // Energy above outer barrier  
  R_E_exc_GS = 0;       // Energy above ground state  
  P_Z_Mean_S0 = 0;      // Mean Z of Mode 1  
  P_Z_Mean_S1 = 52.8;   // Mean Z of Mode 1  
  P_Z_Mean_S2 = 55;     // Mean Z of Mode 2  
  P_Z_Mean_S3 = 65;     // Mean Z of Mode 3  
  P_Z_Mean_S4 = 42.05;  // Mean Z of Mode 4  
  NC_Mode_0 = 0;        // Mean N of symm. Mode  
  NC_Mode_1 = 0;        // Mean N of Mode 1  
  NC_Mode_2 = 0;        // Mean N of Mode 2  
  NC_Mode_3 = 0;        // Mean N of Mode 3  
  NC_Mode_4 = 0;
  B_S1 = 0;             // Barrier S1, relative to SL  
  B_S2 = 0;             // Barrier S2, relative to SL  
  B_S3 = 0;             // Barrier S3, relative to SL  
  B_S4 = 0;
  B_S11 = 0;            // Barrier S11, relative to SL  
  B_S22 = 0;            // Barrier S22, relative to SL  
  DES11ZPM = 0;         // Mod. of eff. barrier due to ZPM in overlap  
  Delta_NZ_Pol = 0;      // Polarization for 132Sn  
  Yield_Mode_0 = 0;     // Relative yield of SL  
  Yield_Mode_1 = 0;     // Relative yield of S1  
  Yield_Mode_2 = 0;     // Relative yield of S2  
  Yield_Mode_3 = 0;     // Relative yield of S3  
  Yield_Mode_4 = 0;     // Relative yield of S4  
  Yield_Mode_11 = 0;    // Relative yield of S11  
  Yield_Mode_22 = 0;    // Relative yield of S22  
  P_POL_CURV_S0 = 0;    // Stiffnes in N/Z  
  T_Coll_Mode_0 = 0;    // Effective collective temperature  
  E_Exc_S0 = 0;         // Energy over barrier of symmetric channel  
  E_Exc_S1 = 0;         // Energy over barrier of S1 channel  
  E_Exc_S2 = 0;         // Energy over barrier of S2 channel  
  E_Exc_S3 = 0;         // Energy over barrier of S3 channel  
  E_Exc_S4 = 0;         // Energy over barrier of S4 channel  
  E_Exc_S11 = 0;        // Energy over barrier of S11 channel  
  E_Exc_S22 = 0;        // Energy over barrier of S22 channel  
  E_POT_SCISSION = 0;   // Potential-energy gain saddle-scission  
  E_diss_Scission = 0;  // Dissipated energy between saddle and scission  
  EINTR_SCISSION = 0;   // Intrinsic excitation energy at scission  
  EeffS1 = 0;           // Governs S1 reduction  
  Sigpol_Mode_0 = 0;    // Width of isobaric Z distribution 
}

void GEF::ReadBEldmTF(void)
{
  char filename[256];
  char *fname;
  Float_t temp;
  for(int i=0;i<204;i++)
    for(int j=0;j<137;j++)
      BEldmTF[i][j]=0;
  int counter=0;

  fname = getenv("NPTOOL");
  sprintf(filename,"%s/NPLib/Physics/PAR_GEF/BEldmTF.dat",fname);

  ifstream file(filename);
  string line;	
  if (!file) 
    cout << "Couldn't open file : " << filename << endl ; 
  else 
  {
    //cout << "File : " << filename << endl; 
    for(int i=0;i<3;i++)
    {
      getline(file,line);
      //cout<<line<<endl;
    }
    while(getline(file,line))
    {
      istringstream read (line);
      //cout<<line<<endl;
      while(read)
      {
        read>>temp;
        if((counter/136+1)<204&&(counter%136+1)<137)
        {
          BEldmTF[counter/136+1][counter%136+1] = temp;
          //cout<<"i: "<< counter/136+1<< " j: "<<counter%136+1<<" "<<BEldmTF[counter/136+1][counter%136+1]<<endl;
        }
        counter++;  
      }
      counter--; //while(read) it is satisifed one more time at the end of each line
    }
  } 
}

void GEF::ReadBEexp(void)
{
  char filename[256];
  char *fname;
  int zz, aa;
  float temp;
  for(int i=0;i<204;i++)
    for(int j=0;j<137;j++)
      BEexp[i][j]=-1.e11;

  fname = getenv("NPTOOL");
  sprintf(filename,"%s/NPLib/Physics/PAR_GEF/BEexp.dat",fname);

  ifstream file(filename);
  string line;	
  if (!file) 
    cout << "Couldn't open file : " << filename << endl ; 
  else 
  {
    //cout << "File : " << filename << endl; 
    for(int i=0;i<3;i++)
    {
      getline(file,line);
      //cout<<line<<endl;
    }
    while(getline(file,line))
    {
      istringstream read (line);
      read>>zz>>aa>>temp;
      BEexp[aa-zz][zz]=temp;
      //cout<<zz<<" "<<aa<<" "<<BEexp[aa-zz][zz]<<endl;
    }
  }
}

void GEF::ReadDEFOtab(void)
{
  char filename[256];
  char *fname;
  int zz, nn;
  float temp1,temp2NotUsed;
  for(int i=0;i<237;i++)
    for(int j=0;j<137;j++)
      DEFOtab[i][j]=0.;

  fname = getenv("NPTOOL");
  sprintf(filename,"%s/NPLib/Physics/PAR_GEF/DEFO.dat",fname);

  ifstream file(filename);
  string line;	
  if (!file) 
    cout << "Couldn't open file : " << filename << endl ; 
  else 
  {
    //cout << "File : " << filename << endl; 
    while(getline(file,line))
    {
      istringstream read (line);
      read>>zz>>nn>>temp1>>temp2NotUsed;
      DEFOtab[nn][zz]=temp1;
      //cout<<zz<<" "<<nn<<" "<<DEFOtab[nn][zz]<<endl;
    }
  } 
}

void GEF::ReadShellMO(void)
{
  char filename[256];
  char *fname;
  Float_t temp;
  for(int i=0;i<204;i++)
    for(int j=0;j<137;j++)
      ShellMO[i][j]=0;
  int counter=0;

  fname = getenv("NPTOOL");
  sprintf(filename,"%s/NPLib/Physics/PAR_GEF/ShellMO.dat",fname);

  ifstream file(filename);
  string line;	
  if (!file) 
    cout << "Couldn't open file : " << filename << endl ; 
  else 
  {
    //cout << "File : " << filename << endl; 
    for(int i=0;i<2;i++)
    {
      getline(file,line);
      //cout<<line<<endl;
    }
    while(getline(file,line))
    {
      istringstream read (line);
      //cout<<line<<endl;
      while(read)
      {
        read>>temp;
        if((counter/136+1)<204&&(counter%136+1)<137)
        {
          ShellMO[counter/136+1][counter%136+1] = temp;
          //cout<<"i: "<< counter/136+1<< " j: "<<counter%136+1<<" "<<ShellMO[counter/136+1][counter%136+1]<<endl;
        }
        counter++;  
      }
      counter--; //while(read) it is satisifed one more time at the end of each line
    }
  } 
}

void GEF::ReadNucProp(void)
{
  char filename[256];
  char *fname;
  int index,zz, aa,iso;
  float spin;
  for(int i=0;i<112;i++)
    for(int j=0;j<280;j++)
      R_SPI[i][j]=0.;

  fname = getenv("NPTOOL");
  sprintf(filename,"%s/NPLib/Physics/PAR_GEF/NucPropNUBASE.dat",fname);

  ifstream file(filename);
  string line;	
  if (!file) 
    cout << "Couldn't open file : " << filename << endl ; 
  else 
  {
    for(int i=0;i<9;i++)
      getline(file,line);
    while(getline(file,line))
    {
      istringstream read (line);
      read>>index>>zz>>aa>>iso>>spin;
      if(iso==0)
        R_SPI[zz][aa]=spin;
    }
  } 
}

void GEF::ReadGEFParameters(void)
{
  char filename[256];
  char *fname;
  char temp[256];
  float PARAM[96];
  for(int i=0;i<96;i++)
    PARAM[i]=0.;

  fname = getenv("NPTOOL");
  sprintf(filename,"%s/NPLib/Physics/PAR_GEF/GEFParameters.dat",fname);

  ifstream file(filename);
  string line;	
  if (!file) 
    cout << "Couldn't open file : " << filename << endl ; 
  else 
  {
    //cout << "File : " << filename << endl; 
    for(int i=0;i<96;i++)
    {
      getline(file,line);
      istringstream read (line);
      read>>temp>>PARAM[i];
      //cout<<PARAM[i]<<endl;
    }
  } 

  Emax_valid = PARAM[0];      
  Eexc_min_multi = PARAM[1];             
  _Delta_S0 = PARAM[2];          
  EOscale = PARAM[3];   
  Emode = int(PARAM[4]);        
  D_Par_Fac =  PARAM[5];             
  _P_DZ_Mean_S1 = PARAM[6];
  _P_DZ_Mean_S2 = PARAM[7];
  _P_DZ_Mean_S3 = PARAM[8];     
  _P_DZ_Mean_S4 = PARAM[9]; 
  ZC_Mode_4L = PARAM[10];   
  _P_Z_Curv_S1 = PARAM[11];
  P_Z_Curvmod_S1 = PARAM[12];       
  _P_Z_Curv_S2 = PARAM[13];      
  _S2leftmod = PARAM[14];        
  P_Z_Curvmod_S2 = PARAM[15];     
  _P_A_Width_S2 = PARAM[16];    
  P_Cut_S2 = PARAM[17];         
  _P_Z_Curv_S3 = PARAM[18]; 
  P_Z_Curvmod_S3 = PARAM[19];     
  P_Z_Curv_SL4 = PARAM[20];
  P_Z_Sigma_SL4 = PARAM[21]; 
  _P_Z_Curv_S4 = PARAM[22];
  P_Z_Curvmod_S4 = PARAM[23];     
  _P_Shell_S1 = PARAM[24];       
  _P_Shell_S2 = PARAM[25];       
  _P_Shell_S3 = PARAM[26];       
  P_Shell_SL4 = PARAM[27];    
  _P_Shell_S4 = PARAM[28];     
  P_S4_NZmod = PARAM[29];      
  PZ_S3_olap_pos = PARAM[30];       
  PZ_S3_olap_curv = PARAM[31]; 
  ETHRESHSUPPS1 = PARAM[32];    
  ESIGSUPPS1 = PARAM[33];      
  Level_S11 = PARAM[34];          
  Shell_fading = PARAM[35];      
  _T_low_S1 = PARAM[36];  
  _T_low_S2 = PARAM[37];         
  _T_low_S3 = PARAM[38];         
  _T_low_S4 = PARAM[39];         
  _T_low_SL = PARAM[40];         
  T_low_S11 = PARAM[41];        
  _P_att_pol = PARAM[42];        
  P_att_pol2 = PARAM[43];  
  P_att_pol3 = PARAM[44]; 
  _P_att_rel = PARAM[45];      
  dE_Defo_S1 = PARAM[46];       
  dE_Defo_S2 = PARAM[47];       
  dE_Defo_S3 = PARAM[48];       
  dE_Defo_S4 = PARAM[49];     
  betaL0 = PARAM[50]; 
  betaL1 = PARAM[51]; 
  betaH0 = PARAM[52];         
  betaH1 = PARAM[53]; 
  kappa = PARAM[54];         
  TCOLLFRAC = PARAM[55];       
  ECOLLFRAC = PARAM[56]; 
  TFCOLL = PARAM[57];   
  TCOLLMIN = PARAM[58]; 
  ESHIFTSASCI_intr = PARAM[59];      
  ESHIFTSASCI_coll = PARAM[60];      
  EDISSFRAC = PARAM[61]; 
  Epot_shift = PARAM[62]; 
  SIGDEFO = PARAM[63];   
  SIGDEFO_0 = PARAM[64]; 
  SIGDEFO_slope = PARAM[65]; 
  SIGENECK = PARAM[66];         
  EexcSIGrel = PARAM[67]; 
  DNECK = PARAM[68];              
  FTRUNC50 = PARAM[69];            
  ZTRUNC50 = PARAM[70];            
  FTRUNC28 = PARAM[71];           
  ZTRUNC28 = PARAM[72];            
  ZMAX_S2 = PARAM[73];             
  NTRANSFEREO = PARAM[74];         
  NTRANSFERE = PARAM[75];          
  Csort = PARAM[76];             
  PZ_EO_symm = PARAM[77];         
  PN_EO_Symm = PARAM[78];         
  R_EO_THRESH = PARAM[79];        
  R_EO_SIGMA = PARAM[80]; 
  R_EO_MAX = PARAM[81];            
  _POLARadd = PARAM[82];           
  POLARfac = PARAM[83];            
  T_POL_RED = PARAM[84];         
  _HOMPOL = PARAM[85];             
  ZPOL1 = PARAM[86];              
  P_n_x = PARAM[87];              
  Tscale = PARAM[88]; 
  Econd = PARAM[89];    
  T_orbital = PARAM[90];          
  _Jscaling = PARAM[91];           
  Spin_odd = PARAM[92];                
  Esort_extend = PARAM[93];        
  Esort_slope = PARAM[94];              
  Esort_slope_S0 = PARAM[95];           

}

float GEF::Max(float a, float b)
{
  if(a>b)
    return a;
  else
    return b;
}
float GEF::Min(float a, float b)
{
  if(a<b)
    return a;
  else
    return b;
}
float GEF::Erf(float x)
{
  /* Sergei Winitzki, 2008: relative accuracy < 1.4E-4 */
  //  Dim As Single a = 0.147
  //  Const As Single b = 1.27324  // 4/pi,
  //  Const As Single pi = 3.14159
  //  Dim As Single Result
  //  Result = Sqr(1.E0 - exp(-x^2 * (b + a * x^2) / (1.E0 - a * x^2) ) )
  //  If x < 0 Then Result = - Result
  //  WiErf = Result
  return 1. - Erfc(x);
}
float GEF::Erfc(float x)
{
  /* Complementary error function from numerical recipes */
  double t,z,r;
  z = abs(x);
  t = 1./(1. + 0.5 * z);
  r = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  if (x < 0)
    return 2. - r;
  else
    return r;
}
float GEF::Tanh(float x)
{
  float Result;
  if( x >= 0)
    Result = (1 - exp(-2*x))/(1 + exp(-2*x));
  else
    Result = (exp(2*x) - 1)/(exp(2*x) + 1);
  return Result;
}
float GEF::Coth(float x)
{
  float Result;
  if(x >= 0)
    Result = (1 + exp(-2*x))/(1 - exp(-2*x));
  else
    Result = (exp(2*x) + 1)/(exp(2*x) - 1);
  return Result;
}
float GEF::Log10(float R)
{
  return log(R) / log(10);  
}
float GEF::Bell(float xpos, float xleft, float xright)
{
  // Bell-shaped curve, maximum = 1, zero for x < xleft and x > xright
  float x,y;
  if((xpos - xleft) > 1.e-3 && (xright - xpos) > 1.e-3)
  {
    x = (xpos - xleft) / (xright - xleft);
    y = (pow(x,2) * pow(1-x,2))/0.0625;
  }
  else
    y = 0;
  return y;  
}
float GEF::Csng(int a)
{
  return 1.*a;
}
int GEF::CInt(float a)
{
  return int(a);
}
float GEF::Mod(float a, float b)
{
  return a - int(a/b)*b;
}
float GEF::FLOOR(float a)
{
  return int(a);
}
float GEF::CEIL(float a)
{
  if (a>FLOOR(a)) 
    return FLOOR(a)+1;
  else
    return a;
}
int GEF::Sgn(float a)
{
  if(a!=0)
    return int(a/abs(a));
  else
    return 0;
}

void GEF::FissionGENERAL()
{
  float ZsqrA_fis;

  ZsqrA_fis = Csng(I_Z_CN^2) / Csng(I_A_CN);
  //ZsqrA_fis = pow(I_Z_CN,2)/I_A_CN;
  //Escission_lim = 0;//900.0 * exp(-ZsqrA_fis/13.0); // from PRC 86 (2012) 034605
  Escission_lim = 900.0 * exp(-ZsqrA_fis/13.0); // from PRC 86 (2012) 034605
  CentralZ_FM();
  I_N_CN = I_A_CN - I_Z_CN;
  MeanDefScission();
  MeanZ_asA();
  RelationsZ_A_FM();
  FissionBarriers();
  PolarizationStiffness();

}

void GEF::CentralZ_FM(void)
{
  /* Central Z values of fission modes */    
  /* Test of barrier height for compact fission */
  /*   Dim As Single Fbarr_shell
       Dim As Single Fbarr_macro
       Fbarr_shell = 1.44 * (I_C_CN/2)^2 / (1.4 * 2 * (I_A_CN/2)^(1/3))
       Fbarr_macro = 1.44 * (I_C_CN/2)^2 / (1.4 * 2 * (I_A_CN/2)^(1/3)) */ 

  /* Fit to positions of fission channels (Boeckstiegel et al., 2008) */
  /* P_DZ_Mean_S1 and P_DZ_Mean_S2 allow for slight adjustments */

  float R_Z_mod, R_A_mod,R_corr2;
  R_Z_mod = Csng(I_Z_CN);
  R_A_mod = Csng(I_A_CN);
  R_corr2 = 0.055 * (R_A_mod - R_Z_mod*236/92) ;
  // * SL (S0) : *
  ZC_Mode_0 = R_Z_mod * 0.5;      /* Central Z value of SL mode */
  // * S1 : *
  ZC_Mode_1 = (55.8 - 54.5) / (1.56 - 1.50) * 
    (pow(R_Z_mod,1.3) / R_A_mod - 1.50) + 51.5 + P_DZ_Mean_S1 + R_corr2;
  // * S2: *          
  ZC_Mode_2 = (55.8 - 54.5) / (1.56 - 1.50) * 
    (pow(R_Z_mod,1.3) / R_A_mod - 1.50) + 54.5 + P_DZ_Mean_S2 + R_corr2;
  // * S3: * 
  //      ZC_Mode_3 = ZC_Mode_2 + 4.5E0 + P_DZ_Mean_S3   
  ZC_Mode_3 = ZC_Mode_2 + 4.87 + P_DZ_Mean_S3;   
  ZC_Mode_3_shift = - 0.015 * (ZC_Mode_3 - ZC_Mode_0);  

  ZC_Mode_3 = ZC_Mode_2 + 5.5 + P_DZ_Mean_S3;   
  ZC_Mode_3_shift = - 0.035 * (ZC_Mode_3 - ZC_Mode_0);  

  ZC_Mode_3 = ZC_Mode_3 + ZC_Mode_3_shift;
  // To account for ZCN-dependent shift of S3 towards minimum of mac. potential
  // * S4: *
  // Do not delete these lines (,because this is a very good fit!):
  //    ZC_Mode_4 = 38.5 + (I_A_CN-I_Z_CN-110)*0.12 - (I_A_CN-I_Z_CN-110)^2 * 0.009 _
  //                - (I_Z_CN-77)*0.34 + P_DZ_Mean_S4 
  // assumption: mode position moves with Z and A (adjusted to exp. data
  // of Itkis and Andreyev et al.
  //  ZC_Mode_4 = - (55.8E0 - 54.5E0) / (1.56E0 - 1.50E0) * _
  //              (R_Z_mod^1.3E0 / I_A_CN - 1.50E0) + 37.5 + P_DZ_Mean_S4  - R_corr2 _
  //                + 0.035 * (I_A_CN- I_Z_CN - 100)   //fits Tl201 (itkis), but not so well Po194,196 (Andreyev) 
  // corresponding P_DZ_Mean_S4 = -0.5 
  ZC_Mode_4 = 35.5 + Csng(I_A_CN - I_Z_CN - 100) * 0.11 + P_DZ_Mean_S4;   // mean Z of light fragment in Mode 4

  // Print "R_corr2: ",R_corr2                
  // Print "Mode 2, Z and N:", ZC_Mode_2, ZC_Mode_2 / P_Z_CN * (P_A_CN-P_Z_CN)                
  // Print "Mode 4, Z and N:", ZC_Mode_4, ZC_Mode_4 / P_Z_CN * (P_A_CN-P_Z_CN) 

  /*>*/
  P_Z_Mean_S0 = ZC_Mode_0;  /* Copy to global parameter */
  P_Z_Mean_S1 = ZC_Mode_1;  /* Copy to global parameter */
  P_Z_Mean_S2 = ZC_Mode_2;  /*             "            */
  P_Z_Mean_S3 = ZC_Mode_3;  /*             "            */
  P_Z_Mean_S4 = ZC_Mode_4;  /*             "            */
  /*<*/
}

void GEF::MeanDefScission(void)
{
  /* Mean deformation at scission as a function of mass */
  /* Mode 0: liquid drop */
  beta1_prev = 0.3;
  beta2_prev = 0.3;
  beta1_opt = beta1_prev;
  beta2_opt = beta2_prev;
  for(int i = 10;i<=(I_Z_CN - 10);i++)
  {
    IZ1 = i;
    Z1 = Csng(IZ1);
    IZ2 = I_Z_CN - IZ1;
    Z2 = Csng(IZ2);
    A1 = Z1 / Csng(I_Z_CN) * Csng(I_A_CN);
    A2 = I_A_CN - A1;

    // Deformed shell below Z=50, valid for S0 at low E*
    Beta[0][1][IZ1] = beta_light(IZ1,betaL0,betaL1) - 0.1;
    // Lower deformation than S1/S2, because Coulomb repulsion from deformed heavy fragment is weaker.
    Beta[0][2][IZ2] = beta_light(IZ2,betaL0,betaL1) - 0.1;     //          "
    E_defo = LyMass(Z1,A1,Beta[0][1][IZ1]) - LyMass(Z1,A1,0.0);
    Edefo[0][1][IZ1] = E_defo;
    E_defo = LyMass(Z2,A2,Beta[0][2][IZ2]) - LyMass(Z2,A2,0.0);
    Edefo[0][2][IZ2] = E_defo;

    Beta_Equi(A1,A2,Z1,Z2,DNECK,beta1_prev,beta2_prev,beta1_opt,beta2_opt);

    // Deformation by macroscopic model, valid for S0 at high E*
    //Print "Mode 0, Z1,Z2,beta1,beta2 ";Z1;" ";Z2;" ";beta1_opt,beta2_opt
    //Print Z1;" ";Z2;" ";beta1_opt,beta2_opt
    Beta[1][1][IZ1] = beta1_opt; /* "light" fragment */
    //      Beta(4,1,IZ1) = beta1_opt
    Beta[1][2][IZ2] = beta2_opt; /* "heavy" fragment */
    //      Beta(4,2,IZ2) = beta2_opt
    beta1_prev = beta1_opt;
    beta2_prev = beta2_opt;
    E_defo = LyMass(Z1,A1,beta1_opt) - LyMass(Z1,A1,0.0);
    Edefo[1][1][IZ1] = E_defo;  /* "light" fragment */
    //      Edefo(4,1,IZ1) = E_defo
    E_defo = LyMass(Z2,A2,beta2_opt) - LyMass(Z2,A2,0.0);
    Edefo[1][2][IZ2] = E_defo;  /* "heavy" fragment */
    //      Edefo(4,2,IZ2) = E_defo
  }

  /* Mode 1: deformed shells (light) and spherical (heavy) */
  for(int i = 10;i<=(I_Z_CN - 10);i++)
  {
    Z1 = Csng(i);
    Z2 = Csng(I_Z_CN) - Z1;
    A1 = (Z1 - 0.5) / Csng(I_Z_CN) * Csng(I_A_CN); /* polarization roughly considered */
    A2 = Csng(I_A_CN) - A1;
    if(I_Z_CN * 0.5 < ZC_Mode_1)
    {
      // Beta_opt_light(A1,A2,Z1,Z2,dneck,0,rbeta_ld)
      /* nu_mean of Cf requires shells in the light fragment: */
      //     rbeta = beta_light(I,betaL0,betaL1) - 0.1 
      // smaller than general deformation of light fragment   
      //        (less neck influence due to spherical heavy fragment)
      rbeta = beta_light(i,betaL0,betaL1);        
      if(rbeta < 0)
        rbeta = 0;
    }
    else
    {
      rbeta = beta_heavy(i,betaH0,betaH1);  // equal to S2 channel
      if(rbeta < 0)
        rbeta = 0;
    }
    Beta[2][1][i] = rbeta;    /* "light" fragment */
    E_defo = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.0);
    Edefo[2][1][i] = E_defo; /* "light" fragment */
  }

  for(int i = 10;i<=(I_Z_CN - 10);i++)
  {
    rbeta = 0;
    Beta[2][2][i] = rbeta;
    Edefo[2][2][i] = 0;   /* "heavy" fragment (at S1 shell) */
  }
  /* Mode 2: deformed shells (light and heavy) */
  for(int i = 10;i<=(I_Z_CN - 10);i++)
  {
    Z1 = Csng(i);
    Z2 = Csng(I_Z_CN) - Z1;
    A1 = (Z1 - 0.5) / Csng(I_Z_CN) * Csng(I_A_CN); /* polarization roughly considered */
    A2 = Csng(I_A_CN) - A1;
    if(I_Z_CN * 0.5 < ZC_Mode_2)
    {
      // Beta_opt_light(A1,A2,Z1,Z2,dneck,beta_heavy(Z2),rbeta_ld)
      rbeta = beta_light(i,betaL0,betaL1);   // general deformation of light fragment
      if(rbeta < 0)
        rbeta = 0;  // negative values replaced by 0
    }
    else
      rbeta = beta_heavy(i,betaH0,betaH1);  // equal to S2 channel
    Beta[3][1][i] = rbeta;
    E_defo = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.0);
    Edefo[3][1][i] = E_defo;
  }
  for(int i = 10;i<=(I_Z_CN - 10);i++)
  {
    rbeta = beta_heavy(i,betaH0,betaH1);   /* "heavy" fragment (at S2 shell)*/
    if(rbeta < 0)
      rbeta = 0;  // negative values replaced by 0  
    Beta[3][2][i] = rbeta;
    Z1 = Csng(i);
    A1 = (Z1 + 0.5) / Csng(I_Z_CN) * Csng(I_A_CN); /* polarization roughly considered */
    E_defo = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.0);
    Edefo[3][2][i] = E_defo;
  }
  /* Mode 3 */
  for(int i = 10;i<=(I_Z_CN - 10);i++)
  {
    Z1 = Csng(i);
    Z2 = Csng(I_Z_CN) - Z1;
    A1 = (Z1 - 0.5) / Csng(I_Z_CN) * Csng(I_A_CN); /* polarization roughly considered */
    A2 = Csng(I_A_CN) - A1;
    rbeta = beta_light(i,betaL0,betaL1); 
    rbeta = Max(rbeta-0.10,0.0);  /* for low nu-bar of lightest fragments */
    //  Beta_opt_light(A1,A2,Z1,Z2,dneck,beta_heavy(Z2,betaH0,betaH1),rbeta)  
    Beta[4][1][i] = rbeta;
    E_defo = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.0);
    Edefo[4][1][i] = E_defo;
  }
  for(int i = 10;i<=(I_Z_CN - 10);i++)
  {
    rbeta = beta_heavy(i,betaH0,betaH1) + 0.3; // Shift from isotopic distributions of S3 nuclei in 239Pu(nth,f)  
    if(rbeta < 0)
      rbeta = 0;
    Beta[4][2][i] = rbeta;
    Z1 = Csng(i);
    A1 = (Z1 + 0.5) / Csng(I_Z_CN) * Csng(I_A_CN); /* polarization roughly considered */
    E_defo = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.0);
    Edefo[4][2][i] = E_defo;
  }
  /* Mode 4: (Channel S4, (fissioning nucleus in the Z=80 region), Zfrag around 38 */
  for(int i = 10;i<=(I_Z_CN - 10);i++)    // heavy fragment
  {
    Z1 = Csng(i);
    A1 = Z1 / Csng(I_Z_CN) * Csng(I_A_CN); /* charge polarization neglected */
    rbeta = Beta[3][1][i];  // Deformation like the light fragment of S2 in the actinides
    if(rbeta < 0)
      rbeta = 0;
    Beta[5][1][i] = rbeta;
    Beta[5][2][i] = rbeta;
    E_defo = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.0);
    Edefo[5][1][i] = E_defo; /* light fragment */
    Edefo[5][2][i] = E_defo; /* heavy fragment */
  }
  /* Mode 5: (Channel ST1 in both fragments) */
  for(int i = 10; i<=(I_Z_CN - 10);i++)
  {
    Z1 = Csng(i);
    Z2 = Csng(I_Z_CN) - Z1;
    rbeta = Beta[2][2][i];
    if(rbeta < 0)
      rbeta = 0;
    Beta[6][1][int(Z1)] = rbeta;
    Beta[6][2][int(Z1)] = rbeta;
  }

  /* Mode 6: (Channel ST2 in both fragments) */
  for(int i=10;i<=(I_Z_CN - 10);i++)
  {
    Z1 = Csng(i);
    Z2 = Csng(I_Z_CN) - Z1;
    rbeta = Beta[3][2][i];
    if(rbeta < 0)
      rbeta = 0;
    Beta[7][1][int(Z1)] = rbeta;
    Beta[7][2][int(Z1)] = rbeta;
  }  
}

void GEF::MeanZ_asA(void)
{
  /* Mean Z as a function of mass */
  /* Mode 0 */
  for(int I = 10;I<=(I_A_CN - 10);I++)
  {
    ZUCD = Csng(I) / Csng(I_A_CN) * Csng(I_Z_CN);
    beta1 = Beta[1][1][int(ZUCD + 0.5)];
    beta2 = Beta[1][2][int(I_Z_CN - ZUCD + 0.5)];
    Z1 = Z_equi(I_Z_CN,I, I_A_CN - I, beta1, beta2, DNECK,0);
    Zmean[0][1][I] = Z1;
    Zshift[0][1][I] = Z1 - ZUCD;
    Zmean[0][2][I_A_CN - I] = I_Z_CN - Z1;
    Zshift[0][2][I_A_CN - I] = ZUCD - Z1;
    ZshiftOriginal[0][1][I]=Zshift[0][1][I];
    ZshiftOriginal[0][2][I_A_CN - I]=Zshift[0][2][I_A_CN - I];     
  }
  /* Mode 1 */
  for(int I = 10;I<=(I_A_CN - 10);I++)
  {
    ZUCD = Csng(I) / Csng(I_A_CN) * Csng(I_Z_CN);
    Z = ZUCD + ZPOL1; /* Charge polarisation is considered in a crude way */
    beta1 = Beta[2][1][CInt(Z)]; /* "light" fragment */
    Z = ZUCD - ZPOL1;
    beta2 = Beta[2][2][CInt(I_Z_CN-Z)]; /* "heavy" fragment  at S1 shell */
    if(Csng(I_Z_CN) * 0.5 < ZC_Mode_1)
    {
      Z1 = Z_equi(I_Z_CN,I, I_A_CN - I, beta1, beta2, DNECK,1); 
      Z1 = Z1 + POLARadd;
    }
    else
      Z1 = Z_equi(I_Z_CN,I, I_A_CN - I, beta1, beta2, DNECK,1);
    Z1 = Z1 + ZPOL1;  /* Charge polarization by shell */
    if((I_Z_CN - Z1) < 50. && (Csng(I_Z_CN) - Z1) > Z1) 
      Z1 = Csng(I_Z_CN) - 50.;    /* Z of mean heavy fragment not below 50 */
    Zmean[1][1][I] = Z1;
    Zshift[1][1][I] = Z1 - ZUCD;     // neutron-deficient
    Zmean[1][2][I_A_CN - I] = Csng(I_Z_CN) - Z1;
    Zshift[1][2][I_A_CN - I] = ZUCD - Z1;  // neutron rich at shell
    ZshiftOriginal[1][1][I]=Zshift[1][1][I];     
    ZshiftOriginal[1][2][I_A_CN - I]=Zshift[1][2][I_A_CN - I];     
  }
  /* Mode 2 */
  for(int I = 10;I<=(I_A_CN - 10);I++)
  {
    ZUCD = Csng(I) / Csng(I_A_CN) * Csng(I_Z_CN);
    Z = ZUCD; /* Charge polarisation is here neglected */
    beta1 = Beta[3][1][CInt(Z)];
    beta2 = Beta[3][2][CInt(I_Z_CN-Z)];
    if(Csng(I_Z_CN) * 0.5 < ZC_Mode_2)
    {
      Z1 = Z_equi(I_Z_CN,I, I_A_CN-I, beta1, beta2, DNECK,2);
      Z1 = Z1 + POLARadd;
      // Polarization caused by N=50 shell (assumption)
      Z1 = Z1 - 0.55 * Bell(Csng(I)-Z1,45,49.5);
    }
    else
      Z1 = Z_equi(I_Z_CN,I, I_A_CN-I, beta1, beta2, DNECK,2);      
    Zmean[2][1][I] = Z1;
    Zshift[2][1][I] = Z1 - ZUCD;        // neutron deficient
    Zmean[2][2][I_A_CN - I] = Csng(I_Z_CN) - Z1;  
    Zshift[2][2][I_A_CN - I] = ZUCD - Z1;  // neutron rich at shell
    ZshiftOriginal[2][1][I]=Zshift[2][1][I];     
    ZshiftOriginal[2][2][I_A_CN - I]=Zshift[2][2][I_A_CN - I];     
  }
  /* Mode 3 */
  for(int I = 10;I<=(I_A_CN - 10);I++)
  {
    ZUCD = Csng(I) / Csng(I_A_CN) * Csng(I_Z_CN);
    Z = ZUCD; /* Charge polarisation is here neglected */
    beta1 = Beta[4][1][CInt(Z)];
    beta2 = Beta[4][2][I_Z_CN-CInt(Z)];
    Z1 = Z_equi(I_Z_CN,I, I_A_CN - I, beta1, beta2, DNECK,3);
    Z1 = Z1 + POLARadd;
    // Polarization caused by N=50 shell (assumption)
    Z1 = Z1 - 0.55 * Bell(Csng(I)-Z1,45,49.5);
    //           POLARadd+0.15,POLARfac)   // Stronger charge polarization in S3, heavy fragment  //!!! KHS
    //           POLARadd+0.4,POLARfac)   // Stronger charge polarization in S3, heavy fragment
    Zmean[3][1][I] = Z1;
    Zshift[3][1][I] = Z1 - ZUCD;
    Zmean[3][2][I_A_CN - I] = Csng(I_Z_CN) - Z1;
    Zshift[3][2][I_A_CN - I] = ZUCD - Z1;
    ZshiftOriginal[3][1][I]=Zshift[3][1][I];     
    ZshiftOriginal[3][2][I_A_CN - I]=Zshift[3][2][I_A_CN - I];     
  }
  /* Mode 4 (Charge polarization of heavy fragment assumed to be equal to light fragment in S2 in the actinides) */
  for(int I = 10;I<=(I_A_CN - 10);I++)   // Loop is over the "second" (heavy) fragment!
  {
    ZUCD = Csng(I) / Csng(I_A_CN) * Csng(I_Z_CN);
    Z = ZUCD; /* Charge polarisation is here neglected */
    beta1 = Beta[5][1][CInt(I_Z_CN-Z)];  // light fragment
    beta2 = Beta[5][2][CInt(Z)];          // heavy fragment      
    Z2 = Z_equi(I_Z_CN,I, I_A_CN-I, beta2, beta1, DNECK,4);
    Z2 = Z2 + POLARadd;
    // Polarization caused by N=50 shell (assumption)
    //  Z2 = Z2 - 0.55 * Bell(Csng(I)-Z1,45,49.5)
    Zmean[4][2][I] = Z2;
    //Print "Z2UCD,Z2",ZUCD,Z2-POLARadd,Z2      
    Zshift[4][2][I] = Z2 - ZUCD;        // neutron deficient
    Zmean[4][1][I_A_CN - I] = Csng(I_Z_CN) - Z2;  
    Zshift[4][1][I_A_CN - I] = ZUCD - Z2;  // neutron rich at shell
    ZshiftOriginal[4][2][I]=Zshift[4][2][I];     
    ZshiftOriginal[4][1][I_A_CN - I]=Zshift[4][1][I_A_CN - I];     
  }  
}

void GEF::RelationsZ_A_FM(void)
{
  // Mode 0
  RZpol = 0;
  for(int I = 1;I<= 3;I++)
  {
    RA = (ZC_Mode_0 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
    RZpol = Zshift[0][2][CInt(RA)];  // heavy fragment
  }
  AC_Mode_0 = (ZC_Mode_0 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN); /// mean position in mass ///
  NC_Mode_0 = AC_Mode_0 - ZC_Mode_0;
  // Mode 1
  RZpol = 0;
  for(int I = 1;I<= 3;I++)
  {
    RA = (ZC_Mode_1 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
    RZpol = Zshift[1][2][CInt(RA)];  // heavy fragment
  }
  AC_Mode_1 = (ZC_Mode_1 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
  NC_Mode_1 = AC_Mode_1 - ZC_Mode_1;
  // Mode 2
  RZpol = 0;
  for(int I = 1;I<= 3;I++)
  {
    RA = (ZC_Mode_2 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
    RZpol = Zshift[2][2][CInt(RA)]; // heavy fragment
  }
  AC_Mode_2 = (ZC_Mode_2 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
  NC_Mode_2 = AC_Mode_2 - ZC_Mode_2;
  // Mode 3    
  RZpol = 0;
  for(int I = 1;I<= 3;I++)
  {
    RA = (ZC_Mode_3 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
    RZpol = Zshift[3][2][CInt(RA)];  // heavy fragment
  }
  AC_Mode_3 = (ZC_Mode_3 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
  NC_Mode_3 = AC_Mode_3 - ZC_Mode_3;
  // Mode 4
  RZpol = 0;
  for(int I = 1;I<= 3;I++)
  {
    RA = (ZC_Mode_4 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);
    RZpol = Zshift[4][1][CInt(RA)]; // light fragment
  }
  AC_Mode_4 = (ZC_Mode_4 - RZpol) * Csng(I_A_CN) / Csng(I_Z_CN);  // light fragment
  NC_Mode_4 = AC_Mode_4 - ZC_Mode_4;   // light fragment
}

void GEF::FissionBarriers(void)
{
  /* Fission barriers -> global parameters */
  B_F = BFTF(Csng(I_Z_CN),Csng(I_A_CN),1);   
  B_F_ld = BFTF(Csng(I_Z_CN),Csng(I_A_CN),0);
  E_B = BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);   
  E_B_ld = BFTFB(Csng(I_Z_CN),Csng(I_A_CN),0);
}

void GEF::PolarizationStiffness(void)
{
  RZ = Csng(I_Z_CN) * 0.5;
  RA = Csng(I_A_CN) * 0.5;
  beta1 = Beta[1][1][CInt(RZ)];
  beta2 = Beta[1][2][CInt(RZ)];
  R_Pol_Curv_S0 = ( LyMass( RZ - 1.E0, RA, beta1 ) + LyMass( RZ + 1.0E0, RA, beta2 ) + 
      LyMass( RZ + 1.0E0, RA, beta1 ) + LyMass( RZ - 1.0E0, RA, beta2 ) + 
      ECOUL( RZ - 1.0E0, RA, beta1, RZ + 1.0E0, RA, beta2, DNECK) + 
      ECOUL( RZ + 1.0E0, RA, beta1,RZ - 1.0E0, RA, beta2, DNECK) - 
      2.0E0*ECOUL( RZ, RA, beta1, RZ, RA, beta2, DNECK) - 
      2.0E0*LyMass( RZ, RA, beta1 ) -  2.0E0*LyMass( RZ, RA, beta2) ) * 0.5;

  P_POL_CURV_S0= R_Pol_Curv_S0;

  // Assumption: stiffenss is dominated (and thus well represented) by the macroscopic potential
  R_Pol_Curv_S1 = R_Pol_Curv_S0;
  R_Pol_Curv_S2 = R_Pol_Curv_S0;
  R_Pol_Curv_S3 = R_Pol_Curv_S0;
  R_Pol_Curv_S4 = R_Pol_Curv_S0;
}









float GEF::Getyield(float E_rel,float E_ref, float T_low, float T_high)
{
  // Erel: Energy relative to the barrier 
  // T_low: Effective temperature below barrier 
  // T_high: Effective temperature above barrier/
  float Exp1;
  float Yield;

  Exp1 = E_rel/T_low - E_ref/0.4;   //energy far below barrier   // Subtraction of E_ref/0.4 to prevent numerical problems.
  if(Exp1<-50)
    Yield = 0.;
  else
    Yield = exp(E_rel / T_high - E_ref/0.4) * 1. / (1. + exp(-E_rel/ (T_high*T_low/(T_high-T_low))));

  return Max(Yield,0.0);
}

float GEF::F1(float Z_S_A)
{
  // Fit to the lower part of the data 
  float Result;
  Result = exp(-9.05 + 4.58 *log(Z_S_A/2.3));
  return Result;
}

float GEF::F2(float Z_S_A)
{
  // Fit to the upper part of the data
  float Result;
  Result = exp(12.08 - 3.27 * log(Z_S_A/2.3));
  return Result;
}

float GEF::Masscurv(float Z, float A, float RL, float kappa)
{
  //  Fit to  Data of Fig. 7 of                                            
  //  "Shell effect in the symmetric-modal fission of pre-actinide nuclei" 
  //  S. I. Mulgin, K.-H. Schmidt, A. Grewe, S. V. Zhdanov               
  //  Nucl. Phys. A 640 (1998) 375 (From fit of the width of the mass distributions.) 
  float  RI, Result1, Result2, Result; 
  float Z_square_over_A;
  float ZsqrA;
  float c_rot = 600.0;

  Z_square_over_A = Z*Z/A;
  RI = (A - 2*Z)/A;
  ZsqrA = Z_square_over_A * (1. - kappa * pow(RI,2)) /(1. - kappa * pow((226. - 2.*91.)/226.,2)) + c_rot * pow(RL,2) / pow(A,7.0/3.0);// Hasse & Myers
  Result1 = F1(ZsqrA);
  Result2 = F2(ZsqrA);
  Result = Min(Result1,Result2);
  return Result;
}

float GEF::Masscurv1(float Z, float A, float RL, float kappa)
{
  //  Fit to  Data of Fig. 7 of                                            
  //  "Shell effect in the symmetric-modal fission of pre-actinide nuclei"  
  //  S. I. Mulgin, K.-H. Schmidt, A. Grewe, S. V. Zhdanov                 
  //  Nucl. Phys. A 640 (1998) 375 (The left part assumed to be valid for the yields of the fission channels.)
  float RI,Result1, Result2, Result ;
  float Z_square_over_A;
  float ZsqrA;
  float c_rot = 600.0;

  Z_square_over_A = Z*Z/A;
  RI = (A - 2*Z)/A;
  ZsqrA = Z_square_over_A * (1. - kappa * pow(RI,2)) / (1. - kappa * pow((226. - 2.*91.)/226.,2))+ c_rot * pow(RL,2) / pow(A,7.0/3.0);  // Hasse & Myers

  if(ZsqrA < 36.0) // adjusted to Y(S2) in light nuclei (80<Z<92)
    ZsqrA = ZsqrA + 0.9 * (36.0 - ZsqrA)  ;
  Result1 = F1(ZsqrA);
  Result2 = F2(ZsqrA);

  return Result1;
}

float GEF::De_Saddle_Scission(float Z_square_over_Athird,float ESHIFTSASCI)
{
  // Energy release between saddle and scission 
  // M. Asghar, R. W. Hasse, J. Physique C 6 (1984) 455 
  float Result;
  Result = (31. - 11.) / (1550. - 1300.) * (Z_square_over_Athird - 1300. + ESHIFTSASCI) + 11.;
  // This formula with ESHIFTSASCI = 0 is the parameterisation of the results
  // of Ashgar and Hasse, JPC 6 (1984) 455, see 
  // F. Rejmund, A. V. Ignatyuk, A. R. Junghans, K.-H. Schmidt
  // Nucl. Phys. A 678 (2000) 215  
  Result = Max(Result,0.0);
  return  Result;
}

float GEF::TEgidy(float A,float DU, float Fred)
{
  // Temperature parameter of the constant-temperature formula for the nuclear level density.
  // Input parameters: A = Mass number of nucleus, DU = Shell effect (corrected for pairing:P=0 for odd-A nuclei)
  //  From "Correlations between the nuclear level density parameters"
  //    Dorel Bucurescu, Till von Egidy
  //    Phys. Rev. C 72 (2005) 067304    and
  //        "Systematics of nuclear level density parameters"
  //     Dorel Bucurescu, Till von Egidy
  //     J. Phys. G: Nucl. Part. Phys. 31 (2005) S1675 and
  //        "Systematics of nuclear level density parameters"
  //     Till von Egidy, Dorel Bucurescu
  //     Phys. Rev. C 72 (2005) 044311 
  float Temp_smooth,Temp,T_Fac;
  Temp_smooth = 1.0 / (0.0570 * pow(A,0.6666667));
  Temp = 1.0 / ( (0.0570 + 0.00193*DU) * pow(A,0.6666667));  // from  PRC 80 (2009) 054310 
  T_Fac = Temp / Temp_smooth;
  Temp = Temp * Fred;  // (For influence of deformation) 
  return Temp;
}

float GEF::TRusanov(float E, float A)
{
  // Fermi-gas level density, parameterisation of Rusanov et al. 
  return sqrt(E / (0.094 * A) );
}

float GEF::LyMass(float Z, float A, float beta)
{ 
  // liquid-drop mass, Myers & Swiatecki, Lysekil, 1967  
  // pure liquid drop, without pairing and shell effects 
  // On input:    Z     nuclear charge of nucleus        
  //              N     number of neutrons in nucleus    
  //              beta  deformation of nucleus           
  // On output:   binding energy of nucleus              

  float N;
  float alpha;
  float XCOM,XVS,XE,EL;

  N = A - Z;
  alpha = sqrt(5./(4.*pi)) * beta;
  XCOM = 1. - 1.7826 * pow((A - 2.*Z)/A,2);
  // factor for asymmetry dependence of surface and volume term 
  XVS = - XCOM * (15.4941*A - 17.9439*pow(A,2./3.)*(1.+0.4*pow(alpha,2)));
  // sum of volume and surface energy 
  XE = Z*Z * (0.7053/pow(A,1./3.)*(1.-0.2*pow(alpha,2))- 1.1529/A);
  EL = XVS + XE;
  // EL = EL + LyPair(Z,A); 
  return EL;  
}

float GEF::LyPair(int Z, int A)
{
  // Calculates pairing energy 
  // odd-odd nucleus:   Lypair = 0 
  // even-odd nucleus:  Lypair = -12/sqr(A) 
  // even-even nucleus: Lypair = -2*12/sqr(A) 
  float E_PAIR;
  E_PAIR = - 12. / sqrt(Csng(A)) * Csng( ( (Z+1)%2 + (A-Z+1)%2) );
  return E_PAIR;
}

float GEF::TFPair(int Z, int A)
{
  // Pairing energy from Thomas-Fermi model of Myers and Swiatecki 
  // Shifted that TFPair is zero for odd-odd nuclei 
  int N;
  float E_Pair;
  N = A - Z;
  if(Z%2==0 && N%2 == 0)// even-even 
    E_Pair = - 4.8 / pow(Z,0.333333) - 4.8 / pow(N,0.333333) + 6.6 / pow(A,0.666666);
  if(Z%2==0 && N%2 == 1) // even Z, odd N 
    E_Pair = - 4.8 / pow(Z,0.333333) + 6.6 / pow(A,0.666666);
  if(Z%2==1 && N%2 == 0) // odd Z, even N 
    E_Pair = - 4.8 / pow(N,0.333333) + 6.6 / pow(A,0.666666);
  if(Z%2==1 && N%2 == 1) // odd N, odd N 
    E_Pair = 0.0;

  return E_Pair;
}

float GEF::Pmass(float Z, float A, float beta)
{
  // Liquid-drop model of Pearson, 2001 
  float N,EA,BE;
  float avol = -15.65;
  float asf = 17.63;
  float r0 = 1.233;
  float asym = 27.72;
  float ass = -25.60;
  float alpha;
  N = A - Z;
  alpha = sqrt(5./(4.*pi)) * beta;
  EA = avol + asf * pow(A,-0.333333)*(1.+0.4*pow(alpha,2)) + 0.6 * 1.44 * Z*Z / (pow(A,1.333333) * r0 )*(1.-0.2*pow(alpha,2)) + (asym + ass * pow(A,2-0.333333)) * pow(N-Z,2) / pow(A,2);
  BE = EA * A;
  return BE;
}

float GEF::FEDEFOP(float Z, float A, float beta)
{
  // According to liquid-drop model of Pearson 2001 
  float asf = 17.63;
  float r0 = 1.233;
  float N,Alpha;     
  N = A - Z;
  Alpha = sqrt(5./(4.*pi)) * beta;
  return asf * pow(A,0.666667)*(0.4*pow(Alpha,2)) - 0.6 * 1.44 * Z*Z / (pow(A,0.333333) * r0 )*(0.2*pow(Alpha,2));
}

float GEF::FEDEFOLys(float Z, float A, float beta)
{
  return LyMass(Z,A,beta) - LyMass(Z,A,0.0);
}

float GEF::LDMass(float Z, float A, float beta)
{
  float N,BEtab;
  N = A - Z;
  BEtab = BEldmTF[CInt(N)][CInt(Z)] + 2.0 * 12.0 / sqrt(Csng(A)) - 0.00001433*pow(Z,2.39);
  // The values in BEtab are the negative binding energies! 
  // Pairing in Thomas Fermi masses is zero for Z,N even !
  if(BEtab == 0.0)
  {
    BEtab = LyMass(Z,A,0.0);
    cout<< "Warning: Binding energy of Z="<<Z<<", A="<<A<<" not in mass table, replaced by LYMASS"<<endl;
    cout<< "I_Mode = "<<I_Mode<<endl;               
  } 
  return BEtab + FEDEFOLys(Z,A,beta);
}

float GEF::AME2012(int IZ, int IA)
{
  // Masses from the 2003 mass evaluation, complemented by TF masses
  // and Lysekil masses.
  float BEexpval;
  float Z,A,N;
  int INeu;
  INeu = IA - IZ;
  A = Csng(IA);
  Z = Csng(IZ);
  N = A - Z;
  BEexpval = BEexp[INeu][IZ]; 
  if(BEexpval > -1.e10)
    return BEexpval;
  else
    return LDMass(Z,A,0.0) + U_SHELL(IZ,IA) + LyPair(IZ,IA);
}

float GEF::U_SHELL(int Z, int A)
{
  int N;
  float Res;
  N = A - Z;
  Res = ShellMO[N][Z];
  if(Res>0.0)
    Res = 0.3 * Res;     // KHS (12. Feb. 2012)
  // The positive shell effects for deformed nuclei seem to be too positive
  // This gives too many high-energetic prompt neutrons.
  return Res;
}

float GEF::U_SHELL_exp(int IZ, int IA)
{
  float Res;
  float Z,A;
  Z = Csng(IZ);
  A = Csng(IA);
  //   Res = 2.0 * ( AME2012(IZ,IA) - Lypair(IZ,IA) - LDMass(Z,A,0.0) )
  //          - 0.25 * ( AME2012(IZ,IA-1) - Lypair(IZ,IA-1) - LDMass(Z,A-1.0,0.0) )
  //          - 0.25 * ( AME2012(IZ,IA+1) - Lypair(IZ,IA+1) - LDMass(Z,A+1.0,0.0) )
  //          - 0.25 * ( AME2012(IZ+1,IA+1) - Lypair(IZ+1,IA+1) - LDMass(Z+1.0,A+1.0,0.0) )
  //          - 0.25 * ( AME2012(IZ-1,IA-1) - Lypair(IZ-1,IA-1) - LDMass(Z-1.0,A-1.0,0.0) )
  Res = 0.5 * ( AME2012(IZ,IA) - LyPair(IZ,IA) - LDMass(Z,A,0.0) )
    + 0.125 * ( AME2012(IZ,IA-1) - LyPair(IZ,IA-1) - LDMass(Z,A-1.0,0.0) )
    + 0.125 * ( AME2012(IZ,IA+1) - LyPair(IZ,IA+1) - LDMass(Z,A+1.0,0.0) )
    + 0.125 * ( AME2012(IZ+1,IA+1) - LyPair(IZ+1,IA+1) - LDMass(Z+1.0,A+1.0,0.0) )
    + 0.125 * ( AME2012(IZ-1,IA-1) - LyPair(IZ-1,IA-1) - LDMass(Z-1.0,A-1.0,0.0) );
  return Res;
}

float GEF::U_SHELL_EO_exp(int IZ, int IA)
{
  // Returns experimental shell and even-odd staggering,
  // just the difference of experimental and macroscopic mass.
  float Res;
  float Z,A;
  Z = Csng(IZ);
  A = Csng(IA);
  Res = AME2012(IZ,IA) - LDMass(Z,A,0.0); 
  return Res;             
}

float GEF::U_MASS(float Z, float A)
{
  // LD + congruence energy + shell (no pairing) 
  float BE;
  if((Z<0) |(A<0))
    cout<<"U_Mass: Z, A: "<<Z<<" "<<A<<endl;

  BE = LDMass(Z,A,0.0)  + U_SHELL(CInt(Z),CInt(A));
  //    BE = AME2012(Cint(Z),Cint(A)) - Lypair(Z,A)
  //    BE = Lymass(Z,A,0.0) + U_Shell(CInt(Z),CInt(A))     
  //    BE = Lymass(Z,A,0.0)  
  return BE;
}

float GEF::ECOUL(float Z1, float A1, float beta1, float Z2, float A2, float beta2, float d)
{
  // Coulomb potential between two nuclei                  
  // surfaces are in a distance of d                       
  // in a tip to tip configuration                          
  // approximate formulation                               
  // On input: Z1      nuclear charge of first nucleus     
  //           A1      mass number of irst nucleus   
  //           beta1   deformation of first nucleus         
  //           Z2      nuclear charge of second nucleus     
  //           A2      mass number of second nucleus  
  //           beta2   deformation of second nucleus       
  //           d       distance of surfaces of the nuclei  

  float N1,N2,REcoul;
  float dtot;
  float r0 = 1.16;

  N1 = A1 - Z1;
  N2 = A2 - Z2;
  dtot = r0 *( pow(Z1+N1,0.3333333) * (1.+0.6666667*beta1)
      + pow(Z2+N2,0.3333333) * (1.+0.6666667*beta2) ) + d;
  REcoul = Z1 * Z2 * 1.44 / dtot;
  return REcoul;
}

float GEF::beta_light(int Z, float betaL0, float betaL1)
{
  // Deformation of light fission fragment for S1 and S2 
  // Systematic correlation Z vs. beta for deformed shells 
  // Z of fission fragment 
  float beta;
  beta = (Z - betaL0) * betaL1/20.; 
  return beta;
}

float GEF::beta_heavy(int Z, float betaH0, float betaH1)
{
  // Deformation of heavy fission fragment for S2 
  // Systematic correlation Z vs. beta for deformed shells 
  // Z of fission fragment 
  float beta;
  beta = (Z - betaH0) * betaH1/20.;
  return beta;
}

float GEF::Z_equi(int ZCN, int A1, int A2, float beta1, float beta2, float d, int Imode)
{
  // Determines the minimum potential of the scission-point configuration
  //   represented by two deformed nuclei divided by a tip distance d.
  //   A1, A2, beta1, beta2, d are fixed, Z1 is searched for and returned on output.  
  // ZCN: Z of fissioning nucleus 
  // A1: A of first fission fragment 
  // A2: A of second fission fragment 
  // beta1: deformation of first fission fragment 
  // beta2: deformation of second fission fragment 
  // d: tip distance 

  float RZ_equi;
  float RA1,RA2,RZCN,RACN;
  float Z1UCD,Z2UCD;
  float re1,re2,re3,eps1,eps2,DZ_Pol; // help variables 
  RA1 = Csng(A1);
  RA2 = Csng(A2);
  RZCN = Csng(ZCN);       
  RACN = RA1 + RA2;
  Z1UCD = RA1 / (RA1 + RA2) * RZCN;
  Z2UCD = RZCN - Z1UCD;
  re1 = LyMass( Z1UCD-1., RA1, beta1 ) +
    LyMass( Z2UCD+1., RA2, beta2 ) +
    ECOUL( Z1UCD-1., RA1, beta1, Z2UCD+1., RA2, beta2, d );
  re2 = LyMass( Z1UCD, RA1, beta1) +
    LyMass( Z2UCD, RA2, beta2) +
    ECOUL( Z1UCD, RA1, beta1, Z2UCD, RA2, beta2, d );
  re3 = LyMass( Z1UCD+1., RA1, beta1 ) +
    LyMass( Z2UCD-1., RA2, beta2 ) +
    ECOUL( Z1UCD+1., RA1, beta1, Z2UCD-1., RA2, beta2, d );
  eps2 = ( re1 - 2.*re2 + re3 ) / 2.;
  eps1 = ( re3 - re1 ) / 2.;
  DZ_Pol = -eps1 / ( 2. * eps2 );
  if((DZ_Pol>2) | (DZ_Pol<-2))
    DZ_Pol = 0.;
  RZ_equi = Z1UCD + DZ_Pol;   
  return RZ_equi;
}

void GEF::Beta_opt_light(float A1, float A2, float Z1, float Z2, float d, float beta2_imposed, float &beta1_opt)
{
  /// Determines the optimum deformation of the light fragment when the deformation of the heavy fragment is imposed.
  float beta1,dbeta1,beta1_prev,beta1_next;
  float Uguess,Uplus,Uminus,Uprev,Unext;
  int I;
  beta1 = 0.5;
  dbeta1 = 0.01;
  Uguess = LyMass(Z1, A1, beta1) + 
    LyMass(Z2, A2, beta2_imposed) + 
    ECOUL(Z1, A1, beta1, Z2, A2, beta2_imposed, d);
  Uplus  = LyMass(Z1, A1, beta1 + dbeta1) + 
    LyMass(Z2, A2, beta2_imposed) + 
    ECOUL(Z1, A1, beta1 + dbeta1, Z2, A2, beta2_imposed, d);
  Uminus = LyMass(Z1, A1, beta1 - dbeta1) + 
    LyMass(Z2, A2, beta2_imposed) + 
    ECOUL(Z1, A1, beta1 - dbeta1, Z2, A2, beta2_imposed, d);
  if(Uplus>Uguess && Uminus>Uguess)
    beta1_opt = beta1;
  else
  {
    if(Uplus < Uguess)
      dbeta1 = 0.01;
    if(Uminus < Uguess)
      dbeta1 = -0.01;
    Unext = Uguess;
    beta1_next = beta1;
    for(int i=1;i<=10000;i++)
    {
      beta1_prev = beta1_next;
      Uprev = Unext;
      beta1_next = beta1_prev + dbeta1;
      Unext = LyMass(Z1, A1, beta1_next) + 
        LyMass(Z2, A2, beta2_imposed) + 
        ECOUL(Z1, A1, beta1_next, Z2, A2, beta2_imposed, d);
      if(Unext>=Uprev)
        break;
      if(i==1000)
        cout<< "Loop overflow in Beta_opt_light"<<endl;
    }
    beta1_opt = beta1_prev;
  }
}

void GEF::Beta_Equi(float A1, float A2, float Z1, float Z2, float d, float beta1prev, float beta2prev, float &beta1opt, float &beta2opt)
{
  // Determines the minimum potential of the scission-point configuration
  // represented by two deformed nuclei, divided by a tip distance d.
  // A1, A2, Z1, Z2, d are fixed, beta1 and beta2 are searched for and returned on output

  int B_analytical = 0;
  // Switch to use the analytical approximation 
  // that replaces the long numerical calculation.
  float x,y,xcoul;
  float xcoul236U = 1369.64;
  float beta1,beta2;
  float U,Uprev,Ulast,Ubest,Uopt;
  float sbeta1 = 0;
  float sbeta2 = 0;
  int N,N1,N2;
  int Nopt = 0;
  float eps = 5.e-4;

  if(B_analytical == 0) // Numerical algorithm
  {
    beta1 = beta1prev;
    beta2 = beta2prev;
    Uprev = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECOUL(Z1,A1,beta1,Z2,A2,beta2,d);
    Uopt = Uprev;

    // Test slope of variation of U 
    beta1 = beta1prev + eps;

    U = 1.e30;
    beta2 = beta2prev;
    for(int i=1;i<=CInt(beta2prev/eps);i++)
    {
      beta2 = beta2 - eps;
      Ulast = U;
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECOUL(Z1,A1,beta1,Z2,A2,beta2,d);
      if(U>Ulast)
        break;
      else
        Ubest = U;
    }
    if(Ubest < Uopt)
    {
      Uopt = Ubest;
      sbeta1 = eps;
      sbeta2 = -eps;
    }

    U = 1.e30;
    beta2 = beta2prev;
    for(int i=1;i<=CInt((1 - beta2prev)/eps);i++)
    {
      beta2 = beta2 + eps;
      Ulast = U;
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECOUL(Z1,A1,beta1,Z2,A2,beta2,d);
      if(U > Ulast)
        break;
      else
        Ubest = U;
    }
    if(Ubest < Uopt)
    {
      Uopt = Ubest;
      sbeta1 = eps;
      sbeta2 = eps;
    }

    beta1 = beta1prev - eps;
    U = 1.e30;
    beta2 = beta2prev;
    for(int i=1;i<=CInt(beta2prev/eps);i++)
    {
      beta2 = beta2 - eps;
      Ulast = U;
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECOUL(Z1,A1,beta1,Z2,A2,beta2,d);
      if(U > Ulast)
        break;
      else
        Ubest = U;
    }
    if(Ubest < Uopt)
    {
      Uopt = Ubest;
      sbeta1 = -eps;
      sbeta2 = -eps;
    }

    U = 1.e30;
    beta2 = beta2prev;
    for(int i=1;i<=CInt((1-beta2prev)/eps);i++)
    {
      beta2 = beta2 + eps;
      Ulast = U;
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECOUL(Z1,A1,beta1,Z2,A2,beta2,d);
      if(U > Ulast)
        break;
      else
        Ubest = U;
    } 
    if(Ubest < Uopt)
    {
      Uopt = Ubest;
      sbeta1 = -eps;
      sbeta2 = eps;
    }

    Ubest = LyMass(Z1,A1,beta1prev) + LyMass(Z2,A2,beta2prev)
      + ECOUL(Z1,A1,beta1prev,Z2,A2,beta2prev,d);
    U = LyMass(Z1,A1,beta1prev+Csng(sbeta1)) +
      LyMass(Z2,A2,beta2prev+Csng(sbeta2)) +
      ECOUL(Z1,A1,beta1prev+sbeta1,Z2,A2,beta2prev+Csng(sbeta2),d);

    for(N=1;N<=1000;N++)
    {
      for(N1 =1;N1<=N;N1++)
      {
        N2 = N-N1;
        beta1 = beta1prev + sbeta1*N1;
        beta2 = beta2prev + sbeta2*N2;
        U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECOUL(Z1,A1,beta1,Z2,A2,beta2,d);
        if(U < Ubest)
        {
          Ubest = U;
          beta1opt = beta1;
          beta2opt = beta2;
          Nopt = N;
        }
      }
      if(N-Nopt > 2)
        break; 
    }
    if(N>998)
      cout<<"Beta_Equi not converged: "<<Z1<<" "<<N<<endl;
  }
  else // Analytical approximation
  {     // Must be adapted if the relevant parameters of GEF are modified!
    xcoul = pow(Z1 + Z2,2) / pow(A1 + A2,1.0/3.0);
    x = pow(Z1 / (Z1 + Z2),xcoul/xcoul236U);
    y = 1.2512e-4 + 0.00122851*x - 0.00267707*pow(x,2) + 0.00372901*pow(x,3) - 0.00219903*pow(x,4);       
    beta1opt = y * xcoul ;
    x = pow(Z2 / (Z1 + Z2),xcoul/xcoul236U);
    y = 1.2512e-4 + 0.00122851*x - 0.00267707*pow(x,2) + 0.00372901*pow(x,3) - 0.00219903*pow(x,4);       
    beta2opt = y * xcoul ;  
  }   
}

void GEF::Eva(int Ilh, float Z_CN, float A_CN, float E_INIT, float T, float J_frag, float &Z_RES, float &A_RES, float &E_FINAL, float *Array_En, float *Array_Tn, float *Array_Eg0)
{
  // Z_CN,A_CN,E_init       Parameters of initial nucleus 
  // Z_RES,A_RES,E_FINAL    after evaporation 
  // T temperature coefficient of level density (not used) 
  // Array_En       kinetic energies of neutrons 
  // Array_Tn       neutron emission time after scission 
  // Array_Eg0      energies of statistical gammas 
  static float E_MIN = 0;   // Final energy for evaporation chain 
  float SN,SNeff,SNmean,SNexp,SNld;        // Neutron separation energy 
  float RShell;
  int ITry,Ifold;
  float Ai,Af,Zi,Zf,Ni,Nf,Ei,Ef;
  float bshift,bshiftm;
  float Tm, Td, Tf, Tmean, TCT;
  float Gamma_n,Gamma_g,Pgamma;
  float rho,rhom;
  float alev;
  float g_koeff;
  float E_kin;
  float Tn,Tn_acc;
  float J_crit;  // critical angular momentum (disappearance of pairing)
  float J_crit_shell;
  float Fred;   // reduction of pairing gap by ang. momentum
  float Fred_shell;  // reduction of shell
  int In_gamma;  // counts gammas for Array_Eg0()
  bool B_Expmass = 1;  // use mass model or empirical masses

  E_MIN = E_FINAL;
  In_gamma = 0;       
  Tn_acc = 0  ;
  Ifold = 1;
  Ai = A_CN;
  Zi = Z_CN;
  Ei = E_INIT;
  J_crit = 12. * sqrt(Ai/100.0); // L. G. Moretto, Nucl. Phys. A 185 (1972) 145
  J_crit_shell = 30.;

  if( J_frag >= J_crit)
    Fred = 0.;
  else
    Fred = 0.65 * sqrt(1.0 - J_frag/J_crit); // L. G. Moretto, Nucl. Phys. A 185 (1972) 145

  if(J_frag >= J_crit_shell) 
    Fred_shell = 0.;
  else
    Fred_shell = 0.9 - J_frag/J_crit_shell;

  // Shell effects included 
  if(B_Expmass == 0)
  {
    SN = U_MASS(Zi,Ai-1.) + Fred * LyPair(Zi,Ai-1.) - (U_MASS(Zi,Ai) + Fred * LyPair(Zi,Ai)) ;
    SNeff = U_MASS(Zi,Ai-1.) - U_MASS(Zi,Ai) ;   // with shells, without pairing
    //   SNmean = 0.5 * (U_MASS(Zi,Ai-2.) - U_MASS(Zi,Ai));
  }
  else
  {
    SNexp = ( AME2012(Zi,Ai-1) - AME2012(Zi,Ai) );  // empirical
    SNeff = U_MASS(Zi,Ai-1.) - U_MASS(Zi,Ai);    // with shells, without pairing
    SNld = LDMass(Zi,Ai-1.,0.0) - LDMass(Zi,Ai,0.0); // no shells, no pairing
    SN = Fred_shell * (Fred * SNexp + (1.0 - Fred) * SNeff) + (1.0 - Fred_shell) * SNld;      
  }  

  // Shell effects excluded 
  //   SN = LDMass(Zi,Ai-1.,0.0) + LyPair(Zi,Ai-1.) - (LDMass(Zi,Ai,0.0) + LyPair(Zi,Ai)) ;
  //   SNeff = LDMass(Zi,Ai-1.,0.0) - LDMass(Zi,Ai,0.0) ; 
  while(Ei-SN>E_MIN)
  {
    // Treat gamma competition 
    Tm = U_Temp(Zi,Ai,Ei,1,1,Tscale,Econd);       // emitting nucleus
    Td = U_Temp(Zi,Ai-1,Ei-SNeff,1,1,Tscale,Econd);
    if(Ilh > 0) // Emission from fragments
    {
      Gamma_g = 0.624 * pow(Ai,1.6) * pow(Tm,5);    // in meV (Ignatyuk, Bologna)
      Gamma_g = Gamma_g * 1.e-9;        // in MeV
    }
    else             // Emission between saddle and scission
      Gamma_g = 0;
    Tmean = (Tm + Td)/2;
    // Gamma_n = pow(Ai-1,0.66667) * 0.0302 * pow(Td,2) / exp(SNeff/Tmean)   // in MeV (Moretto)
    Gamma_n = pow(Ai-1,0.66667) * 0.13 * pow(Td,2) / exp(SNeff/Td);  // in MeV (Mor. PRC 54,3062)
    // Cut the distribution at the ground-state energy:    
    Gamma_n = Gamma_n - pow(Ai-1,0.66667) * 0.13 * pow(Td,2) / exp((SNeff+Ei)/Td);
    Tn = PExp(0.658 / Gamma_n);   // in units of 10^-21 s (hbar=0.658zs*MeV)
    // ? Tn = Tn * 2;      // Due to pre-exponential factor of Maxwellian, Tn is about 2 times to large   
    Tn_acc = Tn_acc + Tn;
    // Influence on lev.dens. of pairing at low E*
    // (Constant temperature assumed)
    if(Ei-SN < abs(Fred * LyPair(Zi,Ai-1)))  // rest energy below mean g.s. of odd-odd nuclei
    {
      if( (Mod(Zi,2) < 0.5) | (Mod(Ai-Zi-1,2) < 0.5))  // even Z or even N
        Gamma_n = exp(-12./sqrt(Ai)/Td) * Gamma_n;
      if((Mod(Zi,2) < 0.5) && (Mod(Ai-Zi-1,2) < 0.5))  //n Z and even N
        // For low level density below pairing gap in even-even nuclei
        Gamma_n = exp(-12./sqrt(Ai)/Td) * Gamma_n;
    }
    // Reduces the even-odd effect in neutron number 
    // due to low level density below the pairing gap 
    Pgamma = Gamma_g / (Gamma_g + Gamma_n);  
    if(rn->Rndm()<Pgamma) //gamma will be emitted
    {
      In_gamma = In_gamma + 1;
      //Scope //I don't konw what is doing this Scope (Diego)
      int N;
      float Eg;

      Eg = P_Egamma_high(Zi,Ai,Ei);      
      Array_Eg0[In_gamma] = Eg;
      // Accumulate E1 gammas  // I comment all this because I don't know what it means, by the moment (Diego)
      //N = CInt(Eg*1000);
      //if( N > 0)
      //{
      //  Ngtot = Ngtot + 1;
      //  if(Ilh == 1)
      //	Nglight = Nglight + 1;
      //  if(Ilh == 2)
      //	Ngheavy = Ngheavy + 1;
      //   Egtot1000 = Egtot1000 + EG*1000;
      //    if( N <= Ubound(Egamma))
      //     {
      //	  Egamma(N) = Egamma(N) + 1;
      //    #Ifdef B_EgammaA 
      //	  StoreEgammaA(N,A_CN);
      //   #Endif  
      //} 
      //}  
      Ei = Ei - Eg;        
      //End Scope
    }

    if( Ei-SN <= E_MIN)
    {
      Zf = Zi;
      Af = Ai;
      Ef = Ei;
      break;
    }

    ITry = 0;
    bool Too_Low=true;
    while (Too_Low)
    {
      Too_Low = false;
      ITry = ITry + 1;
      if(ITry < 99)
      {
        Td = U_Temp(Zi,Ai-1,Ei-SNeff,1,1,Tscale,Econd);
        // maximum residual energy of daughter nucleus (for En_kin = 0)
        E_kin = PMaxwellMod(Td,Ai-1);   // Maxwell, with partial 1/v behaviour 
        Tf = U_Temp(Zi,Ai-1,Ei-E_kin-SNeff,1,1,Tscale,Econd);
        // final energy of daughter nucleus with En_kin considered 

        // if(Ei-E_kin-2*Td > 10)   // In Fermi-gas regime
        if(Ei-E_kin-Td > 5)   // to avoid Tf at negative energies
        {
          //  if(rn->Rndm() > sqrt( exp(E_kin/Td)/ exp(E_kin/Tf) )) 
          //  if (rn->Rndm() > ( Exp(E_kin/Td)/ Exp(E_kin/Tf) )^0.33333 )  // last option
          //  if (rn->Rndm() > ( Exp(E_kin/Td)/ Exp(E_kin/Tf) )^0.25 )
          // Modified Maxwell that adapts to the Fermi-gas regime
          if( rn->Rndm() > pow( pow(exp(E_kin/Td),2)/ exp(E_kin/Tf),0.25) )
            Too_Low=true;
          // adjusted to data, justification not clear
        }
      }
      else
      {
        // E_kin too high after several attemps
        // no neutron emitted 
        Af = Ai;
        Zf = Zi;
        Ef = Ei;
        break;
      }

      if(E_kin > Ei-SN)
      {
        // E_kin from PMaxwell is not available 
        // Try again 
        Too_Low=true;
      }
    }
    Af = Ai - 1;
    Zf = Zi;
    Ef = Ei - E_kin - SN;

    // ANAL(EN,E_kin);
    //ANAL(ENM(I_MODE),E_kin); 

    // Shell effects included //  
    if( B_Expmass == 0 )
    {
      SN = (U_MASS(Zf,Af-1.) + Fred * LyPair(Zf,Af-1.)) - (U_MASS(Zf,Af) + Fred * LyPair(Zf,Af));
      SNeff = U_MASS(Zf,Af-1.) - U_MASS(Zf,Af);  
      //  SNmean = 0.5 * (U_MASS(Zf,Af-2.) - U_MASS(Zf,Af));
    }
    else
    {
      SNexp = ( AME2012(Zf,Af-1) - AME2012(Zf,Af) );  // empirical
      SNeff = U_MASS(Zf,Af-1.) - U_MASS(Zf,Af);    // with shells, without pairing
      SNld = LDMass(Zi,Ai-1.,0.0) - LDMass(Zi,Ai,0.0); // no shells, no pairing
      SN = Fred_shell * (Fred * SNexp + (1.0 - Fred) * SNeff) + (1.0 - Fred_shell) * SNld;      
    }

    // Shell effects excluded   
    //   SN = LDMass(Zf,Af-1.,0.0) + Lypair(Zf,Af-1.)  - (LDMass(Zf,Af,0.0) + Lypair(Zf,Af)) ;
    //   SNeff = LDMass(Zf,Af-1.,0.0) - LDMass(Zf,Af,0.0);          

    Ai = Af;
    Zi = Zf;
    Ei = Ef;

    Array_En[Ifold] = E_kin * (Af-1.) / Af;
    Array_Tn[Ifold] = Tn_acc;
    Ifold = Ifold + 1;  
  }


  A_RES = Af;
  Z_RES = Zf;
  E_FINAL = Max(Ef,0.0); 

}

float GEF::u_accel(float A1, float Z1, float A2, float Z2, float TKE , float E0, float Tn)
{
  // returns the velocity of the fragment 1 after time Tn
  //Acceleration of fission fragments by their Coulomb field    

  // natural constants
  static float e2 = 1.44;   // MeV
  static float u = 931.5;  // MeV / c^2 
  static float hbarc = 197;  // MeV fm

  // variables
  float Ared,d0,v0;
  float d, t, dt, a, v, v1, E;
  float vinf;         // relative velocity at infinity

  Ared = A1 * A2 / (A1 + A2);
  vinf = sqrt(TKE/Ared);         // sqrt(E/A), Ekin asymmptotic in MeV
  if((t > 100) | (TKE < E0))
    v = vinf;
  else
  {
    d0 = e2 * Z1 * Z2 / (TKE-E0); // fm   
    v0 = sqrt(E0);   // sqr(E/A), Ekin at scission in MeV

    d = d0;
    v = 0;
    dt = 0.01;
    t=0;
    while(t<=1)   // in  10^-21 s
    {
      if (t >= Tn)
        break;
      E = E0 + e2 * Z1 * Z2 * (1/d0 - 1/d);   // MeV
      v = sqrt(E/Ared);         // sqrt(E/A), E in MeV
      d = d + v * 14 * dt;             // in fm
      t += dt;
    }
    dt = 0.1;
    t=1.1;
    while(t<=10)  // in  10^-21 s
    {
      if(t >= Tn)
        break;
      E = E0 + e2 * Z1 * Z2 * (1/d0 - 1/d);   // MeV
      v = sqrt(E/Ared);         // sqrt(E/A), E in MeV
      d = d + v * 14 * dt;             // in fm
      t += dt;
    }
    dt = 1;
    t=11;
    while(t<100)  // in  10^-21 s
    {
      if(t >= Tn)
        break;
      E = E0 + e2 * Z1 * Z2 * (1/d0 - 1/d);   // MeV
      v = sqrt(E/Ared);         // sqrt(E/A), E in MeV
      d = d + v * 14 * dt;             // in fm
    }
  }
  v1 = v * A2/(A1+A2);
  return v1;  
}

float GEF::P_Egamma_low(float Zi, float Ai, float Ei)
{
  // Random function, returns gamma energy in MeV
  // For energies below Sn: no competition with neutrons
  float Rres;
  int N;
  float sigMax, Eg, Erest, rhorest, fg;
  float xran, yran;
  float Tm;
  float GammaExp;

  float betadef = 0;
  float gammadef = 0;

  float G0[4],E0[4];

  float alev;

  N = CInt(Ei*10 + 0.5); 
  if( N <= 0)
    N = 1;
  float sigma[N];  // sigma is not normalized!
  // (Normalization is done by Monte-Carlo procedure.)

  Tm = U_Temp(Zi,Ai,Ei,1,1,Tscale,Econd);

  betadef = DEFOtab[CInt(Ai-Zi)][CInt(Zi)];       // ground-state deformation
  gammadef = -120 * betadef + 47.4;  // A. Junghans

  // Print Tm
  // Tm = Tm / 2 ' + 0.1 * (1.0 / U_I_Shell(Zi,Ai) - 1)
  // Print Tm
  // sleep  
  // For eventual specific behaviour of gamma strength in magic nuclei
  //betadef = 0
  //gammadef = 0   
  if(betadef == 0 && gammadef == 0)
  {
    E0[2] = E0_GDR(Zi,Ai);
    G0[2] = Width_GDR(E0[2]);

    // Establish distribution
    sigMax = 0;
    Eg = 0.1;
    while(Eg <= Ei)
    {
      N = CInt(Eg*10);
      //       sigma(N) = exp(-Eg/Tm) * Eg^2; // for testing shape of gamma-strength function
      sigma[N] = exp(-Eg/Tm) *  pow(Eg,3) * (G0[2] * Eg) / ( pow(pow(Eg,2) - pow(E0[2],2),2) + pow(G0[2],2) * pow(E0[2],2) ) ; 
      // + 0.7 * G0(2) * 4 * pi^2 * Tm^2 / E0(2)^5 )  
      // last line: correction for low gamma energy (PRC 41 (190) 1941)
      if(sigma[N] > sigMax)
        sigMax = sigma[N];
      Eg += 0.1;
    }
    //Print "Tm,E0,G0)";Tm,E0(2),G0(2)       
    //For Eg = 0.1 To Ei Step 0.1
    //         N = CInt(Eg*10)
    //         Print Eg,sigma(N) / sigMax
    //Next Eg
  }
  else
  {
    sigMax = 0;
    for(int i=1;i<=3;i++)
    {
      E0[i] = E0_GDR(Zi,Ai) * Efac_def_GDR(betadef,gammadef,i);
      G0[i] = Width_GDR(E0[i]);
      // Establish distribution
      Eg = 0.1;
      while(Eg<=Ei)
      {
        N = CInt(Eg*10);
        // sigma(N) = sigma(N) + exp(-Eg/Tm) * Eg^3 *_
        //     (G0(I) * Eg) / ( (Eg^2 - E0(I)^2)^2 + G0(I)^2 * E0(I)^2 ) 
        sigma[N] = sigma[N] + exp(-Eg/Tm) * pow(Eg,3) *
          ( (G0[i] * Eg) / ( pow(pow(Eg,2) - pow(E0[i],2),2) + pow(G0[i],2) * pow(E0[i],2) )  
            /* +  0.001 * exp(-0.5)     */
            /* + 0.7 * G0[i] * 4 * pi^2 * Tm^2 / E0[i]^5*/ ) ;
        // exponential: M1 strength at low energy, PRL 111, 232504 (2013)    
        // last line: correction for low gamma energy (PRC 41 (1990) 1941)
        if (sigma[N] > sigMax)
          sigMax = sigma[N];
        Eg += 0.1;
      }
    }
  }  
  // Dice gamma energy from distribution
  bool diceagain_gamma_low = true;
  if(diceagain_gamma_low)
  {
    diceagain_gamma_low = false;
    xran = rn->Rndm() * Ei * 10;   // in units of 100 keV
    yran = rn->Rndm() * sigMax;
    if(yran > sigma[CInt(xran)])
      diceagain_gamma_low = true;
  }
  return xran/10;     // convert to MeV   
}

float GEF::P_Egamma_high(float Zi, float Ai, float Ei)
{
  // Random function, returns gamma energy in MeV 
  // From PRL 49 (1982) 434
  // For energies above Sn: competition with neutrons included
  int N;
  float sigMax, Eg;
  float xran, yran;
  float Tm;
  float G0,E0;

  N = CInt(Ei*10 + 0.5); 
  if(N <= 0)
    N = 1;
  float sigma[N];  // sigma is not normalized
  // (Normalization is done by Monte-Carlo procedure.)

  E0 = E0_GDR(Zi,Ai);
  G0 = Width_GDR(E0);
  Tm = U_Temp(Zi,Ai,Ei,1,1,Tscale,Econd);

  //Establish distribution
  sigMax = 0;
  Eg = 0.1;
  while (Eg<=Ei)
  {
    N = CInt(Eg*10.0);
    sigma[N] = pow(Eg,3) / pow(Tm,2) * exp(-Eg/Tm) * (G0 * Eg) / ( pow(pow(Eg,2) - pow(E0,2),2) + pow(G0,2) * pow(E0,2));
    if(sigma[N] > sigMax)
      sigMax = sigma[N];
    Eg += 0.1;
  }
  // Dice gamma energy from distribution
  bool diceagain_gamma_high=true;
  if(diceagain_gamma_high)
  {
    diceagain_gamma_high=false;
    xran = rn->Rndm() * Ei * 10;   // in units of 100 keV
    yran = rn->Rndm() * sigMax;
    if(yran > sigma[CInt(xran)])
      diceagain_gamma_high=true;
  }
  return xran/10;     // convert to MeV   
}

float GEF::U_Ired(float Z, float A)
{
  // Effective moment of inertia by pairing with correction for excitation energy
  float I_rigid_spher,IfragEff;

  I_rigid_spher = pow(1.16,2) * pow(A,1.6667) / 103.8415; 
  //   IfragEff = I_rigid_spher + 0.003 * A^(4.0/3.0) * U_shell(Cint(Z),Cint(A))
  //   IfragEff = I_rigid_spher + 0.005 * A^(4.0/3.0) * U_shell(Cint(Z),Cint(A))
  // reduction due to shell (Deleplanque et al. PRC 69 (2004) 044309)
  IfragEff = 0.45 * I_rigid_spher; // Effect of superfluidity 
  //   IfragEff = 0.65 * IfragEff   // Average effect of superfluidity and deformation 

  return IfragEff;     
}

float GEF::U_IredFF(float Z, float A)
{
  // Effective moment of inertia by pairing with correction for excitation energy
  // of final fission fragments
  return  U_Ired(Z,A) * U_I_Shell(Z,A);    
}

float GEF::U_I_Shell(float Z, float A)
{
  int N_shells[7];
  // Shell effect on the effective moment of inertia
  float dNmin, dZmin, dNsubmin;
  float Inv_add = 0;
  float I_inv_add_Z = 0;
  float I_inv_add_N = 0;
  float I_inv_add_Nsub = 0;
  N_shells[1] = 20;
  N_shells[2] = 28;
  N_shells[3] = 50;
  N_shells[4] = 82;
  N_shells[5] = 126;
  N_shells[6] = 56;
  dNmin = 100;
  dZmin = 100;
  dNsubmin = 100;
  for(int i=1;i<=5;i++)
    dZmin = Min(dZmin,abs(N_shells[i] - Z));
  for(int i=1;i<=5;i++)      
    dNmin = Min(dNmin,abs(N_shells[i] - (A-Z))); 
  dNsubmin = abs(N_shells[6] - (A-Z));

  // Effect of shells:
  if(dZmin < 10.0)
  {
    //    I_inv_add_Z = 0.33 * (6.0 * sqr(A/140.) - dZmin) * sqr(140./A)
    I_inv_add_Z = 0.33 * (6.0 * sqrt(A/140.0) - dZmin) * pow(140.0/A,1.5);
    // A^(-1/3) dependence: "A simple phenomenology for 2gamma+ states",
    // N. V. Zamfir, D. Bucurescu, R. F. Casten, M. Ivascu,
    // Phys. Lett. B 241 (1990) 463
    I_inv_add_Z = Max(I_inv_add_Z,0.0);
  }
  if(dNmin < 10.0)
  {
    //    I_inv_add_N = 0.42 * (8.0 * sqr(A/140.) - dNmin) * sqr(140./A)
    I_inv_add_N = 0.42 * (8.0 * sqrt(A/140.0) - dNmin) * pow(140.0/A,1.5);
    I_inv_add_N = Max(I_inv_add_N,0.0);
  }    
  if(dNsubmin < 6.0)
  {
    //   I_inv_add_Nsub = 1.7 * (4.0 - dNsubmin) * (1.0 - 0.32 * Abs(40.0-Z))
    I_inv_add_Nsub = 1.7 * (4.0 - dNsubmin) * (1.0 - 0.18 * abs(40.0-Z));
    // N = 56 subshell only around Z = 40
    I_inv_add_Nsub = Max(I_inv_add_Nsub,0.0);
  }     
  return 1.0 / (1.0 + Max(I_inv_add_N,I_inv_add_Nsub) + I_inv_add_Z);
  //Print "*",I_inv_add_Z, I_inv_add_N, I_inv_add_Nsub,1.0 / (1.0 + Max(I_inv_add_N,I_inv_add_Nsub) + I_inv_add_Z)    
}

float GEF::U_alev_ld(float Z, float A)
{
  //  return 0.073 * A + 0.095 * A^0.666667  //Ignatyuk (1970's)
  return 0.078 * A + 0.115 * pow(A,0.6666667);  // Ignatyuk (Bologna 2000) 
  //  return = 0.089 * A    // only volume term
}

double GEF::U_levdens(int Z, int A, float E, int Ishell, float Ipair, float Tscale, float Econd, float af_an)
{
  // Comment: The normalization of the CT level density to the FG level density
  //          reduces the jump in Pf around 16 MeV and leads to a more regular
  //          evolution of Pf with N_CN and Z_CN. The absolute values of Pf are not
  //          modified very much.
  float Etrans = 8.0;  // Transition from CT to Fermi gas
  double rho, rho1, rho2;
  if(E < Etrans)
  {
    rho1 = U_levdens_FG(Z,A,Etrans,Ishell,Ipair,Tscale,Econd,af_an);
    rho2 = U_levdens_Egidy(Z,A,Etrans,Ishell,Ipair,Tscale,Econd,af_an);
    rho = U_levdens_Egidy(Z,A,E,Ishell,Ipair,Tscale,Econd,af_an) * rho1 / rho2;
    //    rho = U_levdens_Egidy(Z,A,E,Ishell,Ipair,Tscale,Econd,af_an) 
  }
  else
  {
    //    rho1 = U_levdens_FG(Z,A,Etrans,Ishell,Ipair,Tscale,Econd,af_an)
    //     rho2 = U_levdens_Egidy(Z,A,Etrans,Ishell,Ipair,Tscale,Econd,af_an)
    //     rho = U_levdens_FG(Z,A,E,Ishell,Ipair,Tscale,Econd,af_an) * rho2 / rho1
    // normalized to FG formula at Etrans
    rho = U_levdens_FG(Z,A,E,Ishell,Ipair,Tscale,Econd,af_an); 
  }
  return rho;                           
}

double GEF::U_levdens_Egidy(int Z, int A, float E, int Ishell, int Ipair, float Tscale, float Econd, float af_an)
{
  // E may be given above ground state or above the macroscopic ground state         
  float Temp, DU, Ered, Rmicro, Rshell;    
  if( Ishell == 1)
  {
    DU = U_SHELL_exp(Z,A);   // For ground state (only shell effect)
    if(Z > 92)  // hypothetical deviation in the macroscopic masses
      DU = DU + 0.1 * (Z - 92);
  }
  else
    DU = 0;   // for barrier (no shell effect at all)    
  if(Ipair == 1)
  {
    Rmicro = - U_SHELL_EO_exp(Z,A);   // microscopic effects (shell and pairing)
    Rshell = - U_SHELL_exp(Z,A);   // only shell effect
    if(Z > 92)  // hypothetical deviation in the macroscopic masses
    {
      Rmicro = Rmicro - 0.1 * (Z - 92);
      Rshell = Rshell - 0.1 * (Z - 92);
    }
    Ered = E + Rshell - Rmicro  + 2.0 * 12.0 / sqrt(A);  // shift from odd-odd to even-even basis (exp)
    // compatibility with Egidy's formula is tested (Ered = E for even-even nuclei)
    //   Ered = E + Lypair(Z,A) + 2.0 * 12.0 / sqr(A)  // shift from odd-odd to even-even basis (schematic)
  }
  else
    Ered = E + 2.0 * 12.0 / sqrt(A);     // energy is given above the macr. ground state
  // If Z > 92 Then    // hypothetical shift of macr. masses -> lower shell effect
  //   DU = DU + 0.33 * (Z - 92)     
  // End If
  // DU = DU + Econd    
  Temp = 1.0 / ( (0.0570 + 0.00193*DU) * pow(A,0.6666667));  // from  PRC 80 (2009) 054310 
  return 1.0 / Temp * exp(Ered/Temp);                  
}

double GEF::U_levdens_FG(int Z, int A, float E, int Ishell, int Ipair, float Tscale, float Econd, float af_an)
{
  // This function calculates the level density in the FG regime,
  // taking a high-energy value as reference,
  // starting from the macroscopic masses with an imposed Econd.
  // This destroys all structural effects at high E* in a well controlled way.  
  // E on input is the excitation energy above the "real" ground state 
  float Ered, DU;  // energy above the fictive macroscopic ground state     
  float alev,Eshift_shell,Rmicro,Rshell; 
  double Rho1;
  float fgamma = 0.055;              
  float F_enhance = 1.3;  // normalization to CT formula
  //  float F_enhance = 0.18;  // normalization to CT formula
  if(Ishell == 1)
  {
    if(Ipair == 0)
    {
      Rshell = - U_SHELL_exp(Z,A);
      if(Z > 92)  // hypothetical deviation in the macroscopic masses
        Rshell = Rshell - 0.1 * (Z - 92);
      Rmicro = Rshell;
    }
    else
    {
      Rmicro = - U_SHELL_EO_exp(Z,A);   // microscopic effects (shell and pairing)
      Rshell = - U_SHELL_exp(Z,A);   // only shell effect
      if(Z > 92)   // hypothetical deviation in the macroscopic masses
      {
        Rshell = Rshell - 0.1 * (Z - 92);
        Rmicro = Rmicro - 0.1 * (Z - 92);
      }
      //     Eshift_shell = Rshell*exp(-fgamma * E)  // shifting the energy scale to include shell at low E*
      //  Print "*";Rmicro,Rshell;"*"  // The values are positive
    }
    DU = Rshell;      
    // If Z > 92 Then     // hypothetical shift of macr. masses -> lower shell effect
    //   DU = DU - 0.33 * (Z - 92)     
    // End If
    // DU = DU - Econd    
    Eshift_shell = DU * (1.0 - exp(-fgamma * E));   // remove shell effect at high energies
    //      Eshift_shell = Rshell * (1.0 - exp(-fgamma * E))   // remove shell effect at high energies
    Ered = E + Rshell - Rmicro - Econd;  // origin of the FG without pairing
    //      Ered = E + Rshell + LyPair(Z,A) - Econd  // origin of the FG without pairing
    alev = U_alev_ld(Z,A);
    //  Rho0 = 1.E0/Ered^1.25 * exp(2.E0 * sqr(alev * af_an * Ered))
    Rho1 = 1./pow(Ered,1.25) * exp(2. * sqrt(alev * af_an * (Ered - Eshift_shell))); 
    //  Rho1 = exp(2.E0 * sqr(alev * an_af * (Ered - Eshift_shell))) 
  } 
  else
  {
    if(Ipair == 0)
    {
      Rshell = - U_SHELL_exp(Z,A);
      if(Z > 92)  // hypothetical deviation in the macroscopic masses
        Rshell = Rshell - 0.1 * (Z - 92);
      Rmicro = Rshell;
    }
    else
    {
      Rmicro = - U_SHELL_EO_exp(Z,A);   // microscopic effects (shell and pairing)
      Rshell = - U_SHELL_exp(Z,A);   // only shell effect
      if(Z > 92)  // hypothetical deviation in the macroscopic masses
      {
        Rshell = Rshell - 0.1 * (Z - 92);
        Rmicro = Rmicro - 0.1 * (Z - 92);
      }
    }  
    //       Ered = E + LyPair(Z,A) - Econd
    Ered = E + Rshell - Rmicro - Econd;  // origin of the FG without pairing   
    alev = U_alev_ld(Z,A);
    Rho1 = 1./pow(Ered,1.25) * exp(2. * sqrt(alev * af_an * Ered));
    //    Rho1 = exp(2.E0 * sqr(alev * af_an * Ered))
  }
  // Print "*";Rmicro,Rshell,Eshift_shell"*"  // The values are positive
  return F_enhance * Rho1;
}

float GEF::U_Temp(float Z, float A, float E, int Ishell, int Ipair, float Tscale, float Econd)
{
  // Temperature (modified Gilbert-Cameron composite level density)    
  // KHS (10. 2. 2012)       
  float alev ; 
  float Eeff0,Eeff1,E1,Rho0,Rho1,TCT,TFG ;
  float fgamma = 0.055;      
  float RShell,RPair,Res;
  // Used global parameters: Tscale
  //   alev = U_alev_ld(Z,A) * 1.1   // Factor adjusted to high-energy prompt neutrons in U235(nth,f)
  alev = U_alev_ld(Z,A) * 0.95;  // " with the correction for non-constant T (FG range)
  //  alev = U_alev_ld(Z,A)

  if(Ishell == 1)
    RShell = U_SHELL(CInt(Z),CInt(A));
  else
    RShell = 0.0;    
  TCT = TEgidy(A,RShell,Tscale);  

  if(Ipair == 1)
    RPair = LyPair(CInt(Z),CInt(A));
  else
    RPair = 0.0; 
  Eeff0 = E - Econd + RPair + RShell*(1.0 - exp(-fgamma * E));

  if(Eeff0 > 0.5)
  {
    //       Eeff1 = Eeff0 + 0.1
    E1 = E + 0.1;
    Eeff1 = E1 - Econd + RPair + RShell*(1.0 - exp(-fgamma * E1));
    Rho0 = 1./pow(Eeff0,1.25) * exp(2. * sqrt(alev * Eeff0));
    Rho1 = 1./pow(Eeff1,1.25) * exp(2. * sqrt(alev * Eeff1));
    //         Rho0 = 1.E0/Eeff0 * exp(2.E0 * sqr(alev * Eeff0))
    //         Rho1 = 1.E0/Eeff1 * exp(2.E0 * sqr(alev * Eeff1))
    TFG = 0.1/ (log(Rho1) - log(Rho0));
  }
  else
    TFG = 0.0;

  Res = TCT;
  if(TFG > Res)
    Res = TFG;
  // If Res > 1.4 ThenRes = 1.4

  return Res;
}

float GEF::U_Temp2(float Z, float A, float E, float Rshell, float  Rpair, float Tscale, float Econd)
{
  // Temperature (modified Gilbert-Cameron composite level density)    
  // KHS (10. 2. 2012)       
  float alev;  
  float Eeff0,Eeff1,Rho0,Rho1,TCT,TFG; 
  static float fgamma = 0.055;      
  float Res;
  // Used global parameters: Tscale
  //  alev = U_alev_ld(Z,A) * 1.1;   // Factor adjusted to high-energy prompt neutrons in U235(nth,f)
  alev = U_alev_ld(Z,A) * 0.95;  //  with the correction for non-constant T (FG range)
  //  alev = U_alev_ld(Z,A);
  TCT = TEgidy(A,Rshell,Tscale);  
  Eeff0 = E - Econd + Rpair + Rshell*(1.0 - exp(-fgamma * E));
  //    Eeff0 = E - Econd + Lypair(CInt(Z),CInt(A)) + Rshell*(1.0 - exp(-fgamma * E))
  if(Eeff0 > 0.5)
  {
    Eeff1 = Eeff0 + 0.1;
    Rho0 = 1./pow(Eeff0,1.25) * exp(2. * sqrt(alev * Eeff0));
    Rho1 = 1./pow(Eeff1,1.25) * exp(2. * sqrt(alev * Eeff1));
    //         Rho0 = 1.E0/Eeff0 * exp(2.E0 * sqr(alev * Eeff0))
    //         Rho1 = 1.E0/Eeff1 * exp(2.E0 * sqr(alev * Eeff1))
    TFG = 0.1 / (log(Rho1) - log(Rho0));
  }
  else
    TFG = 0.0;

  Res = TCT;
  if(TFG > Res)
    Res = TFG;
  return Res;
}

float GEF::E0_GDR(float Z, float A)
{
  // Calculates the centroid energy of the GDR for spherical nucleus
  // according to the FRDM (ADNDT 59 (1995) 185 and PLB 670 (2008) 200)
  static float epsilon = 0.0768;
  static float J = 32.7;
  static float Q = 29.2;
  static float R0 = 1.16;
  static float mstar = 874;
  static float hbar = 197.3;
  float Aonethird,u,N,E0;

  // according to [9] in Phys. Lett. B 690 (2010) 473:
  return 18.0/pow(A,0.333333) + 25.0 /pow(A,0.1666667);

  // according to the FRDM (ADNDT 59 (1995) 185 and PLB 670 (2008) 200):
  //    Aonethird = A^0.333333
  //    N = A - Z
  //    u = (1-epsilon)/Aonethird * 3*J/Q
  //    E0_GDR = hbar /(R0*Aonethird)*sqr(8*J*A^2/ (mstar*4*N*Z) ) * _
  //      (1 + u - epsilon * (1+epsilon+3*u)/(1+epsilon+u))^(-1/2) 
}

float GEF::Width_GDR(float E0)
{
  // Spreading width of the GDR (Nucl. Phys. A 531 (1991) 27)
  return 1.99 * pow(E0/10,1.6);
}

float GEF::Efac_def_GDR(float Beta, float Gamma, float K)
{
  // Modification factors of the resonance energy due to triaxial deformation
  // Hill-Wheeler parameterisation (PRC 89 (1953) 1102)
  // Possible values for K:  K-2 = -1, 0, 1
  if(Beta == 0 && Gamma == 0)
    return 1.;
  else 
    return 1./(exp(sqrt(5./(4*pi))*Beta * cos(Gamma - 0.666667*(K-2)*pi)));
}     

float GEF::GgGtot(float Z, float A, float E, float Egamma)
{
  // From PRL 49 (1982) 434
  // Probability to emit a gamma of energy Egamma in competition to neutron emission
  float EG, GG, T, SN;
  EG = E0_GDR(Z,A);
  GG = Width_GDR(EG);
  T = U_Temp(Z,A,E,1,1,Tscale,Econd);
  SN = U_MASS(Z,A-1.) + LyPair(Z,A-1.) - (U_MASS(Z,A) + LyPair(Z,A)) ;
  return pow(Egamma,3) / pow(T,2) * exp((SN-Egamma)/T) * GG * EG / (pow(pow(Egamma,2) - pow(EG,2),2) + pow(GG,2) * pow(EG,2));
}

float GEF::E_next(float T1, float T2, float E1, float E2, float A1, float A2)
{
  // Samples the energy transfer in one step between two nuclei 
  // in thermal contact 
  // The energy transfer is only determined by the available phase space. 
  // Only one kind of nucleons considered! 
  // T1,T2 Temperatures of the two nuclei 
  // E1,E2 Initial energies of the two nuclei 
  // A1, A2 Mass numbers of the two nuclei 
  float E12;
  float Delta1,Delta2; /// Pairing gaps 
  float Delta_E1,Delta_E2;
  float E1final;       // Energy 1 after transfer 
  float E1mod,E2mod;

  // Assumed level densities: 
  // Even number of nucleons:
  // 1 ground state at energy E = - 2 Delta not considered
  // continuous level density above E = 0 : rho1,2 = a1,2 * exp(E1,2/T1,2) 

  E12 = E1 + E2; // Total energy 

  E1mod = E1;
  E2mod = E2;
  if(E1mod > E2mod)
  {
    Delta_E1 = Pexplim(-1./T1,0.0,E1mod);
    E1mod = E1mod - Delta_E1;
    E2mod = E2mod + Delta_E1;
    Delta_E2 = Pexplim(-1./T2,0.0,E2mod);
    E2mod = E2mod - Delta_E2;
    E1mod = E1mod + Delta_E2;
  }
  else
  {
    Delta_E2 = Pexplim(-1./T2,0.0,E2mod);
    E2mod = E2mod - Delta_E2;
    E1mod = E1mod + Delta_E2;
    Delta_E1 = Pexplim(-1./T1,0.0,E1mod);
    E1mod = E1mod - Delta_E1;
    E2mod = E2mod + Delta_E1;
  }
  E1final = E1mod;

  /* Select;
     When (E1 > E2) Do;
L3:
Delta_E1 = Pexplim(-1.E0/T1,0.0,E1);
E1final = E1 - Delta_E1;
Delta_E1 = Pexplim(-1.E0/T2,0.0,E12-E1final);
E1final = E1final + Delta_E1;
End;
When (E1 <= E2) Do;
L4:
Delta_E1 = Pexplim(-1.E0/T2,0.0,E12-E1);
E1final = E1 + Delta_E1;
Delta_E1 = Pexplim(-1.E0/T1,0.0,E1final);
E1final = E1final - Delta_E1;
End;
Otherwise Do;
List('This should not happen.');
End;
End; */

  return E1final;
}

float GEF::Pexplim(float R_lambda, float xmin, float xmax)
{
  // random number from an exponential between xmin and xmax 
  // decay constant: f(x) = exp(lambda * x) !!! 
  // xmin, xmax: limits for sampling 
  float umin, umax;  // xmin, xmax transformed 
  float u;  // help variable 
  float R_res;  // sampled value 

  if(abs(R_lambda) < 1.e-30)
    R_res = xmin + rn->Rndm() * (xmax-xmin);
  else
  {
    umin = exp(xmin*R_lambda);
    umax = exp(xmax*R_lambda);
    u = umin + rn->Rndm() * (umax-umin);
    R_res = 1./R_lambda * log(u);
  }
  return R_res;
}

float GEF::U_Even_Odd(int I_Channel, float PEO)
{
  //Creates even-odd fluctuations 
  float R;
  if(I_Channel%2 == 0)
    R = 1.0 + PEO;
  else
    R = 1.0 - PEO;
  return R;   
}

int GEF::EVEN_ODD(float R_ORIGIN,float R_EVEN_ODD)
{
  // Procedure to calculate I_OUT from R_IN in a way that         ///
  // on the average a flat distribution in R_IN results in a      ///
  // fluctuating distribution in I_OUT with an even-odd effect as ///
  // given by R_EVEN_ODD                                          ///
  // ------------------------------------------------------------ ///
  // EXAMPLES :                                                   ///
  // ------------------------------------------------------------ ///
  //    If R_EVEN_ODD = 0 :                                       ///
  //           CEIL(R_IN)  ----                                   ///
  //              R_IN ->                                         ///
  //            (somewhere in between CEIL(R_IN) and FLOOR(R_IN)) ///
  //           FLOOR(R_IN) ----       --> I_OUT                   ///
  // ------------------------------------------------------------ ///
  //    If R_EVEN_ODD > 0 :                                       ///
  //      The interval for the above treatment is                 ///
  //         larger for FLOOR(R_IN) = even and                    ///
  //         smaller for FLOOR(R_IN) = odd                        ///
  //    For R_EVEN_ODD < 0 : just opposite treatment              ///
  // ------------------------------------------------------------ ///
  // ------------------------------------------------------------ ///
  // On input:   R_ORIGIN    nuclear charge (real number)         ///
  //             R_EVEN_ODD  requested even-odd effect            ///
  // Intermediate quantity: R_IN = R_ORIGIN + 0.5                 ///
  // On output:  I_OUT       nuclear charge (integer)             ///
  // ------------------------------------------------------------ ///

  float R_IN,R_REST,R_HELP;
  float R_FLOOR;
  float R_MIDDLE;
  int I_OUT;

  R_EVEN_ODD = Min(R_EVEN_ODD,1.);
  R_IN = R_ORIGIN + 0.5;
  R_FLOOR = FLOOR(R_IN);
  if(abs(R_EVEN_ODD) < 1.e-3)
    I_OUT = R_FLOOR;	
  else
  {
    R_REST = R_IN - R_FLOOR;
    R_MIDDLE = R_FLOOR + 0.5;
    if(Mod(R_FLOOR, 2) == 0)  // even before modif. 
    {
      R_HELP = R_MIDDLE + (R_REST - 0.5) * 1. / Max(0.01,(1. + R_EVEN_ODD));
      R_HELP = Min(R_HELP,R_MIDDLE+1);
      R_HELP = Max(R_HELP,R_MIDDLE-1);
    }
    else  // odd before modification 
    {
      R_HELP = R_MIDDLE + (R_REST - 0.5) * 1.E0 / Max(0.01,(1. - R_EVEN_ODD));
      R_HELP = Min(R_HELP,R_MIDDLE+1);
      R_HELP = Max(R_HELP,R_MIDDLE-1);
    }
    I_OUT = FLOOR(R_HELP);
  }
  return I_OUT;  
}

float GEF::BFTF(float RZ, float RA, int I_Switch)
{
  // Fission barriers from Myers and Swiatecki, Thomas-Fermi model 
  //  I_Switch: 0: liquid-drop; 1: with shells and pairing, 
  //    2: averaged over pairing, 3: with shell and pairing + pairing gap at barrier 
  // 4: liquid-drop + g.s. shell, no Z correction
  float RN,RI,Rkappa,RS,RF,RX;
  float RX0 = 48.5428;
  float RX1 = 34.15;
  float RB ;
  int IZ,IA;
  float bftf;

  IZ = CInt(RZ);
  IA = CInt(RA);
  RN = RA - RZ;
  RI = (RN-RZ) / RA;
  Rkappa = 1.9 + (RZ - 80.) / 75.;
  RS = pow(RA,0.666667) * (1. - Rkappa * pow(RI,2));
  RX = pow(RZ,2) / (RA * (1. - Rkappa * pow(RI,2)));
  if(RX < 30)   // out of range 
    RF = 1.e10;
  if(RX > RX0)   // out of range 
    RF = 0.0;
  if(RX < RX1 && RX > 30) 
    RF = 0.595553 - 0.124136 * (RX - RX1);
  if(RX >= RX1 && RX <= RX0) 
    RF = 0.000199749 * pow(RX0 - RX,3);
  RB = RF * RS;

  switch(I_Switch)
  {
    case 0: //no shell, no pairing
      {
        bftf = RB;
      }break;
    case 1: // including even-odd staggering due to increased pairing strength at barrier
      {
        // Tentative modification from comparison with experimental fission barriers
        // (shell correction at the barrier?)
        if(RZ > 86.5)
          RB = RB - 0.15 * (RZ - 86.5);
        //    If RZ > 90 Then RB = RB + 0.3 * (RZ - 90.0);
        //    If RZ > 98 Then RB = RB - 0.15 * (RZ - 98.0); 
        if(RZ > 90)
          RB = RB + 0.35 * (RZ - 90.0);
        if(RZ > 93)
          RB = RB + 0.15 * (RZ - 93.0);
        if( RZ > 95)
          RB = RB - 0.25 * (RZ - 95.0);
        //    bftf = RB - U_Shell(IZ,IA)
        //    bftf = RB - U_Shell_exp(IZ,IA)
        bftf = RB - U_SHELL_EO_exp(IZ,IA) + LyPair(IZ,IA) * 14.0/12.0;
        //   bftf = RB - U_Shell_EO_exp(RZ,RA) - 14.E0 / sqr(Csng(RA))
        //       * Csng( ( (RZ+1) Mod 2 + (RA-RZ+1) Mod 2) )
      }break;
    case 2: // averaged over even-odd staggering
      {
        if(RZ > 86.5)
          RB = RB - 0.15 * (RZ - 86.5);
        if(RZ > 90)
          RB = RB + 0.35 * (RZ - 90.0);
        if(RZ > 93)
          RB = RB + 0.15 * (RZ - 93.0);
        if(RZ > 95)
          RB = RB - 0.25 * (RZ - 95.0); 
        bftf = RB - U_SHELL_exp(IZ,IA);
      }break;
    case 3: // like Case 1 but without increased pairing gap at barrier
      {
        if(RZ > 86.5)
          RB = RB - 0.15 * (RZ - 86.5);
        if(RZ > 90)
          RB = RB + 0.35 * (RZ - 90.0);
        if(RZ > 93)
          RB = RB + 0.15 * (RZ - 93.0);
        if(RZ > 95)
          RB = RB - 0.25 * (RZ - 95.0);
        bftf = RB - U_SHELL_EO_exp(IZ,IA);
      }break;
    case 4: // like case 3 but without Z correction
      {
        // This is the direct description from the topographic theorem.
        bftf = RB - U_SHELL_exp(IZ,IA);
      }break;
    default:
      {
        cout<< "Undefined option in BFTF"<<endl;
      }break;
  }
  /*  if(I_Switch == 0) 
      bftf = RB;
      else
      {
  // Tentative modification from comparison with experimental fission barriers
  // (shell correction at the barrier?)
  if(RZ > 86.5) RB = RB - 0.15 * (RZ - 86.5);
  //    if(RZ > 90) RB = RB + 0.3 * (RZ - 90.0);
  //    if(RZ > 98)  RB = RB - 0.15 * (RZ - 98.0); 
  if(RZ > 90) RB = RB + 0.35 * (RZ - 90.0);
  if(RZ > 93) RB = RB + 0.15 * (RZ - 93.0);
  if(RZ > 95) RB = RB - 0.25 * (RZ - 95.0);
  //    bftf = RB - U_SHELL(IZ,IA);
  //    bftf = RB - U_SHELL_exp(IZ,IA);
  bftf = RB - U_SHELL_EO_exp(IZ,IA) + LyPair(IZ,IA) * 14.0/12.0;
  }*/
  return bftf;
}
float GEF::BFTFA(float RZ, float RA, int I_Switch)
{
  // inner barrier height 
  float EA,BF0,Z4A,Z3A,DB;
  float coeff = 0.5;
  BF0 = BFTF(RZ,RA,I_Switch);
  // Z4A = RZ^4 / RA
  //  EB - EA from fit to Smirenkin barriers:
  //  V. M. Kupriyanov, K. K. Istekov, B. I. Fursov, G. N. Smirenkin
  //  Sov. J. Nucl. Phys. 32 (1980) 184
  //  DB = -10.3517 + 1.6027E-5 * Z4A + 5.4945E-11 * Z4A^2  // EA - EB  
  //  EB - EA from fit to data from Dahlinger et al. (KHS, 21. Dec. 2012)
  Z3A = pow(RZ,3) / RA;
  DB = -(5.40101 - 0.00666175*Z3A + 1.52531e-6*pow(Z3A,2));
  if(DB > 0.0)
    EA = BF0 - DB;
  else
    EA = BF0 ;
  return EA;
}

float GEF::BFTFB(float RZ, float RA, int I_Switch)
{
  // outer barrier height 
  float EB,BF0,Z4A,Z3A,DB;
  float coeff = 0.5;
  BF0 = BFTF(RZ,RA,I_Switch);
  // Z4A = RZ^4 / RA
  //  EB - EA from fit to Smirenkin barriers:
  //  V. M. Kupriyanov, K. K. Istekov, B. I. Fursov, G. N. Smirenkin
  //  Sov. J. Nucl. Phys. 32 (1980) 184
  //   DB = -10.3517 + 1.6027E-5 * Z4A + 5.4945E-11 * Z4A^2  // EA - EB
  //  EB - EA from fit to data from Dahlinger et al. (KHS, 21. Dec. 2012)
  Z3A = pow(RZ,3) / RA;
  DB = -(5.40101 - 0.00666175*Z3A + 1.52531e-6*pow(Z3A,2));  
  if(DB < 0.0)
    EB = BF0 + DB;
  else
    EB = BF0; 
  return EB;
}
float GEF::Gaussintegral(float R_x, float R_sigma)
{
  // Smoothed step function. Grows from 0 to 1 around R_x
  //    with a Gauss-integral function with given sigma
  float R_ret;
  // Note: The variable R_sigma = standard deviation / sqr(2) !
  R_ret = 0.5 + 0.5 * Erf(R_x / R_sigma);
  return R_ret;
}

float GEF::U_Box(float x, float sigma, float width)
{
  float y;
  // Note: The variable sigma = standard deviation / sqr(2) !
  y = Gaussintegral(x+0.5*width,sigma) - Gaussintegral(x-0.5*width,sigma);
  return y/width;
}

float GEF::U_Box2(float x, float sigma1, float sigma2, float width)
{
  float y;
  // Note: The variable sigma = standard deviation / sqr(2) !
  y = Gaussintegral(x+0.5*width,sigma2) - Gaussintegral(x-0.5*width,sigma1);
  return y/width;
}

float GEF::U_Gauss(float x, float sigma)
{
  float y;
  y = 1.0 / (sqrt(2.0 * pi) * sigma) * exp(-pow(x,2)/ ( 2.0 * pow(sigma,2) ) );
  return y;  
}

float GEF::U_Gauss_abs(float x, float sigma)
{
  float y;

  y = exp(-pow(x,2)/ ( 2.0 * pow(sigma,2) ) );
  return y;
}

float GEF::U_Gauss_mod(float x,float sigma)
{
  // Gaussian with Sheppard correction
  float y;
  float sigma_mod;      
  sigma_mod = sqrt(pow(sigma,2) + 1./12.);

  y = 1.0 / (sqrt(2.0 * pi) * sigma_mod) * exp(-pow(x,2)/ ( 2.0 * pow(sigma_mod,2) ) );
  return y;  
}

float GEF::PBox(float Mean, float Sigma, float Bottom)
{
  // Rectangular distribution folded with a Gaussian distribution   
  float R;
  R = PGauss(Mean,Sigma);
  R = R + (rn->Rndm()-0.5)*Bottom;
  return R;
}

float GEF::PBox2(float Mean, float Sigma1, float Sigma2,  float Bottom)
{
  // Rectangular distribution folded with a Gaussian distribution. 
  // One wing is steeper.  
  // Sigma1 = lower side, Sigma2 = upper side
  float Sigma,R;
  Sigma = Max(Sigma1,Sigma2);

  R = PGauss(Mean,Sigma);
  R = R + (rn->Rndm()-0.5)*Bottom;
  if(Sigma1 < Sigma2)
    if(R < Mean - 0.5*Bottom)
      if( rn->Rndm() > exp( -pow(R - Mean + 0.5*Bottom,2) / (2.0 * pow(Sigma1,2)) )/ exp( -pow(R - Mean + 0.5*Bottom,2) / (2.0 * pow(Sigma2,2)) )) 
        R = Mean - 0.5*Bottom + (Mean - 0.5*Bottom - R);

  if(Sigma2 <= Sigma1)
    if( R > Mean + 0.5*Bottom)
      if( rn->Rndm() > exp( -pow(R - Mean - 0.5*Bottom,2) / (2.0 * pow(Sigma2,2)) ) / exp( -pow(R - Mean - 0.5*Bottom,2) / (2.0 * pow(Sigma1,2)) ) )
        R = Mean + 0.5*Bottom - (R - Mean - 0.5*Bottom);

  return R;
}

float GEF::PPower(int Order, float Rmin, float Rmax)
{
  // Random generator of a power function: (y = x^Order -> x_random = RND^(1/(Order+1))
  // PPower = 0 at Rmin to PPower = Ymax at Rmax 
  float R;
  R = Rmin + (Rmax-Rmin) * pow(rn->Rndm(),1.0/(Order+1));  
  return R; 
}

float GEF::PPower_Griffin_v(int Order, float Rmin, float Rmax)
{
  // Random generator of a power function: (y = x^Order -> x_random = RND^(1/(Order+1))
  // PPower = 0 at Rmin to PPower = Ymax at Rmax 
  float R,v_particle,RRND;
  bool Repeat_Griffin=true;
  while(Repeat_Griffin)
  {
    Repeat_Griffin=false;
    R = Rmin + (Rmax-Rmin) * pow(rn->Rndm(),1.0/(Order));  
    v_particle = sqrt(abs((R-Rmax)/(Rmin-Rmax)));
    RRND = rn->Rndm();
    if(RRND > v_particle)
      Repeat_Griffin=true;
  }  
  return R;
}

float GEF::PPower_Griffin_E(int Order, float Rmin,float  Rmax)
{
  // Random generator of a power function: (y = x^Order -> x_random = RND^(1/(Order+1))
  // PPower = 0 at Rmin to PPower = Ymax at Rmax 
  float R,E_particle;
  bool Repeat_Griffin = true;
  while(Repeat_Griffin)
  {
    Repeat_Griffin = false;
    R = Rmin + (Rmax-Rmin) * pow(rn->Rndm(),1.0/(Order));  
    E_particle = (R-Rmax)/(Rmin-Rmax);
    if(rn->Rndm() > E_particle)
      Repeat_Griffin = true;
  }
  return R;
}

float GEF::PGauss(float Mean, float Sigma)
{
  // Box-Mueller method
  static int ISet = 0;
  float V1,V2,R,Fac,GasDev,Result;
  static float GSet;
  if( ISet == 0)
  {
    bool Repeat=true;
    while(Repeat)
    {
      Repeat = false;
      V1 = 2. * rn->Rndm() - 1.;
      V2 = 2. * rn->Rndm() - 1.;
      R = pow(V1,2) + pow(V2,2);
      if((R >= 1.) | (R == 0.0))
        Repeat = true;
    }
    Fac = sqrt(-2. * log(R)/R);
    GSet = V1 * Fac;
    GasDev = V2 * Fac;
    ISet = 1;
  }
  else
  {
    GasDev = GSet;
    ISet = 0;
  }
  Result = Sigma * GasDev;
  return Mean + Result;
}

float GEF::PLinGauss(float R_Sigma)
{
  // Random-number generator for linear * Gaussian function 
  // Distribution of nuclear angular momenta 
  float R_Res,B_rms;
  B_rms = R_Sigma / sqrt(2.0); // Because the sum of two PGauss functions increases the width.
  R_Res = abs(PGauss(0,B_rms)) + abs(PGauss(0,B_rms));
  R_Res = R_Res + B_rms/4. * (1. - exp(-R_Res/R_Sigma));
  // correction of shape (approximative)
  return R_Res;
}

float GEF::U_LinGauss(float x, float R_Sigma)
{
  // Gaussian times a linear function 
  // Not normalized! 
  float R_Res;
  if(R_Sigma > 0.0)
    R_Res = x * exp(-pow(x,2)/(2.0 * pow(R_Sigma,2)));
  else
    R_Res = 0.0;
  return R_Res;
}

float GEF::PExp(float R_Tau)
{
  // Random-number generator for an exponential distribution 
  float X1,R_Res;
  bool Again=true;
  while(Again)
  {
    Again = false;
    X1 = rn->Rndm();
    if(X1 > 1.e-10 && X1 < 0.99999)  // for avoiding numerical problems
      R_Res = - R_Tau * log(X1);
    else
      Again = true;
  }

  return R_Res;
}

float GEF::PMaxwell(float R_T)
{
  // Random-number generator for a surface Maxwell distribution 
  // y = x * exp(-x/T) 
  double R_Res,R_T_int;
  R_T_int = R_T;
  R_Res = -R_T_int * (log(rn->Rndm()) + log(rn->Rndm()));
  return R_Res;
}

float GEF::PMaxwellv(float R_T)
{
  // Random generator according to a distribution similar to a 
  // Maxwell distribution with quantum-mech. x-section for neutrons  
  // (approximation by KHS) 
  // Y = SQRT(X) * EXP(-X/T) 
  float EN;
  EN = 2. * R_T * sqrt(log(rn->Rndm()) * log(rn->Rndm()));
  return EN;     
}

float GEF::PMaxwellMod(float R_T, float R_A)
{
  // Random generator according to a distribution similar to a 
  // Maxwell distribution with quantum-mech. x-section for neutrons 
  // (approximation by KHS) 
  // Y = SQRT(X) * EXP(-X/T) 
  float EN;
  if(rn->Rndm() < 3.3 / sqrt(R_A))   //according to PR-116-683 (Dostrowsky et al.)
    EN = 2. * R_T * sqrt(log(rn->Rndm()) * log(rn->Rndm()));
  else
    EN = -R_T * (log(rn->Rndm()) + log(rn->Rndm()));
  return EN;
}

float GEF::Round(float R, int N)
{
  // R   Input value
  // N   Number of significant digits
  float RN10, Rred, Rextended, Rrounded, Rout, Rabs;
  int N10,Isign;
  if(R == 0)
    return 0;
  else
  {
    Isign = Sgn(R);
    Rabs = abs(R);
    N10 = CInt(Log10(Rabs));
    RN10 = pow(10,N10);
    Rred = Rabs / RN10;
    Rextended = Rred * pow(10,N-1);
    Rrounded = CInt(Rextended + 0.5);
    Rout = Rrounded / pow(10,N-1) * RN10;
    return Isign * Rout;
  } 
}

long int GEF::Modulo(unsigned long int I, unsigned long int J)
{
  unsigned long int Iratio,Iresult;
  Iratio = I / J;
  Iresult = I - J * Iratio;
  return Iresult;
}
/*
   int GEF::PLoss(unsigned long int IL)
   {
// Extracts and returns number of prompt protons 
Dim As string Ctest;
int NP;
Ctest = Oct(IL);
NP = 0;
for(int i=1;i<=Len(Ctest);i++)
if((Mid(Ctest,i,1) == "2") | (Mid(Ctest,i,1) == "4"))
NP = NP + 1;
return NP;
}
 */
bool GEF::IsValid(int Z, int A)
{
  bool Ivalid;
  Ivalid = 1;
  //   If I_A / I_Z < 210.E0/90.E0 
  if((Z < 76) || (Z > 120))
    Ivalid = 0;  
  return Ivalid;  
  //return 1;
}


bool GEF::U_Valid(int I_Z, int I_A)
{
  bool Ivalid;
  Ivalid = 1;
  //   If I_A / I_Z < 210.E0/90.E0 
  if((I_A / I_Z < 160. / 76.) | (I_A / I_Z > 250./90.)) 
    Ivalid = 0;
  if((I_Z < 76) | (I_Z > 120))
    Ivalid = 0;  
  return Ivalid;  
  //return 1;
}

float GEF::U_Delta_S0(int I_Z, int I_A)
{
  // I_Z and I_A refer to the fissioning nucleus
  float Delta;
  Delta = 0;

  if(I_Z == 89)
    if(I_A == 226)
      Delta = -0.3;

  if(I_Z == 90)
  {
    if(I_A == 228)
      Delta = 0.2;
    if(I_A == 229)
      Delta = 0.4;
    if(I_A == 230)
      Delta = 0.7;
    if(I_A == 231)
      Delta = 0.8;
    if(I_A == 232)
      Delta = 0.9;
    if(I_A == 233)
      Delta = 0.9; 
  }

  if(I_Z == 92)
    Delta = 0.2;    //x
  if(I_Z == 92 && I_A == 233)
    Delta = 0.4;    //x
  if(I_Z == 92 && I_A == 234)
    Delta = 0.4;    //x 

  if(I_Z >= 93)
    Delta = -0.3;  //x

  return Delta;    
}


void GEF::PotCurv_FM(void)
{
  /* Potential curvatures of fission modes */
  // For the width of the mass distribution (potential between saddle and scission):
  // Print Spin_pre_fission,  P_I_rms_CN 
  R_Z_Curv_S0 = 8. / pow(Csng(I_Z_CN),2) * Masscurv(Csng(I_Z_CN), Csng(I_A_CN), Spin_pre_fission, kappa);
  // For the yields of the fission channels (potential near saddle):
  R_Z_Curv1_S0 = 8. / pow(Csng(I_Z_CN),2) * Masscurv1(Csng(I_Z_CN), Csng(I_A_CN), 0.0, kappa);
  R_A_Curv1_S0 = 8. / pow(Csng(I_A_CN),2) * Masscurv1(Csng(I_Z_CN), Csng(I_A_CN), 0.0, kappa);
}

void GEF::EnergyTrans(void)
{
  /*
     switch(Emode)
     {
     case 0:   // Energy above outer barrier given
     {
     R_E_exc_Eb = R_E_exc_used;
     R_E_exc_GS = R_E_exc_used + BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
     }break;
     case 1:  // Energy above ground state given
     {
     R_E_exc_Eb = R_E_exc_used - BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
     R_E_exc_GS = R_E_exc_used;
     }break;
     case 3:   // Energy above ground state given
     {
     R_E_exc_Eb = R_E_exc_used - BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
     R_E_exc_GS = R_E_exc_used;
     }break;
     case -1:   // Energy above ground state given
     {
     R_E_exc_Eb = R_E_exc_used - BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
     R_E_exc_GS = R_E_exc_used;
     }break;
     case 2:     // kinetic energy of neutron given (SN = neutron separation energy)
     {
  //    SN = (U_Mass(Csng(I_Z_CN),Csng(I_A_CN-1)) + Lypair(I_Z_CN,I_A_CN-1)) _
  //       -(U_Mass(Csng(I_Z_CN),Csng(I_A_CN)) + Lypair(I_Z_CN,I_A_CN))
  //    R_E_exc_GS = R_E_exc_used + SN 
  SN = AME2012(I_Z_CN,I_A_CN-1) - AME2012(I_Z_CN,I_A_CN);
  R_E_exc_GS = R_E_exc_used * ((P_A_CN-1) / P_A_CN) + SN;                                            //           target CN           
  R_E_exc_Eb = R_E_exc_GS - BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
  }break;
  case 12:     // kinetic energy of proton given (Sprot = proton separation energy)
  {
  //    Sprot = (U_Mass(Csng(I_Z_CN-1),Csng(I_A_CN-1)) + Lypair(I_Z_CN-1,I_A_CN-1)) _
  //       -(U_Mass(Csng(I_Z_CN),Csng(I_A_CN)) + Lypair(I_Z_CN,I_A_CN))
  //    R_E_exc_GS = R_E_exc_used + Sprot 
  Sprot = AME2012(I_Z_CN-1,I_A_CN-1) - AME2012(I_Z_CN,I_A_CN);
  R_E_exc_GS = R_E_exc_used * ((P_A_CN-1) / P_A_CN) + Sprot   ; 
  R_E_exc_Eb = R_E_exc_GS - BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
  }break;
  case 13:   // list of energies from file
  {
  R_E_exc_Eb = R_E_exc_used - BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
  R_E_exc_GS = R_E_exc_used ;
  }break;
  }
   */
  R_E_exc_Eb = R_E_exc_used - BFTFB(Csng(I_Z_CN),Csng(I_A_CN),1);
  R_E_exc_GS = R_E_exc_used;

}

void GEF::BarriersEx_FM(void)
{

  E_exc_S0_prov = R_E_exc_Eb;

  /* Additional influence of N=82 assumed */
  Delta_NZ_Pol = 82./50. - Csng(I_N_CN)/Csng(I_Z_CN);
  R_Shell_S1_eff = P_Shell_S1 * (1.0 - P_att_rel * P_att_pol * abs(Delta_NZ_Pol));

  //   R_Shell_S1_eff = P_Shell_S1 * _
  //           max(1.0 - P_att_rel,(1.0 - P_att_rel* _
  //                   ( Abs(Delta_NZ_Pol)/P_att_Pol  + (Delta_NZ_Pol/P_att_Pol2)^2 _
  //                   + Abs(Delta_NZ_Pol/P_att_Pol3)^3)))
  // Print "4335 "; max(1.0 - P_att_rel,(1.0 - P_att_rel* _
  //                   ( Abs(Delta_NZ_Pol)/P_att_Pol  + (Delta_NZ_Pol/P_att_Pol2)^2 _
  //                   + Abs(Delta_NZ_Pol/P_att_Pol3)^3)))
  // Print  Abs(Delta_NZ_Pol)/P_att_Pol, _
  //         (Delta_NZ_Pol/P_att_Pol2)^2, _
  //        Abs(Delta_NZ_Pol/P_att_Pol3)^3
  //sleep                     

  /* In Pu, the Z=50 shell meets Z=44 in the light fragment. */
  /* A deformed shell at Z=44 is assumed to explain the enhancement _ 
     of the S1 channel around Pu */
  /* This very same shell automatically produces the double-humped */
  /* mass distribution in 180Hg */   
  //    S1_enhance = P_Shell_SL4 + _
  //            (Csng(I_Z_CN) - ZC_Mode_1 - ZC_Mode_4L)^2 * P_Z_Curv_SL4 
  // 50 instead of ZC_Mode_1, to eliminate the influenc of the mass 
  //(in agreement with experiment, e.g. 238U(nfast,f) ):
  S1_enhance = P_Shell_SL4 + pow(Csng(I_Z_CN) - 50.0 - ZC_Mode_4L,2) * P_Z_Curv_SL4; 

  //   S1_enhance = P_Shell_SL4 * _
  //       U_Gauss_abs(Csng(I_Z_CN) - 50.0 - 0.3 - ZC_Mode_4L,P_Z_Sigma_SL4)
  //Print "4396: U_Gauss_abs";S1_enhance/P_Shell_SL4,R_SHell_S1_eff
  //sleep        

  //   S1_enhance = P_Shell_SL4 + _
  //           (Csng(I_Z_CN) - ZC_Mode_1 - ZC_Mode_4L)^2 * P_Z_Curv_SL4 
  if(S1_enhance > 0)
    S1_enhance = 0;

  if(P_Z_CN == 91)
    S1_enhance = S1_enhance + 0.3;
  if(P_Z_CN == 90)
    S1_enhance = S1_enhance + 0.6;

  //Print "4384 "; P_Shell_SL4, U_Gauss(Csng(I_Z_CN) - 50.0 - ZC_Mode_4L,0.3), S1_enhance   
  //sleep

  // Print "3933"
  // Print "ZC_Mode_1,ZC_Mode_4",ZC_Mode_1,ZC_Mode_4
  // Print "Delta-Z S1-S4, S1_enhance",I_Z_CN-ZC_Mode_1 - ZC_Mode_4L, S1_enhance      
  // Print "3951"
  R_Shell_S1_eff = R_Shell_S1_eff + S1_enhance;
  ;
  // Print I_Z_CN-ZC_Mode_1-ZC_Mode_4L, S1_enhance,R_Shell_S1_eff

  /* The high TKE of S1 in 242Pu(sf) (and neighbours) is obtained by assuming */
  /* that the Z=44 shell reduces the deformation of the light fragment. */
  for(int I = 10;I<=(I_Z_CN - 10);I++)
  {
    Z1 = Csng(I);
    A1 = (Z1 - 0.5) / Csng(I_Z_CN) * Csng(I_A_CN); /* polarization roughly considered */
    //    Beta(1,1,Z1) = Beta(1,1,Z1) + 0.15 * S1_enhance   /* "light" fragment */
    Beta[2][1][I] = exp(S1_enhance) * Beta[2][1][I] + (1.-exp(S1_enhance)) * (Beta[2][1][I]-0.25);
    Beta[2][1][I] = Max(Beta[2][1][I],0.0);
    E_defo = LyMass(Z1,A1,Beta[2][1][I]) - LyMass(Z1,A1,0.0);
    Edefo[2][1][I] = E_defo ;/* "light" fragment */
  } 

  // Influence of S2 shell in complementary fragment
  // May be called "S12 fission channel"
  T_asym_Mode_2 = 0.5;
  SigZ_Mode_2 = sqrt(0.5 * T_asym_Mode_2/(P_Z_Curv_S2));
  SigA_Mode_2 = SigZ_Mode_2 * Csng(I_A_CN) / Csng(I_Z_CN);
  S1_enhance_S2 = P_Shell_S2 * U_Box(Csng(P_A_CN) - AC_Mode_2 - AC_Mode_1, SigA_Mode_2,P_A_Width_S2) *P_A_Width_S2;

  if(S1_enhance_S2 < 0.01)
    //  Print "S1_enhance_S2 ";S1_enhance_S2  
    R_Shell_S1_eff = R_Shell_S1_eff + S1_enhance_S2;   

  R_Shell_S2_eff = P_Shell_S2;     

  // Overlap of S3 and shell in light fragment  
  R_Shell_S3_eff = P_Shell_S3 * (1.0 - 0.8* PZ_S3_olap_curv * pow(Csng(I_Z_CN) - ZC_Mode_3 - PZ_S3_olap_pos,2));
  //  Print "4450 "; P_Shell_S3 * (1.0 - 0.8 * PZ_S3_olap_curv _
  //         * (Csng(I_Z_CN) - ZC_Mode_3 - PZ_S3_olap_pos)^2) , ZC_Mode_3 
  //  sleep            
  //  R_Shell_S3_eff = -5.605
  //        * (Csng(I_Z_CN) - 60.5E0 - PZ_S3_olap_pos)^2)
  R_Shell_S3_eff = Min(R_Shell_S3_eff,0.0);    

  // Additional empirical dependence on N/Z  
  // R_Shell_S3_eff = R_Shell_S3_eff - _
  //       1 * ( (I_A_CN-I_Z_CN)/I_Z_CN - (236-92)/92)  
  //       5 * ( (I_A_CN-I_Z_CN)/I_Z_CN - (236-92)/92)  

  //   R_Shell_S4_eff = 2.0 * (P_Shell_S4 + P_Z_Curv_S4*(ZC_Mode_4 - ZC_Mode_0)^2)
  R_Shell_S4_eff = 2.0 * (P_Shell_S4 + P_Z_Curv_S4 * pow(ZC_Mode_4 - ZC_Mode_0,2));     
  // overlap of S4 in both fragments       
  if(R_Shell_S4_eff > P_Shell_S4)
    R_Shell_S4_eff = P_Shell_S4; 
  // no overlap at large distance

  E_LD_S1 = R_A_Curv1_S0 * pow(Csng(I_A_CN)/Csng(I_Z_CN)*(ZC_Mode_1 - ZC_Mode_0) ,2);
  B_S1 = E_LD_S1 + R_Shell_S1_eff;
  E_exc_S1_prov = E_exc_S0_prov - B_S1;

  E_LD_S2 = R_A_Curv1_S0 * pow(Csng(I_A_CN)/Csng(I_Z_CN)*(ZC_Mode_2 - ZC_Mode_0) ,2);
  B_S2 = E_LD_S2 + R_Shell_S2_eff;
  E_exc_S2_prov = E_exc_S0_prov - B_S2;   

  E_LD_S3 = R_A_Curv1_S0 * pow(Csng(I_A_CN)/Csng(I_Z_CN)*(ZC_Mode_3 - ZC_Mode_0) ,2);
  B_S3 = E_LD_S3 + R_Shell_S3_eff;
  E_exc_S3_prov = E_exc_S0_prov - B_S3; 


  if(I_A_CN < 220)  // Only here S4 is close enough to symmetry to have a chance
  {
    E_LD_S4 = R_A_Curv1_S0 * pow(Csng(I_A_CN)/Csng(I_Z_CN)*(ZC_Mode_4 - ZC_Mode_0) ,2);
    //     R_Shell_S4_eff = R_Shell_S4_eff * (1.0 + P_S4_NZmod * (Csng(I_A_CN)/Csng(I_Z_CN) - (180/80)) ) // variation with A/Z
    R_Shell_S4_eff = R_Shell_S4_eff * (1.0 + P_S4_NZmod * (Csng(I_A_CN-I_Z_CN) - (112)) ); // variation with N
    //R_Shell_S4_eff = R_Shell_S4_eff * (1.0 - 0.09 * (Csng(I_Z_CN) - 80) ) // variation with Z      
    B_S4 = E_LD_S4 + R_Shell_S4_eff;
    E_exc_S4_prov = E_exc_S0_prov - B_S4;
  }
  else
  {
    B_S4 = 9999;
    E_exc_S4_prov = - 9999;  
  }

  /* Mode 11 (overlap of channel 1 in light and heavy fragment */
  /* Potential depth with respect to liquid-drop potential: B_S11 */
  //    B_S11 = 2.E0 * (R_Shell_S1_eff + De_Defo_S1 _
  //             + P_Z_Curv_S1 * (ZC_Mode_1 - ZC_Mode_0)^2 ) - De_Defo_S1 
  B_S11 = 2. * (R_Shell_S1_eff + P_Z_Curv_S1 * pow(ZC_Mode_1 - ZC_Mode_0,2) );  
  // Sum of S1 shells in both fragments exact at symmetry    

  // Print "4475 ";R_Shell_S1_eff, B_S11
  // Print ZC_Mode_0, ZC_Mode_1, P_Z_Curv_S1 * (ZC_Mode_1 - ZC_Mode_0)^2

  // If B_S11 (see above) is higher than the shell at symmetry from only one fragment
  //  If B_S11 > R_Shell_S1_eff + P_Z_Curv_S1 * (ZC_Mode_1 - ZC_Mode_0)^2 Then
  //    B_S11 = Min(B_S11,R_Shell_S1_eff + P_Z_Curv_S1 * (ZC_Mode_1 - ZC_Mode_0)^2 )
  //  End If   

  DES11ZPM = 0;             
  // The S1 shells in the two fragments must be rather close to form  one pocket:
  if(B_S11 < (R_Shell_S1_eff + Level_S11))  
    // Lowering of the zero-point motion grows with the width of the potential pocket:
    //   DES11ZPM = -0.6 * Abs(ZC_Mode_1 - ZC_Mode_0)
    DES11ZPM = -0.8 * abs(ZC_Mode_1 - ZC_Mode_0);  // Fits the mass distr. of 258Fm(sf)

  /* Lowering of effective barrier by lower ZPM due to larger width in
     partial overlap region (shells in light and heavy fragment) */
  //   DES11ZPM = Level_S11 * Min(Abs(ZC_Mode_1 - ZC_Mode_0),4.E0*P_Z_Curv_S1)
  //   DES11ZPM = -0.2 * Abs(ZC_Mode_1 - ZC_Mode_0)  

  // Print "4473: "; R_Shell_S1_eff, B_S11, DES11ZPM
  // Sleep   

  B_S11 = B_S11 + DES11ZPM;

  //  If B_S11 > R_Shell_S1_eff + 0.5E0 Then 
  //   If B_S11 > R_Shell_S1_eff + Level_S11 Then
  //     B_S11 = 100   // S1 and S11 are exclusive
  //   Else
  //     B_S11 = Min(B_S11,R_Shell_S1_eff)  
  //   End If  

  E_exc_S11_prov = E_exc_S0_prov - B_S11;

  /* Mode 22 (overlap of channel 2 in light and heavy fragment */
  /* Potential depth with respect to liquid-drop potential: B_S22 */

  //   B_S22 = 2.E0 * (E_LD_S2 + P_Shell_S2) _
  //       + 2.E0 * P_Z_Curv_S2 * (ZC_Mode_2 - ZC_Mode_0)^2   /* Parabola */
  //Print E_LD_S2,P_Shell_S2,P_Z_Curv_S2,ZC_Mode_2,ZC_Mode_0   
  B_S22 = 2. * R_Shell_S2_eff  * U_Box(Csng(P_A_CN)/2.0 - AC_Mode_2,SigA_Mode_2,P_A_Width_S2) * P_A_Width_S2;
  // The integral of U_Box is normalized, not the height! 
  //    If Abs((P_A_CN/2.E0) - AC_Mode_2) > P_A_Width_S2 Then B_S22 = 9999   
  if(P_A_CN < 226)
    B_S22 = 9999; 

  E_exc_S22_prov = E_exc_S0_prov - B_S22;

  E_Min_Barr = Min(0.0,B_S1);
  E_Min_Barr = Min(E_Min_Barr,B_S2);
  E_Min_Barr = Min(E_Min_Barr,B_S3);
  E_Min_Barr = Min(E_Min_Barr,B_S4);
  E_Min_Barr = Min(E_Min_Barr,B_S11);
  E_Min_Barr = Min(E_Min_Barr,B_S22);

  /* Energy minus the height of the respective fission saddle */
  E_Exc_S0 = E_exc_S0_prov + E_Min_Barr - Delta_S0;
  E_Exc_S1 = E_exc_S1_prov + E_Min_Barr;
  E_Exc_S2 = E_exc_S2_prov + E_Min_Barr;
  E_Exc_S3 = E_exc_S3_prov + E_Min_Barr;
  E_Exc_S4 = E_exc_S4_prov + E_Min_Barr;
  E_Exc_S11 = E_exc_S11_prov + E_Min_Barr;
  E_Exc_S22 = E_exc_S22_prov + E_Min_Barr;

  /* Energy above the lowest fission saddle */
  E_exc_Barr = Max(E_Exc_S0,E_Exc_S1);
  E_exc_Barr = Max(E_exc_Barr,E_Exc_S2);
  E_exc_Barr = Max(E_exc_Barr,E_Exc_S3);
  E_exc_Barr = Max(E_exc_Barr,E_Exc_S4);
  E_exc_Barr = Max(E_exc_Barr,E_Exc_S11);
  E_exc_Barr = Max(E_exc_Barr,E_Exc_S22);
}

void GEF::TCollective(void)
{
  /* Collective temperature used for calculating the widths
     in mass asymmetry and charge polarization */

  if(E_Exc_S0 < 0)
    E_tunn = -E_Exc_S0;
  else
    E_tunn = 0;
  R_E_exc_eff = Max(0.1,E_Exc_S0);
  //  T_Coll_Mode_0 = TFCOLL * R_E_exc_eff + _  // empirical, replaced by TRusanov 
  T_Coll_Mode_0 = TCOLLFRAC * (De_Saddle_Scission(pow(Csng(I_Z_CN),2) / pow(Csng(I_A_CN),0.33333),ESHIFTSASCI_coll) - E_tunn);
  T_Coll_Mode_0 = Max(T_Coll_Mode_0,0.0);

  //Print "4596: De_SS, E_tunn, T_Coll ";De_Saddle_Scission(I_Z_CN^2/I_A_CN^0.3333,ESHIFTSASCI_coll),E_tunn,T_Coll_Mode_0    

  /* Temperature description fitting to the empirical systematics of Rusanov et al. */
  /* Here from Ye. N. Gruzintsev et al., Z. Phys. A 323 (1986) 307 */    
  /* Empirical description of the nuclear temperature according to the */
  /* Fermi-gas description. Should be valid at higher excitation energies */
  float T_Rusanov;
  T_Rusanov = TRusanov(R_E_exc_eff,Csng(I_A_CN)); 
  //Print "Temperatures, (GEF, Total, Rusanov): ", T_Coll_Mode_0, TFCOLL * R_E_exc_eff, T_Rusanov
  //Print "R_E_exc_eff ",R_E_exc_eff

  T_Coll_Mode_0 = Max(T_Coll_Mode_0,T_Rusanov);
  /* Transition vom const. temp. to Fermi gas occurs around 20 MeV by MAX function */
  //    T_Pol_Mode_0 = T_Pol_Red * T_Coll_Mode_0

  // Application of the statistical model, intrinsic temperature at saddle
  T_Pol_Mode_0 = U_Temp(0.5 * Csng(I_Z_CN),0.5 *Csng(I_A_CN), R_E_exc_eff, 0, 0, Tscale, Econd);
  //    T_asym_Mode_0 = Sqr(T_Coll_Mode_0^2 + (6E0*TCOLLMIN)^2)  
  T_asym_Mode_0 = sqrt(pow(T_Coll_Mode_0,2) + pow(1.0*TCOLLMIN,2));
  //Print "4124: T_Coll_Mode_0"; T_Coll_Mode_0
  //sleep  

  E_POT_SCISSION = (De_Saddle_Scission(pow(Csng(I_Z_CN),2) / pow(Csng(I_A_CN),0.33333),ESHIFTSASCI_intr) - E_tunn) + Epot_shift; 
  E_diss_Scission = EDISSFRAC * E_POT_SCISSION;        
  //Print "4054:";EDISSFRAC,E_POT_SCISSION,E_diss_Scission                      

  /* Suppression of S1 fission channel at very low excitation energy at scission */
  /* The idea behind is that the binding energy at scission is such that the
     scission configuration cannot be reached with the available excitation energy. */
  //   EeffS1 = Max(E_Exc_S1,0.0) + EDISSFRAC * E_POT_SCISSION
  //   EeffS1 = Max(0.0,EeffS1)

  // Print "4104", U_Mass(I_Z_CN,I_A_CN); _
  //       2 * U_Mass(I_Z_CN/2.0,I_A_CN/2.0) + 1.44*(I_Z_CN/2.0)^2 / _
  //               (1.5 * ( (I_A_CN/2)^0.333333 + (I_A_CN/2)^0.333333) + DNECK ); EeffS1, _
  //        - 2 * U_Mass(I_Z_CN/2.0,I_A_CN/2.0) - 1.44*(I_Z_CN/2.0)^2 / _
  //               (1.5 * ( (I_A_CN/2)^0.333333 + (I_A_CN/2)^0.333333) + DNECK ) + _
  //         + Max(E_Exc_S1,0.0) + EDISSFRAC * E_POT_SCISSION _     
  //         + U_Mass(I_Z_CN,I_A_CN)           
  //   EeffS1 = - 2 * U_Mass(I_Z_CN/2.0,I_A_CN/2.0) - 1.44*(I_Z_CN/2.0)^2 / _
  //              (1.6 * ( (I_A_CN/2)^0.333333 + (I_A_CN/2)^0.333333) + DNECK ) + _
  //          Max(E_Exc_S1,0.0) + EDISSFRAC * E_POT_SCISSION _     
  //        + U_Mass(I_Z_CN,I_A_CN)     
  //   If EeffS1 < ETHRESHSUPPS1 Then
  //                + 2.E0 * ESIGSUPPS1 Then 
  //     E_Exc_S1 = E_Exc_S1 + EeffS1 - ETHRESHSUPPS1
  //        0.5E0 * 1.5 * 12.E0 / Sqr(132.E0) * Gaussintegral(ETHRESHSUPPS1 - EeffS1,ESIGSUPPS1)
  ////         0.5E0 * 4.E0 * 12.E0 / Sqr(132.E0) * Gaussintegral(ETHRESHSUPPS1 - EeffS1,ESIGSUPPS1)
  //   End If
  //   If EeffS2 < ETHRESHSUPPS1 + 2.E0 * ESIGSUPPS1 Then
  //     E_Exc_S1 = E_Exc_S1 - _
  //        0.5E0 * 1.5 * 12.E0 / Sqr(132.E0) * Gaussintegral(ETHRESHSUPPS1 - EeffS2,ESIGSUPPS1)
  //        0.5E0 * 4.E0 * 12.E0 / Sqr(132.E0) * Gaussintegral(ETHRESHSUPPS1 - EeffS2,ESIGSUPPS1)
  //   EndIf

  T_low_S1_used = T_low_S1;

  T_Coll_Mode_1 = TFCOLL * Max(E_Exc_S1,0.) + TCOLLFRAC * (De_Saddle_Scission(pow(I_Z_CN,2) / pow(Csng(I_A_CN),0.33333),ESHIFTSASCI_coll) - E_tunn);
  T_Coll_Mode_1 = Max(T_Coll_Mode_1,0.0);
  //    T_Pol_Mode_1 = T_Pol_Red * T_Coll_Mode_1
  T_Pol_Mode_1 = T_Pol_Mode_0;
  T_asym_Mode_1 = sqrt(pow(T_Coll_Mode_1,2) + pow(4.0*TCOLLMIN,2));  // TCOLLMIN for ZPM

  T_Coll_Mode_2 = TFCOLL * Max(E_Exc_S2,0.) + TCOLLFRAC * (De_Saddle_Scission(pow(Csng(I_Z_CN),2) / pow(Csng(I_A_CN),0.33333),ESHIFTSASCI_coll) - E_tunn);
  T_Coll_Mode_2 = Max(T_Coll_Mode_2,0.0);
  //    T_Pol_Mode_2 = T_Pol_Red * T_Coll_Mode_2
  T_Pol_Mode_2 = T_Pol_Mode_0;
  //    T_asym_Mode_2 = Sqr(T_Coll_Mode_2^2 + TCOLLMIN^2)
  T_asym_Mode_2 = sqrt(pow(T_Coll_Mode_2,2) + 4*pow(TCOLLMIN,2));

  //    Dim As Single T_asym_Mode_2_dyn  // Collective dynamical effect ?  
  //    T_asym_Mode_2_dyn = 0.009 * (I_Z_CN^2/(I_A_CN^0.333333) - 92.0^2/(238^0.333333) )
  //    T_asym_Mode_2 = Sqr(T_asym_Mode_2^2 + T_asym_Mode_2_dyn^2)

  /*    T_Coll_Mode_3 = TFCOLL * Max(E_Exc_S3,0.E0) + _
        TCOLLFRAC * (De_Saddle_Scission(Csng(I_Z_CN)^2 / _ 
        Csng(I_A_CN)^0.33333E0,ESHIFTSASCI_coll) - E_tunn)
        Print 4954, TFCOLL * Max(E_Exc_S3,0.E0),TCOLLFRAC * (De_Saddle_Scission(Csng(I_Z_CN)^2 / _ 
        Csng(I_A_CN)^0.33333E0,ESHIFTSASCI_coll) - E_tunn),  _
        TCOLLFRAC * (De_Saddle_Scission(Csng(I_Z_CN)^2 / _ 
        Csng(I_A_CN)^0.33333E0,ESHIFTSASCI_coll) ),E_Exc_S3, _
        TCOLLFRAC * 0.03 * (De_Saddle_Scission(Csng(I_Z_CN)^2 / _ 
        Csng(I_A_CN)^0.33333E0,ESHIFTSASCI_coll) )^2           
        sleep */  
  //T_coll_Mode_3 = 0.2     // for 239Pu(nth,f)
  //T_coll_Mode_3 = 0.7     // for 252Cf(sf)
  // Fit to 239Pu(nth,f) and 252Cf(sf) ( unexpectedly large variation with Z^2/A^(1/3) )
  T_Coll_Mode_3 = TFCOLL * Max(E_Exc_S3,0.) + TCOLLFRAC * 0.03 * pow(De_Saddle_Scission(pow(Csng(I_Z_CN),2)/pow(Csng(I_A_CN),0.33333),ESHIFTSASCI_coll),2);           
  T_Coll_Mode_3 = Max(T_Coll_Mode_3,0.0);
  //    T_Pol_Mode_3 = T_Pol_Red * T_Coll_Mode_3
  T_Pol_Mode_3 = T_Pol_Mode_0;
  T_asym_Mode_3 = sqrt(pow(T_Coll_Mode_3,2) + pow(TCOLLMIN,2));   //!!!
  //    Dim As Single T_asym_Mode_3_dyn    
  // Adjusted to the width of Mode 3 in 252Cf(sf)
  // May be, this is a collective dynamic effect along the fission path
  //    T_asym_Mode_3_dyn = 0.009 * (I_Z_CN^2/(I_A_CN^0.333333) - 92.0^2/(238^0.333333) )
  // T_asym_Mode_3_dyn = 0   
  //    T_asym_Mode_3 = Max(T_asym_Mode_3,T_asym_Mode_3_dyn)
  //    T_asym_Mode_3 = Sqr(T_asym_Mode_3^2 + T_asym_Mode_3_dyn^2)

  //Print "4619: ";T_Coll_Mode_3,TCOLLMIN,T_asym_Mode_3    
  //sleep

  T_Coll_Mode_4 = TFCOLL * Max(E_Exc_S4,0.) + TCOLLFRAC * (De_Saddle_Scission(pow(Csng(I_Z_CN),2) /  pow(Csng(I_A_CN),0.33333),ESHIFTSASCI_coll) - E_tunn);
  T_Coll_Mode_4 = Max(T_Coll_Mode_4,0.0);
  //    T_Pol_Mode_4 = T_Pol_Red * T_Coll_Mode_4
  T_Pol_Mode_4 = T_Pol_Mode_0;
  T_asym_Mode_4 = sqrt(pow(T_Coll_Mode_4,2) + 4.0*pow(TCOLLMIN,2));  // ZPM like S1
}

void GEF::MeanValues_FM(void)
{
  /* Mean values and standard deviations of fission modes */    
  SigZ_Mode_0 = sqrt(0.5 * T_asym_Mode_0/R_Z_Curv_S0);
  //Print "4214: SIGZ_Mode_0, T_asym_Mode_0, R_Z_Curv_S0)";SigZ_Mode_0,T_asym_Mode_0,R_Z_Curv_S0
  //sleep    
  if(T_Pol_Mode_0 > 1.e-2)
    Sigpol_Mode_0 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S0 / Tanh(HOMPOL/(2. * T_Pol_Mode_0)));
  else
    Sigpol_Mode_0 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S0);
  /* including influence of zero-point motion */

  R_E_intr_S1 = Max(E_Exc_S1+LyPair(I_Z_CN,I_A_CN),0.0);
  R_Att[1] = exp(-R_E_intr_S1/Shell_fading);
  R_Att[5] = R_Att[1];
  R_Att_Sad[1] = exp(-R_E_intr_S1/Shell_fading);
  R_Att_Sad[5] = R_Att_Sad[1];
  SigZ_Mode_1 = sqrt(0.5 * T_asym_Mode_1/(P_Z_Curv_S1*sqrt(R_Att[1])));
  if(T_Pol_Mode_1 > 1.e-2)
    Sigpol_Mode_1 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S1 / Tanh(HOMPOL/(2. * T_Pol_Mode_1)));
  else
    Sigpol_Mode_1 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S1);

  R_E_intr_S2 = Max(E_Exc_S2+LyPair(I_Z_CN,I_A_CN),0.0);
  R_Att[2] = exp(-R_E_intr_S2/Shell_fading);
  R_Att[6] = R_Att[2];
  R_Att_Sad[2] = exp(-R_E_intr_S2/Shell_fading);
  R_Att_Sad[6] = R_Att_Sad[2];
  SigZ_Mode_2 = sqrt(0.5 * T_asym_Mode_2/(P_Z_Curv_S2*sqrt(R_Att[2])));
  SigZ_SL4 = sqrt(0.5 * T_asym_Mode_2/(P_Z_Curv_SL4*sqrt(R_Att[2])));

  if(T_Pol_Mode_2 > 1.e-2)
    Sigpol_Mode_2 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S2 / Tanh(HOMPOL/(2. * T_Pol_Mode_2)));
  else
    Sigpol_Mode_2 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S2);

  R_E_intr_S3 = Max(E_Exc_S3+LyPair(I_Z_CN,I_A_CN),0.0);
  R_Att[3] = exp(-R_E_intr_S3/Shell_fading);
  R_Att_Sad[3] = exp(-R_E_intr_S3/Shell_fading);
  SigZ_Mode_3 = sqrt(0.5 * T_asym_Mode_3/(P_Z_Curv_S3*sqrt(R_Att[3])));
  if(T_Pol_Mode_3 > 1.e-2)
    Sigpol_Mode_3 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S3 / Tanh(HOMPOL/(2. * T_Pol_Mode_3)));
  else
    Sigpol_Mode_3 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S3);

  R_E_intr_S4 = Max(E_Exc_S4+LyPair(I_Z_CN,I_A_CN),0.0);
  R_Att[4] = exp(-R_E_intr_S4/Shell_fading);
  R_Att_Sad[4] = exp(-R_E_intr_S4/Shell_fading);
  SigZ_Mode_4 = sqrt(0.5 * T_asym_Mode_4/(P_Z_Curv_S4*sqrt(R_Att[4])));
  if(T_Pol_Mode_4 > 1.e-2)
    Sigpol_Mode_4 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S4 /Tanh(HOMPOL/(2. * T_Pol_Mode_4)));
  else
    Sigpol_Mode_4 = sqrt(0.25 * HOMPOL / R_Pol_Curv_S4);
}

void GEF::EnergyDependence(void)
{
  /* Energy-dependent shift of fission channels */
  float DZ_S1,DZ_S2,DZ_S3,DZ_S4;
  float EtotS0, EtotS2;
  float P_Z_Curv_S1_eff;
  P_Z_Curv_S1_eff = P_Z_Curv_S1 * P_Z_Curvmod_S1;
  float P_Z_Curv_S2_eff;
  P_Z_Curv_S2_eff = P_Z_Curv_S2 * P_Z_Curvmod_S2;     
  float P_Z_Curv_S3_eff;
  P_Z_Curv_S3_eff = P_Z_Curv_S3 * P_Z_Curvmod_S3;     
  float P_Z_Curv_S4_eff;
  P_Z_Curv_S4_eff = P_Z_Curv_S4 * P_Z_Curvmod_S4;     

  DZ_S1 = P_Z_Mean_S1 * 
    (P_Z_Curv_S1_eff*R_Att[1] / (R_Z_Curv_S0 + P_Z_Curv_S1_eff*R_Att[1]) 
     - (P_Z_Curv_S1_eff / (R_Z_Curv_S0 + P_Z_Curv_S1_eff) ) );
  DZ_S2 =  P_Z_Mean_S2 * 
    (P_Z_Curv_S2_eff*R_Att[2] / (R_Z_Curv_S0 + P_Z_Curv_S2_eff*R_Att[2]) 
     - (P_Z_Curv_S2_eff / (R_Z_Curv_S0 + P_Z_Curv_S2_eff) ) );  
  DZ_S3 =  P_Z_Mean_S3 * 
    (P_Z_Curv_S3_eff*R_Att[3] / (R_Z_Curv_S0 + P_Z_Curv_S3_eff*R_Att[3]) 
     - (P_Z_Curv_S3_eff / (R_Z_Curv_S0 + P_Z_Curv_S3_eff) ) );
  DZ_S4 = Sgn(P_Z_Mean_S4 - P_Z_Mean_S0) * P_Z_Mean_S4 * 
    (P_Z_Curv_S4_eff*R_Att[4] / (R_Z_Curv_S0 + P_Z_Curv_S4_eff*R_Att[4]) 
     - (P_Z_Curv_S4_eff / (R_Z_Curv_S0 + P_Z_Curv_S4_eff) ) ) ; 

  // Empirical shift of S2 channel at low excitation energy at scission 
  // for better reproduction of 238U(s,f) and some data for Th isotopes.
  // Does not solve the problem of 229Th(nth,f).    
  //   EtotS2 = Max(E_Exc_S2 + E_diss_Scission,0.0)
  //   If EtotS2 < 5.E0 Then
  //     DZ_S2 = DZ_S2 + (5.E0 - EtotS2) * 0.1
  //   End If             

  //   DZ_S1 = 0
  //   DZ_S2 = 0
  //   DZ_S3 = 0
  //   DZ_S4 = 0

  ZC_Mode_0 = P_Z_Mean_S0;
  ZC_Mode_1 = P_Z_Mean_S1 + DZ_S1;  
  ZC_Mode_2 = P_Z_Mean_S2 + DZ_S2;  
  ZC_Mode_3 = P_Z_Mean_S3 + DZ_S3;
  //   ZC_Mode_4 = P_Z_Mean_S4 + DZ_S4  
  // shift is very small, because S4 exists only close to symmetry
  ZC_Mode_4 = P_Z_Mean_S4; 

  /* Energy dependence of charge polarization */
  /* Due to washing out of shells */

  for(int I = 10;I<=(I_A_CN - 10);I++)   // mass number
    for(int J = 1;J<=4;J++)    // fission channel
      for(int K = 1;K<=2;K++)    // light - heavy group
        Zshift[J][K][I] = ZshiftOriginal[0][K][I] + (ZshiftOriginal[J][K][I] - ZshiftOriginal[0][K][I])*R_Att[J];

  /* Energy dependence of shell-induced deformation */
  /* Due to washing out of shells */
  /* (Under development) */
  /*For I = 10 To I_Z_CN - 10  // mass number
    For J = 1 To 4           // fission channel
    For K = 1 To 2         // light - heavy group
    beta(J,K,I) = beta(0,K,I) + (beta(J,K,I) - beta(0,K,I))*R_Att_Sad(J)
    if beta(J,K,I) < 0 Then 
    beta(J,K,I) = 0
    End If  
    Z1 = I
    Z2 = I_Z_CN - Z1
    A1 = Z1 / Csng(I_Z_CN) * Csng(I_A_CN)
    A2 = I_A_CN - A1
    E_defo = Lymass(Z1,A1,beta(J,K,I)) - Lymass(Z1,A1,0.0)
    Edefo(J,K,I) = E_defo
    Next
    Next    
    Next  */  
}

void GEF::Yields_FM(void)
{
  //  Yield_Mode_0 = Getyield(E_Exc_S0,E_Exc_S0,T_low_SL,TEgidy(Csng(I_A_CN),0.E0,Tscale))
  Yield_Mode_0 = Getyield(E_Exc_S0,E_Exc_S0,T_low_SL,TEgidy(Csng(I_A_CN),Delta_S0,Tscale));
  Yield_Mode_1 = Getyield(E_Exc_S1,E_Exc_S0,T_low_S1_used,TEgidy(Csng(I_A_CN),R_Shell_S1_eff + dE_Defo_S1,Tscale));
  /*  - Getyield(E_Exc_S0 - E_ld_S1,T_low,T_high) */
  Yield_Mode_2 = Getyield(E_Exc_S2,E_Exc_S0,T_low_S2,TEgidy(Csng(I_A_CN),R_Shell_S2_eff + dE_Defo_S2,Tscale));
  /*  - Getyield(E_Exc_S0 - E_ld_S2,T_low,T_high) */
  Yield_Mode_3 = Getyield(E_Exc_S3,E_Exc_S0,T_low_S3,TEgidy(Csng(I_A_CN),R_Shell_S3_eff + dE_Defo_S3,Tscale));
  /*  - Getyield(E_Exc_S0 - E_ld_S3,T_low,T_high) */
  Yield_Mode_4 = Getyield(E_Exc_S4,E_Exc_S0,T_low_S4,TEgidy(Csng(I_A_CN),R_Shell_S4_eff + dE_Defo_S4,Tscale));  
  /*   - Getyield(E_Exc_S0 - E_ld_S4,T_low,T_high) */ 
  //Print TEgidy(Csng(I_A_CN),0.E0,Tscale), TEgidy(Csng(I_A_CN),R_Shell_S2_eff + dE_Defo_S2,Tscale), de_Defo_S2 
  //sleep
  if(B_S11 > B_S1) 
    Yield_Mode_11 = 0.0;
  else
    Yield_Mode_11 = Getyield(E_Exc_S11,E_Exc_S0, T_low_S11,TEgidy(Csng(I_A_CN),R_Shell_S1_eff + 2.* dE_Defo_S1,Tscale));
  if(B_S22 > B_S2)
    Yield_Mode_22 = 0.0;
  else
    Yield_Mode_22 = Getyield(E_Exc_S22,E_Exc_S0, T_low_S2, TEgidy(Csng(I_A_CN),R_Shell_S2_eff,Tscale));     
  Yield_Norm = Yield_Mode_0 + Yield_Mode_1 + Yield_Mode_2 + Yield_Mode_3 + Yield_Mode_4 + Yield_Mode_11 + Yield_Mode_22;
  Yield_Mode_0 = Yield_Mode_0 / Yield_Norm;
  Yield_Mode_1 = Yield_Mode_1 / Yield_Norm;
  Yield_Mode_2 = Yield_Mode_2 / Yield_Norm;
  Yield_Mode_3 = Yield_Mode_3 / Yield_Norm;
  Yield_Mode_4 = Yield_Mode_4 / Yield_Norm;
  Yield_Mode_11 = Yield_Mode_11 / Yield_Norm;
  Yield_Mode_22 = Yield_Mode_22 / Yield_Norm;
  //cout<< B_S11<<" "<<B_S1<<" "<<E_Exc_S0<<"-- "<<Yield_Mode_0<<" "<<Yield_Mode_1<<" "<<Yield_Mode_2<<" "<<Yield_Mode_3<<" "<<Yield_Mode_4<<" "<<Yield_Mode_11<<" "<<Yield_Mode_22<<endl;
}

void GEF::MassWidths_FM(void)
{
  SigA_Mode_0 = SigZ_Mode_0 * Csng(I_A_CN) / Csng(I_Z_CN); /* width in mass */
  SigA_Mode_1 = SigZ_Mode_1 * Csng(I_A_CN) / Csng(I_Z_CN);
  SigA_Mode_1 = Min(SigA_Mode_1,SigA_Mode_0);  // not broader than liquid-drop
  SigA_Mode_2 = SigZ_Mode_2 * Csng(I_A_CN) / Csng(I_Z_CN);
  SigA_Mode_2 = Min(SigA_Mode_2,SigA_Mode_0);  // not broader than liquid-drop
  SigA_Mode_3 = SigZ_Mode_3 * Csng(I_A_CN) / Csng(I_Z_CN);
  SigA_Mode_3 = Min(SigA_Mode_3,SigA_Mode_0);
  SigA_Mode_4 = SigZ_Mode_4 * Csng(I_A_CN) / Csng(I_Z_CN);
  SigA_Mode_4 = Min(SigA_Mode_4,SigA_Mode_0);
  SigA_Mode_11 = SigZ_Mode_1 * sqrt(2.E0) * Csng(I_A_CN) / Csng(I_Z_CN);
  SigA_Mode_11 = Min(SigA_Mode_11,SigA_Mode_0);
  SigA_Mode_22 = SigZ_Mode_2 * sqrt(2.E0) * Csng(I_A_CN) / Csng(I_Z_CN);
  SigA_Mode_22 = Min(SigA_Mode_22,SigA_Mode_0);
}

void GEF::ShellEff_FM(void)
{
  /* This is the "real" microscopic shell effect, not the effective shell-correction energy */
  /* EShell acts on the level density and determines the T parameter */
  for(int I = 1;I<=(I_A_CN - 1);I++)
  {
    for(int J = 0;J<= 4;J++)
      EShell[J][1][I] = 0;   /* Shells in "light" fragment assumed to be zero */
    DU0 = 0;
    EShell[0][2][I] = 0; /* Shell = 0 in symmetric mode */
    DU1 = R_Shell_S1_eff + dE_Defo_S1; /* + R_A_Curv1_S1 * (AC_Mode_1 - Float(I,6))**2; */
    DU1 = Min(DU1,0.);  /* Technical limit */
    EShell[1][2][I] = DU1;
    DU2 = R_Shell_S2_eff + dE_Defo_S2; /* + R_A_Curv1_S2 * (AC_Mode_2 - Float(I,6))**2; */
    DU2 = Min(DU2,0.);  /* Technical limit */
    EShell[2][2][I] = DU2;
    DU3 = R_Shell_S3_eff + dE_Defo_S3; /* + R_A_Curv1_S3 * (AC_Mode_3 - Float(I,6))**2; */
    DU3 = Min(DU3,0.);  /* Technical limit */
    EShell[3][2][I] = DU3;
    DU4 = R_Shell_S4_eff + dE_Defo_S4; /* + R_A_Curv1_S4 * (AC_Mode_4 - Float(I,6))**2; */
    DU4 = Min(DU4,0.);  /* Technical limit */
    EShell[4][2][I] = DU4;
  }
}

void GEF::IntrinsicT(void)
{
  /* Intrinsic temperatures of fragments at scission */
  /* Mean values */
  T_intr_Mode_0 = TEgidy(AC_Mode_0,0.0,0.8);
  T_intr_Mode_1_heavy = TEgidy(AC_Mode_1,R_Shell_S1_eff + dE_Defo_S1,Tscale);
  T_intr_Mode_1_light = TEgidy(Csng(I_A_CN) - AC_Mode_1,0.0,Tscale);
  T_intr_Mode_2_heavy = TEgidy(AC_Mode_2,R_Shell_S2_eff + dE_Defo_S2,Tscale);
  T_intr_Mode_2_light = TEgidy(Csng(I_A_CN) - AC_Mode_2,0.0,Tscale);
  T_intr_Mode_3_heavy = TEgidy(AC_Mode_3,R_Shell_S3_eff + dE_Defo_S3,Tscale);
  T_intr_Mode_3_light = TEgidy(Csng(I_A_CN) - AC_Mode_3,0.0,Tscale);
  T_intr_Mode_4_heavy = TEgidy(AC_Mode_4,R_Shell_S4_eff + dE_Defo_S4,Tscale);
  T_intr_Mode_4_light = TEgidy(Csng(I_A_CN) - AC_Mode_4,0.0,Tscale);

  /* Mass-dependent values of individual fragments */
  /* Mode 0 */
  for(int I = 1;I<=(I_A_CN - 1);I++)
  {
    T = TEgidy(Csng(I),EShell[0][1][I],Tscale);
    Temp[0][1][I] = T; /* "light" fragment at freeze-out (somewhere before scission) */
    T = TEgidy(Csng(I),EShell[0][2][I],Tscale);
    Temp[0][2][I] = T; /* "heavy" fragment at freeze-out (somewhere before scission) */

    T = TEgidy(Csng(I),0.0,1.0);
    TempFF[0][1][I] = T;       // FF in their ground state
    TempFF[0][2][I] = T;       // FF in their ground state 
  }
  /* Mode 1 */
  for(int I = 1;I<=(I_A_CN - 1);I++)
  {
    T = TEgidy(Csng(I),EShell[1][1][I],Tscale);
    Temp[1][1][I] = T;  /* "light" fragment */
    T = TEgidy(Csng(I),EShell[1][2][I],Tscale);
    Temp[1][2][I] = T;  /* "heavy" fragment */

    T = TEgidy(Csng(I),0.0,1.0);
    TempFF[1][1][I] = T;       // FF in their ground state
    TempFF[1][2][I] = T;       // FF in their ground state
  }
  /* Mode 2 */
  for(int I = 1;I<=(I_A_CN - 1);I++)
  {
    T = TEgidy(Csng(I),EShell[2][1][I],Tscale);
    Temp[2][1][I] = T; /* "light" fragment */
    T = TEgidy(Csng(I),EShell[2][2][I],Tscale);
    Temp[2][2][I] = T; /* "heavy" fragment */
    /* The next section is introduced, because energy sorting is not strong enough,
       when shells are only introduced in the heavy fragment.
       Ad hoc assumption: For Mode 2 there are shells in both fragments of about
       equal size. Technically, we neglect the shells in both fragments.
       This has about the same effect for the energy sorting. */
    T = TEgidy(Csng(I),0.0,Tscale);   // FF at scssion
    Temp[2][1][I] = T; /* "light" fragment */
    T = TEgidy(Csng(I),0.0,Tscale);   // FF at scission
    Temp[2][2][I] = T; /* "heavy" fragment */

    T = TEgidy(Csng(I),0.0,1.0);    // shell effect neglected
    TempFF[2][1][I] = T;    // FFs in their ground state
    TempFF[2][2][I] = T;    // FFs in their ground state
  }
  /* Mode 3 */
  for(int I = 1;I<=(I_A_CN -1);I++)
  {
    T = TEgidy(Csng(I),0.0,Tscale);
    Temp[3][1][I] = T;
    T = TEgidy(Csng(I),0.0,Tscale);
    Temp[3][2][I] = T;

    T = TEgidy(Csng(I),0.0,1.0);
    TempFF[3][1][I] = T;       // FF in their ground state
    TempFF[3][2][I] = T;       // FF in their ground state
  }
  /* Mode 4 */
  for(int I = 1;I<=(I_A_CN -1);I++)
  {
    T = TEgidy(Csng(I),0.0,Tscale);
    Temp[4][1][I] = T;
    T = TEgidy(Csng(I),0.0,Tscale);
    Temp[4][2][I] = T;

    T = TEgidy(Csng(I),0.0,1.0);
    TempFF[4][1][I] = T;       // FF in their ground state
    TempFF[4][2][I] = T ;      // FF in their ground state
  } 
}

void GEF::IntrinsicEx(void)
{
  /*** Intrinsic excitation energy at saddle and at scission as well as   ***/
  /*** Even-odd effect in proton and neutron number for each fission mode ***/
  for(I_Mode = 0; I_Mode<=6;I_Mode++)
  {
    E_coll_saddle[I_Mode] = 0;
    if(I_Mode == 0) Etot = E_Exc_S0;
    if(I_Mode == 1) Etot = E_Exc_S1;
    if(I_Mode == 2) Etot = E_Exc_S2;
    if(I_Mode == 3) Etot = E_Exc_S3;
    if(I_Mode == 4) Etot = E_Exc_S4;
    if(I_Mode == 5) Etot = E_Exc_S11;
    if(I_Mode == 6) Etot = E_Exc_S22;

    if(( Mod(I_Z_CN,2) + Mod(I_N_CN,2)) == 0)   /* Even-even CN */      
      if(Etot > 0 && Etot < (2. * 14./sqrt(Csng(I_A_CN))) )
      {
        E_coll_saddle[I_Mode] = Etot;
        Etot = 0;
        /* Excitation below the pairing gap in even-even CN goes into collective excitations */
      }
    //    If I_Z_CN Mod 2 + I_N_CN Mod 2 = 0 Then    // even-even
    //      Ediff = Min(Etot, 14.0/sqr(Csng(I_A_CN)))
    //    End If
    //    If I_Z_CN Mod 2 + I_N_CN Mod 2 = 1 Then    // even-odd or odd-even
    //       Ediff = Min(Etot, 2.0 * 14.0/sqr(Csng(I_A_CN)))
    //    End If
    //    Ediff = Max(Ediff,0.0) 
    //    Etot = Etot - Ediff
    if(Etot < 0)
      E_tunn = -Etot;
    else
      E_tunn = 0;
    Etot = Max(Etot,0.0);
    E_POT_SCISSION = (De_Saddle_Scission(pow(Csng(I_Z_CN),2) / pow(Csng(I_A_CN),0.33333),ESHIFTSASCI_intr) ) ;
    E_diss_Scission = EDISSFRAC * (E_POT_SCISSION - E_tunn) + Epot_shift ; 
    if(E_diss_Scission < 0)
      E_diss_Scission = 0;
    Etot = Etot + E_diss_Scission;
    /* All excitation energy at saddle and part of the potential-energy gain to scission
       go into intrinsic excitation energy at scission */
    if(I_Mode == 2)
      EINTR_SCISSION = Etot; /* (For Mode 2) Global parameter */

    for(IA1 = 40;IA1<=(I_A_CN - 40);IA1++)
    {
      IA2 = I_A_CN - IA1;
      if(I_Mode <= 4)
      {
        T1 = Temp[I_Mode][1][IA1];
        T2 = Temp[I_Mode][2][IA2];
      }
      if(I_Mode == 5)
      {
        T1 = Temp[1][2][IA1];
        T2 = Temp[1][2][IA2];
      }
      if(I_Mode == 6)
      {
        T1 = Temp[2][2][IA1];
        T2 = Temp[2][2][IA2];
      }
      DT = abs(T2 - T1);

      /* Even-odd effect */
      if(Mod(I_Z_CN,2) == 0)
        Rincr1P = exp(-Etot/PZ_EO_symm);
      else
        Rincr1P = 0;  
      if( Mod(I_N_CN,2) == 0)
        Rincr1N = exp(-Etot/PN_EO_Symm);
      else
        Rincr1N = 0;

      PEOZ[I_Mode][1][IA1] = Rincr1P;
      PEOZ[I_Mode][2][IA2] = Rincr1P;
      PEON[I_Mode][1][IA1] = Rincr1N;
      PEON[I_Mode][2][IA2] = Rincr1N;

      if(Etot > 0)
        Rincr2 = Gaussintegral(DT/Etot-R_EO_THRESH,R_EO_SIGMA*(DT+0.0001));
      /* even-odd effect due to asymmetry */
      else   // even-odd effect is already 100%
        Rincr2 = 0;            
      Rincr2P = (R_EO_MAX - Rincr1P) * Rincr2;
      Rincr2N = (R_EO_MAX - Rincr1N) * Rincr2;     

      if(IA1 <= IA2)  // A1 is lighter or equal to A2
      {
        PEOZ[I_Mode][1][IA1] = PEOZ[I_Mode][1][IA1] + Rincr2P;
        if(Mod(I_Z_CN,2) == 0)
          PEOZ[I_Mode][2][IA2] = PEOZ[I_Mode][2][IA2] + Rincr2P;
        else
          PEOZ[I_Mode][2][IA2] = PEOZ[I_Mode][2][IA2] - Rincr2P;
        PEON[I_Mode][1][IA1] = PEON[I_Mode][1][IA1] + Rincr2N;
        if(Mod(I_N_CN,2) == 0)
          PEON[I_Mode][2][IA2] = PEON[I_Mode][2][IA2] + Rincr2N;
        else
          PEON[I_Mode][2][IA2] = PEON[I_Mode][2][IA2] - Rincr2N;
      }
      else
      {
        PEOZ[I_Mode][1][IA1] = PEOZ[I_Mode][2][IA1];
        PEON[I_Mode][1][IA1] = PEON[I_Mode][2][IA1];
        PEOZ[I_Mode][2][IA2] = PEOZ[I_Mode][1][IA2];
        PEON[I_Mode][2][IA2] = PEON[I_Mode][1][IA2];
      }            
      /*  Else
          PEOZ(I_Mode,2,IA2) = _
          PEOZ(I_Mode,1,IA2) + Rincr2P
          IF I_Z_CN Mod 2 = 0 Then
          PEOZ(I_Mode,1,IA1) = _
          PEOZ(I_Mode,1,IA1) + Rincr2P
          Else
          PEOZ(I_Mode,1,IA1) = _
          PEOZ(I_Mode,1,IA1) - Rincr2P
          End if
          PEON(I_Mode,2,IA2) = _
          PEON(I_Mode,2,IA2) + Rincr2N
          IF I_N_CN Mod 1 = 0 Then
          PEON(I_Mode,1,IA1) = _
          PEON(I_Mode,1,IA1) + Rincr2N
          Else
          PEON(I_Mode,1,IA1) = _
          PEON(I_Mode,1,IA1) - Rincr2N
          End if
          End If  */
      PEOZ[I_Mode][1][IA1] = PEOZ[I_Mode][1][IA1] * EOscale;
      PEOZ[I_Mode][2][IA2] = PEOZ[I_Mode][2][IA2] * EOscale;
      PEON[I_Mode][1][IA1] = PEON[I_Mode][1][IA1] * EOscale;
      PEON[I_Mode][2][IA2] = PEON[I_Mode][2][IA2] * EOscale;

      /*  For T1 = 0.5 To 1 Step 0.05
          T2 = 1 - T1
          Print T1,T2, 1.E - Gaussintegral(T2-T1,0.1)
          Next T1
          Sleep  */

      /* Energy sorting */
      /* E1 = Etot * Gaussintegral(T2-T1,0.03); */
      if(abs(T1-T2) < 1.e-6)
        E1 = 0.5E0 * Etot;
      else
      {
        //   E1ES = Csort * T1 * T2 / ( Abs(T1 - T2) )
        if(I_Mode == 0)
          E1ES = Etot * Gaussintegral(T2-T1,Esort_slope_S0);
        else
          E1ES = Etot * Gaussintegral(T2-T1,Esort_slope);
        E1ES = Min(E1ES,0.5E0*Etot);
        /* Asymptotic value after "complete" energy sorting */
        E1FG = Etot * IA1 / I_A_CN;  /* in Fermi-gas regime */

        /*  The following section assumes T=const in the superfluid regime
            If Etot < 13 Then E1 = E1ES  // complete energy sorting
            If Etot >= 13 and Etot <= 20 Then  // transition region
            E1 = E1ES + (Etot-13)/7*(E1FG-E1ES)
            End If
            If Etot > 20 Then E1 = E1FG   // Fermi-gas regime   */

        /* The following section extends energy sorting to higher energies */
        if(Etot < 13 * Esort_extend)       // complete energy sorting
          E1 = E1ES; 
        else if(Etot >= 13 * Esort_extend && Etot <= 20 * Esort_extend)  // transition region
          E1 = E1ES + (Etot-13 * Esort_extend)/7*(E1FG-E1ES);
        else if(Etot > 20 * Esort_extend)
          E1 = E1FG;   // Fermi-gas regime   */
      }
      E2 = Etot - E1;
      EPART[I_Mode][1][IA1] = Max(E1,0.0);  /* Mean E* in light fragment */
      EPART[I_Mode][2][IA2] = Max(E2,0.0);  /* Mean E* in heavy fragment */
    }
  }
}

void GEF::FF_angularMomentum(void)
{
  // RMS angular momentum of fission fragments 
  // Following Naik et al., EPJ A 31 (2007) 195 and  
  // S. G. Kadmensky, Phys. At. Nucl. 71 (2008) 1193  

  float AUCD;   // UCD fragment mass 
  float I_rigid_spher;  // I rigid for spherical shape 
  float I_rigid;        // I rigid for deformed scission shape 
  float I_eff;          // I with reduction due to pairing 
  float alph;           // deformation parameter 
  float E_exc;          // Excitation energy 
  float J_rms;          // rms angular momentum 
  int ZT,AT;         // Z and A of target nucleus 
  float I_ISO;          // ISO number 

  // CN spin 
  ZT = P_Z_CN;
  //  AT = I_A_CN
  AT = P_A_CN;
  if(Emode == 2)
    AT = P_A_CN -1;
  if(Emode == 12)
  {
    AT = P_A_CN -1;
    ZT = P_Z_CN -1;
  }

  //    Spin_CN = P_J_CN
  //    P_I_rms_CN = P_J_CN

  int I_E_iso=0;  //Only Ground STATE
  if(I_E_iso == 0) // fissioning nucleus or target in ground state
  {
    if(P_I_rms_CN == 0)
    {
      if(AT == 0)
        Spin_CN = 0;
      else  
        Spin_CN =  R_SPI[ZT][AT];
    }
    else
      Spin_CN = P_I_rms_CN;
  }
  /*
     else // fissioning nucleus or target in isomeric state
     {
     if (N_ISO_MAT(I_MAT) < I_E_iso) //??????????????????
     {
     cout<<"The isomer is not in the table of nuclear properties."<<endl;
     cout<<"Z, A, #iso "<<ZT<<" "<<AT<<" "<<I_E_iso<<endl;
     cout<< "Please restart GEF."<<endl;
     }
     Spin_CN = NucTab(I_MAT + I_E_iso).R_SPI;     //??????????????????
     }*/

  //   Print "ZT, AT, I_MAT",ZT, AT, I_MAT
  //   Print "SPIN_CN",Spin_CN
  //   Sleep 
  Spin_pre_fission = Spin_CN;  // target or CN ground-state spin
  // Incoming neutron (spin + orbital) 
  if(Emode == 2)
    // 2/3 * 1.16 * sqr(2 * 939.65) / 197.33 = 0.1699 
    if(R_E_exc_used > 0) 
      Spin_pre_fission = sqrt(pow(Spin_CN,2) + pow(0.5,2) + pow(0.1699 * pow(AT,0.333333) * sqrt(R_E_exc_used),2));
  if(Emode == 12)  // preliminary (because Coulomb barrier neglected) 
    if(R_E_exc_used > 0)
      Spin_pre_fission = sqrt(pow(Spin_CN,2) + pow(0.5,2) + pow(0.1699 * pow(AT,0.333333) * sqrt(R_E_exc_used),2)); 

  for(IZ1 = 10;IZ1<=(I_Z_CN - 10);IZ1++)
  {
    AUCD = CInt(Csng(IZ1) * Csng(I_A_CN) / Csng(I_Z_CN));
    for(IA1 = CInt(AUCD - 15);IA1<= CInt(AUCD + 15);IA1++)
    {
      IN1 = IA1 - IZ1;
      if(IA1 - IZ1 >= 10)
      {
        // Rigid momentum of inertia for spherical nucleus 
        I_rigid_spher = pow(1.16,2) * pow(Csng(IA1),1.6667) / 103.8415;
        // unit: hbar^2/MeV 
        for(I_Mode = 0;I_Mode<=6;I_Mode++)  
        {
          // First (normally light) fission fragment: 
          beta1 = Beta[I_Mode+1][1][IZ1]; //I_Mode+1 because there is I_Mode = -1 for Beta[][][]
          alph = beta1 / sqrt(4. * pi / 5.);
          I_rigid = I_rigid_spher * (1. + 0.5*alph + 9./7.*pow(alph,2));
          // From Hasse & Myers, Geometrical Relationships ... 
          E_exc = EPART[I_Mode][1][IA1];
          if(E_exc < 0)
            E_exc = 0;
          T = U_Temp(Csng(IZ1),Csng(IA1),E_exc,1,1,Tscale,Econd);          
          //   T = sqr(T^2 + 0.8^2)       // For ZPM
          //   T = T_orbital
          //   T =  sqr(T^2 + T_orbital^2)
          if(T_orbital > 0.1)
            T = T_orbital / Tanh(T_orbital/T);  // T_orbital represents the ZPM

          I_eff = I_rigid * (1. - 0.8 * exp(-0.693 * E_exc / 5.));
          J_rms = sqrt(2. * I_eff * T);  

          J_rms = J_rms * Jscaling; 

          if((Mod(IZ1,2) == 1) | (Mod(IN1,2) == 1) )
            J_rms = J_rms + Spin_odd * pow(Csng(IA1)/140.0,0.66667);
          //                 Max(0,1 - (E_exc-1)/9) // empirical 
          // Additional angular momentum of unpaired proton. 
          // See also Tomar et al., Pramana 68 (2007) 111 

          // Print Z1,I_Mode,beta1,T,E_exc,Spin_CN         
          // Print " ",I_rigid_spher,I_rigid,I_eff,J_rms

          J_rms = sqrt(pow(J_rms,2) + pow(IA1/I_A_CN * Spin_pre_fission,2));

          SpinRMSNZ[I_Mode][1][IA1-IZ1][IZ1] = J_rms;

          //     Print
          //     Print IA1,T,E_exc,I_rigid_spher,I_rigid,I_eff,J_rms

          // Second (normally heavy) fission fragment: 

          beta2 = Beta[I_Mode+1][2][IZ1];
          alph = beta2 / sqrt(4. * pi / 5.);
          I_rigid = I_rigid_spher * (1. + 0.5*alph + 9./7.*pow(alph,2));
          // From Hasse & Myers, Geometrical Relationships ... 
          E_exc = EPART[I_Mode][2][IA1];
          if(E_exc < 0)
            E_exc = 0;
          T = U_Temp(Csng(IZ1),Csng(IA1),E_exc,1,1,Tscale,Econd);          
          //    T = sqr(T^2 + 0.8^2)       // For ZPM
          //    T = T_orbital
          //    T =  sqr(T^2 + T_orbital^2)
          if(T_orbital > 0.1)
            T = T_orbital / Tanh(T_orbital/T);  // T_orbital represents the ZPM
          I_eff = I_rigid * (1. - 0.8 * exp(-0.693 * E_exc / 5.));
          J_rms = sqrt(2. * I_eff * T);

          J_rms = J_rms * Jscaling; 

          if((Mod(IZ1,2) == 1) | (Mod(IN1,2)== 1)) 
            J_rms = J_rms + Spin_odd * pow(Csng(IA1)/140.0,0.66667);  
          //                  Max(0,1 - (E_exc-1)/9) // empirical 
          // Additional angular momentum of unpaired proton. 
          // See also Tomar et al., Pramana 68 (2007) 111 

          J_rms = sqrt(pow(J_rms,2) + pow(IA1/I_A_CN * Spin_pre_fission,2));

          SpinRMSNZ[I_Mode][2][IA1-IZ1][IZ1] = J_rms;

          //      Print IA1,T,E_exc,I_rigid_spher,I_rigid,I_eff,J_rms          

        }
      }
    }
  }  
}


void GEF::FissionModeChoise(void)
{
  /* Chosing fission mode*/
  float R_Sum_0,R_Sum_1,R_Sum_2,R_Sum_3,R_Sum_4,R_Sum_5,R_Sum_6;
  float R_Choice,Rincr;
  R_Choice = rn->Rndm();
  R_Sum_0 = Yield_Mode_0;
  R_Sum_1 = R_Sum_0 + Yield_Mode_1;
  R_Sum_2 = R_Sum_1 + Yield_Mode_2;
  R_Sum_3 = R_Sum_2 + Yield_Mode_3;
  R_Sum_4 = R_Sum_3 + Yield_Mode_4;
  R_Sum_5 = R_Sum_4 + Yield_Mode_11;
  R_Sum_6 = R_Sum_5 + Yield_Mode_22;   
  I_Mode = 6;
  if(R_Choice < R_Sum_0) I_Mode = 0;
  if(R_Choice >= R_Sum_0 && R_Choice < R_Sum_1) I_Mode = 1;
  if(R_Choice >= R_Sum_1 && R_Choice < R_Sum_2) I_Mode = 2;
  if(R_Choice >= R_Sum_2 && R_Choice < R_Sum_3) I_Mode = 3;
  if(R_Choice >= R_Sum_3 && R_Choice < R_Sum_4) I_Mode = 4;
  if(R_Choice >= R_Sum_4 && R_Choice < R_Sum_5) I_Mode = 5;
  if(R_Choice >= R_Sum_5 && R_Choice < R_Sum_6) I_Mode = 6;
  // Mode_Events(I_Mode) = Mode_Events(I_Mode) + 1  
  //	Mode_Events(10) = Mode_Events(10) + 1  //For normalization
}

void GEF::ZAChoise(void)
{ 
  bool DiceA=true;
  while(DiceA)
  {
    DiceA = false;
    switch(I_Mode)
    {
      case 0:
        {
          R_A_heavy = PGauss(AC_Mode_0,SigA_Mode_0); /* random choice of mass */
          RZpol = Zshift[0][2][CInt(R_A_heavy)]; /* local polarization */
          RZ = R_A_heavy * I_Z_CN / I_A_CN + RZpol; /* mean position in Z for given mass */
          R_Z_Heavy = PGauss(RZ,Sigpol_Mode_0); /* random choice of Z */
        }break;
      case 1:
        {
          R_A_heavy = PGauss(AC_Mode_1,SigA_Mode_1);
          RZpol = Zshift[1][2][CInt(R_A_heavy)]; /* local polarization */
          RZ = R_A_heavy * I_Z_CN / I_A_CN + RZpol;          
          R_Z_Heavy = PGauss(RZ,Sigpol_Mode_1);
        }break;
      case 2:
        {
          Iguess = 0;
          bool Next2 = true;
          while(Next2)
          {
            Next2=false;
            Iguess = Iguess + 1;

            R_A_heavy = PBox2(AC_Mode_2,SigA_Mode_2*S2leftmod,SigA_Mode_2,P_A_Width_S2);
            RZpol = Zshift[2][2][CInt(R_A_heavy)];
            RZ = R_A_heavy * I_Z_CN / I_A_CN + RZpol    ;   
            RN = R_A_heavy - RZ;   
            //    if(Iguess < 1.e3)
            //      if(Rtest > Gaussintegral(RN-82,1.5*SigZ_Mode_2))
            //        Next2 = true;
            //  if (Iguess < 1.e3)
            //   if (Rtest > Gaussintegral(RZ-ZTRUNC50,FTRUNC50*SigZ_Mode_2)) 
            //    if ( (Rtest > 0.5 * ERF((RZ-ZTRUNC50)/(FTRUNC50*SigZ_Mode_2)) + 0.5) |( Rtest > 0.5 * ERF((I_Z_CN - RZ - ZTRUNC28)/(FTRUNC28*SigZ_Mode_2)) + 0.5)) 
            //      Next2 = true;
          }
          /* truncation below Z = 35 and below Z = 50 due to properties of deformed shells */
          R_Z_Heavy = PGauss(RZ,Sigpol_Mode_2);
        }break;
      case 3:
        {
          R_A_heavy = PGauss(AC_Mode_3,SigA_Mode_3);
          RZpol = Zshift[3][2][CInt(R_A_heavy)];
          RZ = R_A_heavy * I_Z_CN / I_A_CN + RZpol;  
          R_Z_Heavy = PGauss(RZ,Sigpol_Mode_3);
        }break;
      case 4:
        {
          if(ZC_Mode_4 > ZC_Mode_0)
            R_A_heavy = PGauss(AC_Mode_4,SigA_Mode_4);
          else  
            R_A_heavy = I_A_CN - PGauss(AC_Mode_4,SigA_Mode_4);
          //AC_Mode_4 refers to the light fragment
          RZpol = Zshift[4][2][CInt(R_A_heavy)];
          RZ = R_A_heavy * I_Z_CN / I_A_CN + RZpol;  
          R_Z_Heavy = PGauss(RZ,Sigpol_Mode_4);
        }break;
      case 5:
        {
          R_A_heavy = PGauss(AC_Mode_0,SigA_Mode_11);
          RZ = R_A_heavy * I_Z_CN / I_A_CN;
          R_Z_Heavy = PGauss(RZ,Sigpol_Mode_0);
        }break;
      case 6:
        {
          R_A_heavy = PGauss(AC_Mode_0,SigA_Mode_22);
          RZ = R_A_heavy * I_Z_CN / I_A_CN;
          R_Z_Heavy = PGauss(RZ,Sigpol_Mode_0);
        }break;
      default:
        {
          cout<< "Wrong I_Mode"<<endl;
        }break;
    }

    R_Z_Light = I_Z_CN - R_Z_Heavy;
    R_A_light = I_A_CN - R_A_heavy;

    if((R_A_heavy < 1) | (R_A_light < 1))
      DiceA = true;
    if(R_A_heavy < R_A_light)
      DiceA = true;
  }
  /*   Diego --- They are only for spectra proposes
  // Pre-neutron Z distribution, without even-odd effect 

  I = CInt(R_Z_Heavy)
  If I >= Lbound(ZPROV) and I <= Ubound(Zprov) Then
  ZPROV(I) = ZPROV(I) + Racc
  ZMPROV(I_Mode,I) = ZMPROV(I_Mode,I) + Racc
  End If 
  I  = CInt(R_Z_Light)
  If I >= Lbound(ZPROV) and I <= Ubound(Zprov) Then 
  ZPROV(I) = ZPROV(I) + Racc
  ZMPROV(I_Mode,I) = ZMPROV(I_Mode,I) + Racc
  End If

  // Provisional mass distribution, pre-neutron,
  //  directly deduced from Z distribution with UCD assumption 

  APROV(CInt(R_Z_Heavy * I_A_CN / I_Z_CN + 0.5)) = _
  APROV(CInt(R_Z_Heavy * I_A_CN / I_Z_CN + 0.5)) + Racc
  APROV(CInt(R_Z_Light * I_A_CN / I_Z_CN + 0.5)) = _
  APROV(CInt(R_Z_Light * I_A_CN / I_Z_CN + 0.5)) + Racc
  AMPROV(I_Mode,CInt(R_Z_Heavy * I_A_CN / I_Z_CN + 0.5)) = _
  AMPROV(I_Mode,CInt(R_Z_Heavy * I_A_CN / I_Z_CN + 0.5)) + Racc
  AMPROV(I_Mode,CInt(R_Z_Light * I_A_CN / I_Z_CN + 0.5)) = _
  AMPROV(I_Mode,CInt(R_Z_Light * I_A_CN / I_Z_CN + 0.5)) + Racc
   */

  /* Round on integer values with even-odd effect */

  R_N_heavy = R_A_heavy - R_Z_Heavy;

  I_N_heavy_sad = EVEN_ODD(R_N_heavy,PEON[I_Mode][2][CInt(R_A_heavy)]);
  I_Z_heavy_sad = EVEN_ODD(R_Z_Heavy,PEOZ[I_Mode][2][CInt(R_A_heavy)]);
  I_A_heavy_sad = I_N_heavy_sad + I_Z_heavy_sad;

  I_N_light_sad = I_N_CN - I_N_heavy_sad;
  I_Z_light_sad = I_Z_CN - I_Z_heavy_sad;
  I_A_light_sad = I_A_CN - I_A_heavy_sad;

  I_A_sad = I_A_CN;
  I_Z_sad = I_Z_CN;
  I_N_sad = I_N_CN;

  Qvalue = AME2012(I_Z_sad,I_A_sad) - (AME2012(I_Z_heavy_sad,I_A_heavy_sad) +  AME2012(I_Z_light_sad,I_A_light_sad));

  E_intr_light = EPART[I_Mode][1][I_A_light_sad];
  E_intr_heavy = EPART[I_Mode][2][I_A_heavy_sad];

}

void GEF::N_Saddle_Scission(void)
{
  float E_Final_Light,J_Frag_light,TlightFF,R_Z_light_sci,R_A_light_sci;
  float E_Final_Heavy,J_Frag_heavy,TheavyFF,R_Z_heavy_sci,R_A_heavy_sci;
  float Array_E_n1_ss[51];
  float Array_E_n2_ss[51];

  TlightFF = 1;  //not used
  TheavyFF = 1;  //not used

  I_nu_ss = 0;

  if(E_intr_light + E_intr_heavy > Escission_lim)
  {
    for(int i = 1; i<=50;i++)
    {
      Array_E_n1_ss[i] = 0;
      Array_Tn[i] = 0;
      Array_E_n1_frag1[i] = 0;
    }
    E_Final_Light = Escission_lim * Csng(I_A_light_sad) / Csng(I_A_sad) - 9;
    J_Frag_light = 0;
  //  Eva(0,I_Z_light_sad,I_A_light_sad,E_intr_light,TlightFF,J_Frag_light,R_Z_light_sci,
  //      R_A_light_sci,E_Final_Light,Array_E_n1_ss,Array_Tn,Array_Eg0_light);             
    Eva(1,I_Z_light_sad,I_A_light_sad,E_intr_light,TlightFF,J_Frag_light,R_Z_light_sci,
        R_A_light_sci,E_Final_Light,Array_E_n1_frag1,Array_Tn,Array_Eg0_light);             
 
    for(int i = 1;i<= 50;i++)
    {
      //cout << Array_E_n1_frag1[i] << endl;
      if(Array_E_n1_ss[i] == 0)
        break;
      I_nu_ss = I_nu_ss + 1;
    }
    I_A_light_sci = R_A_light_sci;
    I_Z_light_sci = R_Z_light_sci;
    I_N_light_sci = I_A_light_sci - I_Z_light_sci;
    E_intr_light =  E_Final_Light;

    for(int i = 1;i<=50;i++)
    {
      Array_E_n2_ss[i] = 0;
      Array_Tn[i] = 0;
      Array_E_n2_frag2[i] = 0;
    }
    E_Final_Heavy = Escission_lim * Csng(I_A_heavy_sad) / Csng(I_A_sad) - 9;
    J_Frag_heavy = 0;
  //  Eva(0,I_Z_heavy_sad,I_A_heavy_sad,E_intr_heavy,TheavyFF,J_Frag_heavy,R_Z_heavy_sci,
  //      R_A_heavy_sci,E_Final_Heavy,Array_E_n2_ss,Array_Tn,Array_Eg0_heavy);                   
    Eva(1,I_Z_heavy_sad,I_A_heavy_sad,E_intr_heavy,TheavyFF,J_Frag_heavy,R_Z_heavy_sci,
        R_A_heavy_sci,E_Final_Heavy,Array_E_n2_frag2,Array_Tn,Array_Eg0_heavy);                   
 
    for(int i = 1;i<=50;i++)
    {
      if(Array_E_n2_ss[i] == 0)
        break;
      I_nu_ss = I_nu_ss + 1;
    }
    I_A_heavy_sci = R_A_heavy_sci;
    I_Z_heavy_sci = R_Z_heavy_sci;
    I_N_heavy_sci = I_A_heavy_sci - I_Z_heavy_sci;
    E_intr_heavy =  E_Final_Heavy;

    I_A_sci = I_A_light_sci + I_A_heavy_sci;
    I_Z_sci = I_Z_light_sci + I_Z_heavy_sci;
    I_N_sci = I_N_light_sci + I_N_heavy_sci;
  }
  else
  {
    I_A_light_sci = I_A_light_sad;
    I_Z_light_sci = I_Z_light_sad;
    I_N_light_sci = I_N_light_sad;
    I_A_heavy_sci = I_A_heavy_sad;
    I_Z_heavy_sci = I_Z_heavy_sad;
    I_N_heavy_sci = I_N_heavy_sad;
    I_A_sci = I_A_light_sad + I_A_heavy_sad;
    I_Z_sci = I_Z_light_sad + I_Z_heavy_sad;
    I_N_sci = I_N_light_sad + I_N_heavy_sad;
  }
}

void GEF::FragmentsEnergy(void)
{
  if(I_Mode == 0)
  {
    // Transition from shell-defined to macroscopic shape with increasing E* 
    // RW_mac = relative macroscopic influence      

    float EtotS0;
    EtotS0 = Max(E_Exc_S0 + E_diss_Scission,0.0);
    RW_mac = Gaussintegral(EtotS0-20.0,5.0);
    /* Only deformation energy: */
    Eexc_light_mean = (1.0 - RW_mac) * Edefo[0][1][I_Z_light_sci] + RW_mac * Edefo[1][1][I_Z_light_sci];
    Eexc_heavy_mean = (1.0 - RW_mac) * Edefo[0][2][I_Z_heavy_sci] + RW_mac * Edefo[1][2][I_Z_heavy_sci];

    RS = SIGDEFO_0 ;
    //     RS = SIGDEFO_0 * Sqr(T_Pol_Mode_0 / TEgidy(P_A_CN,0.0,1.0))
    // Scaling with Sqr(intrinsic temperature / const. T value)    
    ESIGDEFOlight = 
      ( LyMass(I_Z_light_sci,I_A_light_sci,Beta[I_Mode+1][1][I_Z_light_sci] + RS) - 
        LyMass(I_Z_light_sci,I_A_light_sci,Beta[I_Mode+1][1][I_Z_light_sci] ));    
    ESIGDEFOheavy = 
      ( LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[I_Mode+1][2][I_Z_heavy_sci] + RS) - 
        LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[I_Mode+1][2][I_Z_heavy_sci] ));
  }
  if(I_Mode > 0 && I_Mode <= 4)
  {
    Eexc_light_mean = Edefo[I_Mode+1][1][I_Z_light_sci]; /* Only deformation energy */
    Eexc_heavy_mean = Edefo[I_Mode+1][2][I_Z_heavy_sci]; /* Only deformation energy */
    //    RS = SIGDEFO/Sqr(R_Att_Sad(I_Mode))
    RS = SIGDEFO;
    ESIGDEFOlight = 
      ( LyMass(I_Z_light_sci,I_A_light_sci,Beta[I_Mode+1][1][I_Z_light_sci] + RS) - 
        LyMass(I_Z_light_sci,I_A_light_sci,Beta[I_Mode+1][1][I_Z_light_sci] )); 
    ESIGDEFOheavy = 
      ( LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[I_Mode+1][2][I_Z_heavy_sci] + RS) - 
        LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[I_Mode+1][2][I_Z_heavy_sci] )); 

    // If beta(I_MOde,1,I_Z_light) > 0.68 Then
    //   ESIGDEFOlight = 0.4*ESIGDEFOlight
    // End If              
  }
  if(I_Mode == 5)
  {
    Eexc_heavy_mean = Edefo[2][2][I_Z_heavy_sci];
    Eexc_light_mean = Edefo[2][2][I_Z_light_sci];  /* Shell effect stored for "heavy" fragment */
    //    RS = SIGDEFO/Sqr(R_Att_Sad(I_Mode))
    RS = SIGDEFO;
    ESIGDEFOheavy = 
      ( LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[2][2][I_Z_heavy_sci] + RS) - 
        LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[2][2][I_Z_heavy_sci] ));
    ESIGDEFOlight = 
      ( LyMass(I_Z_light_sci,I_A_light_sci,Beta[3][2][I_Z_light_sci] + RS) - 
        LyMass(I_Z_light_sci,I_A_light_sci,Beta[3][2][I_Z_light_sci] ));
  }
  if(I_Mode == 6)
  {
    Eexc_heavy_mean = Edefo[3][2][I_Z_heavy_sci];
    Eexc_light_mean = Edefo[3][2][I_Z_light_sci];  /* Shell effect stored for "heavy" fragment */
    //     RS = SIGDEFO/Sqr(R_Att_Sad(I_Mode))
    RS = SIGDEFO;
    ESIGDEFOheavy = 
      ( LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[3][2][I_Z_heavy_sci] + RS) - 
        LyMass(I_Z_heavy_sci,I_A_heavy_sci,Beta[3][2][I_Z_heavy_sci] ));
    ESIGDEFOlight = 
      ( LyMass(I_Z_light_sci,I_A_light_sci,Beta[3][2][I_Z_light_sci] + RS) - 
        LyMass(I_Z_light_sci,I_A_light_sci,Beta[3][2][I_Z_light_sci] ));
  }
  if(Eexc_heavy_mean < 0)
    Eexc_heavy_mean = 0;
  if(Eexc_light_mean < 0)
    Eexc_light_mean = 0;

  /* Deformation of heavy and light fragment assumed to be uncorrelated */

  TKEmin = 1.44 * I_Z_light_sci * I_Z_heavy_sci / (3.0 * (pow(I_A_light_sci,0.33333) + pow(I_A_heavy_sci,0.3333)));
  /* TKEmin for limiting the excitation energies of the fragments */
  if(Qvalue-TKEmin < 0)
  {
    cout<<"<W> Estimated TKEmin > Qvalue: TKEmin = "<<TKEmin<<", Qvalue = "<<Qvalue<<endl;
    cout<< "    A_heavy = "<<I_A_heavy_sci<<", A_light; "<<I_A_light_sci<<", Mode = "<<I_Mode<<endl;
    cout<< "    Z_heavy = "<<I_Z_heavy_sci<<", Z_light; "<<I_Z_light_sci<<", Mode = "<<I_Mode<<endl;
    TKEmin = Qvalue - 1;  
  }
  // Limitation of E* to avoid impossible values: Switch from shell-shape to macroscopic shape
  if(Eexc_light_mean > (Qvalue-TKEmin) * I_A_heavy_sci / I_A_sci)
    Eexc_light_mean = Edefo[1][1][I_Z_light_sci];
  if(Eexc_heavy_mean > (Qvalue-TKEmin) * I_A_light_sci / I_A_sci)
    Eexc_heavy_mean = Edefo[1][2][I_Z_heavy_sci];

  // Contribution to width of TKE and TXE by fluctuation of neck length
  //   ESIGDEFOlight = sqr(ESIGDEFOLight^2 + (SIGENECK * R_Z_Light * R_Z_heavy)^2)    
  //   ESIGDEFOheavy = sqr(ESIGDEFOheavy^2 + (SIGENECK * R_Z_light * R_Z_heavy)^2)

  Eexc_light = -1.;
  while((Eexc_light < 0.) | (Eexc_light > (Qvalue-TKEmin) * I_A_heavy_sci / I_A_sci))
    Eexc_light = PGauss(Eexc_light_mean,ESIGDEFOlight);
  Eexc_heavy = -1.;
  while((Eexc_heavy < 0.) | (Eexc_heavy > (Qvalue-TKEmin) * I_A_light_sci / I_A_sci)) 
    Eexc_heavy = PGauss(Eexc_heavy_mean,ESIGDEFOheavy);
  // If Eexc_light_mean > 50 or Eexc_heavy_mean > 50 Then    
  // Print "Mode,Zlight,Alight,Eexc_light_mean",I_Mode,I_Z_light,I_A_light,Eexc_light_mean
  // Print "Mode,Zheavy,Aheavy,Eexc_heavy_mean",I_Mode,I_Z_heavy,I_A_heavy,Eexc_heavy_mean
  // End If

  /* Assumption: width in TKE is the */
  /* quadratic sum of widths in defo and coll. */
  /* Remark of caution: The width of the TKE for fixed mass contains often
     several fission modes. The width in Lang et al. for fixed Z contains several A,
     which contributes already with about 3% to the width. Therefore, the
     width in TXE (or TKE) for fixed A and fixed mode may be much smaller than
     4 or 5 percent! */


  /* Temperatures of fragments */

  if (I_Mode <= 4)
  {
    Tlight = Temp[I_Mode][1][I_A_light_sci];
    Theavy = Temp[I_Mode][2][I_A_heavy_sci];
  }
  if(I_Mode == 5)
  {
    Tlight = Temp[1][2][I_A_light_sci];
    Theavy = Temp[1][2][I_A_heavy_sci];
  }
  if(I_Mode == 6)
  {
    Tlight = Temp[2][2][I_A_light_sci];
    Theavy = Temp[2][2][I_A_heavy_sci];
  }

  /* Intrinsic excitation energies of fragments */

  float Delta_E_Q;
  int I_E_Q;
  float E_intr_tot, Sigma_E_intr;

  E_intr_light_mean = EPART[I_Mode][1][I_A_light_sci];
  E_intr_heavy_mean = EPART[I_Mode][2][I_A_heavy_sci];

  if(I_Mode == 0)
  {  
    // At high E*(CN), where the SL (S0) mode is dominant,
    // TKE is only a macroscopic quantity, microscopic contributions
    // must be suppressed.
    // TKE is derived from empirical formula.

    I_E_Q = 0;
    bool Repeat_E_Q=true;
    while(Repeat_E_Q)
    {
      Repeat_E_Q=false;
      I_E_Q = I_E_Q + 1;
      E_intr_light = -1.;
      E_intr_tot = E_intr_light_mean + E_intr_heavy_mean;
      if(E_intr_tot < 10)
        // Fit function: Fluctuation of energy division from
        // width of phase-space factor with Fermi-gas level density 
        // (below 10 MeV: constant-T level density).
        // Fluctuation of E_intr_tot not changed.
        Sigma_E_intr = E_intr_tot * 0.47 * exp(-sqrt(10/160)) / I_E_Q / 2.35;
      else
        Sigma_E_intr = E_intr_tot * 0.47 * exp(-sqrt(E_intr_tot/160)) / I_E_Q / 2.35;

      E_intr_light = PGauss(E_intr_light_mean, Sigma_E_intr);
      E_intr_heavy = E_intr_tot - E_intr_light; 
      // Print
      // Print "At scission ";I_A_light;I_A_heavy;E_intr_light,E_intr_heavy     

      // Add E* gain after scission due to shell and pairing of fragments.
      // This way, they are removed from TKE (further down in the code).
      // This is reasonable, because fragment shell and pairing energies appear after scisssion.
      // Thus, they do not appear in the TKE.
      E_intr_light = E_intr_light - 
        AME2012(I_Z_light_sci,I_A_light_sci) + LDMass(I_Z_light_sci,I_A_light_sci,0.)   // shell and pairing
        - 2.0 * 12.0 / sqrt(I_A_light_sci);     // general shift    
      E_intr_heavy = E_intr_heavy - 
        AME2012(I_Z_heavy_sci,I_A_heavy_sci) + LDMass(I_Z_heavy_sci,I_A_heavy_sci,0.)  // shell and pairing
        - 2.0 * 12.0 / sqrt(I_A_heavy_sci);     // general shift    
      float Epre_mean,TKE_mac,TKE1,TKE2,dmod;
      Epre_mean = 0.4 * De_Saddle_Scission(pow(I_Z_CN,2)/pow(I_A_CN,0.333333),ESHIFTSASCI_coll);  
      // Mean pre-scission kinetic energy 
      dmod = (DNECK + 1); //  adjusted to TKE of 250Fm(E*=45 MeV) (in conflict with Loveland 1970?)
      TKE_mac = 1.44*I_Z_light_sci*I_Z_heavy_sci/  
        (1.16*(pow(I_A_light_sci,0.33333)*(1.0+0.66667*Beta[1][1][I_Z_light_sci]) +
               pow(I_A_heavy_sci,0.33333)*(1.0+0.66667*Beta[1][2][I_Z_heavy_sci])) +  dmod) 
        - 1.44*pow(I_Z_sci/2,2) /    
        (1.16*( pow(I_A_sci/2,0.33333)*(1.0+0.66667*Beta[1][1][I_Z_sci/2]) +
                pow(I_A_sci/2,0.33333)*(1.0+0.66667*Beta[1][2][I_Z_sci/2])) + DNECK);
      // (Variation of TKE from macroscopic formula of Wilkins et al.
      //  relative to value at symmetry.)
      Delta_E_Q = - AME2012(I_Z_sci,I_A_sci) + LDMass(I_Z_sci,I_A_sci,0.) 
        + LyMass(I_Z_light_sci,I_A_light_sci,0.) + LyMass(I_Z_heavy_sci,I_A_heavy_sci,0.) 
        - LyMass(I_Z_sci,I_A_sci,0) 
        - (LyMass(I_Z_sci/2,I_A_sci/2,0.) + LyMass(I_Z_sci/2,I_A_sci/2,0.) - LyMass(I_Z_sci,I_A_sci,0)) 
        + TKE_mac   // replace macr. variation of Q value by Wilkins formula              
        + Epre_mean;  // add pre-scission TKE  
      // Print "Delta_E_Q",Delta_E_Q,I_Z_sci,I_A_sci,TKE_mac,Epre_mean                
      // positives Delta_E_Q (ohne TKE_mac) erhöht TKE in den Flanken (mehr als TKE_mac)   
      // positives Epre_mean erhöht TKE
      // positives TKE_mac vermindert TKE in den Flanken   

      E_intr_light = E_intr_light - I_A_light_sci/I_A_sci * (Delta_E_Q);
      E_intr_heavy = E_intr_heavy - I_A_heavy_sci/I_A_sci * (Delta_E_Q); 
      E_intr_light_S0_mac = Max(0.0,E_intr_light);
      E_intr_heavy_S0_mac = Max(0.0,E_intr_heavy);

      // Print, "After fission: ";E_intr_light,E_intr_heavy

      if((E_intr_light_S0_mac < 0) | (E_intr_heavy_S0_mac < 0))
      {
        if(I_E_Q < 3)
        {
          cout<< "Reapeat_E_Q"<<endl;
          Repeat_E_Q=true;
        } 
        cout<<"6326: E_intr_light/heavy = "<<E_intr_light<<"/"<<E_intr_heavy<<endl;
        cout<<"6327: Mean values: "<<E_intr_light_mean<<" "<<E_intr_heavy_mean<<endl;
        E_intr_light_S0_mac = Max(0.0,E_intr_light);
        E_intr_heavy_S0_mac = Max(0.0,E_intr_heavy);
      }   
    }
    //-------------------------------------------------------- 
    float TXE_shift,EtotS0;
    TXE_shift=EtotS0=0;
    if(Mod(I_Z_CN,2) == 0)
      // Even-odd fluctuation of TXE acc. to Lang et al. NPA 345 (1980) 34
      // Lower TXE for completely paired configuration at scission.
      // Only the even-odd effect at symmetry is considered, because the
      // energy gain by the asymmetry-driven even-odd effect goes to E_intr
      // of the heavy fragment due to energy sorting.
      if(Mod(I_Z_light_sci,2) == 0)
        if(rn->Rndm() < PEOZ[I_Mode][1][int(0.5*I_A_sci)])
          TXE_shift = -12./sqrt(I_A_sci);

    E_intr_light_S0_mic = E_intr_light_mean + TXE_shift * E_intr_light_mean / (E_intr_light_mean + E_intr_heavy_mean);
    if(E_intr_light_S0_mic < 0.)
      E_intr_light_S0_mic = 0.;
    E_intr_heavy_S0_mic = E_intr_heavy_mean + TXE_shift * E_intr_heavy_mean / (E_intr_heavy_mean + E_intr_light_mean);             
    if(E_intr_heavy_S0_mic < 0.)
      E_intr_heavy_S0_mic = 0.;      

    E_intr_light_S0_mic = E_intr_light_S0_mic - LyPair(I_Z_heavy_sci,I_A_heavy_sci);  
    E_intr_heavy_S0_mic = E_intr_heavy_S0_mic - LyPair(I_Z_heavy_sci,I_A_heavy_sci);  

    //-------------------------------------------------------- 

    EtotS0 = Max(E_Exc_S0 + E_diss_Scission,0.0);
    RW_mac = Gaussintegral(EtotS0-20.0,5.0);
    E_intr_light_mean = (1.0 - RW_mac) * E_intr_light_S0_mic + RW_mac * E_intr_light_S0_mac;
    E_intr_heavy_mean = (1.0 - RW_mac) * E_intr_heavy_S0_mic + RW_mac * E_intr_heavy_S0_mac;

    //-------------------------------------------------------- 

    E_intr_light = -1.;
    while(E_intr_light < 0.)
      E_intr_light = PGauss(E_intr_light_mean,EexcSIGrel* sqrt(2) * 0.5*sqrt(pow(E_diss_Scission,2)+pow(SIGENECK,2))+0.5);	  
    E_intr_heavy = -1.;
    while (E_intr_heavy < 0.)
      E_intr_heavy = PGauss(E_intr_heavy_mean,EexcSIGrel* sqrt(2) * 0.5*sqrt(pow(E_diss_Scission,2)+pow(SIGENECK,2))+0.5);

    //-------------------------------------------------------- 

    if((E_intr_light < 0)|(E_intr_heavy < 0)) 
      cout<< "E_intr < 0"<<E_intr_light<<" "<<E_intr_heavy<<endl;        
  }
  else // If I_Mode <> 0
  {
    float TXE_shift=0;
    if(Mod(I_Z_CN,2) == 0)  
      // Even-odd fluctuation of TXE acc. to Lang et al. NPA 345 (1980) 34
      // Lower TXE for completely paired configuration at scission.
      // Only the even-odd effect at symmetry is considered, because the
      // energy gain by the asymmetry-driven even-odd effect goes to E_intr
      // of the heavy fragment due to energy sorting.
      if(Mod(I_Z_light_sci,2) == 0)
        if(rn->Rndm() < PEOZ[I_Mode][1][int(0.5*I_A_sci)])
          TXE_shift = -12./sqrt(I_A_sci);

    E_intr_light_mean = E_intr_light_mean + TXE_shift * E_intr_light_mean / (E_intr_light_mean + E_intr_heavy_mean);
    if(E_intr_light_mean < 0.)
      E_intr_light_mean = 0.;
    E_intr_heavy_mean = E_intr_heavy_mean + TXE_shift * E_intr_heavy_mean / (E_intr_heavy_mean + E_intr_light_mean);            
    if(E_intr_heavy_mean < 0.)
      E_intr_heavy_mean = 0.;

    E_intr_light = -1.;
    while(E_intr_light < 0.)
    {
      E_intr_light = PGauss(E_intr_light_mean,EexcSIGrel* sqrt(2) * 0.5*sqrt(pow(E_diss_Scission,2)+pow(SIGENECK,2))+0.5);
      E_intr_light = E_intr_light - LyPair(I_Z_light_sci,I_A_light_sci);
    }
    E_intr_heavy = -1.;
    while(E_intr_heavy < 0.)
    {
      E_intr_heavy = PGauss(E_intr_heavy_mean,EexcSIGrel* sqrt(2) * 0.5*sqrt(pow(E_diss_Scission,2)+pow(SIGENECK,2))+0.5);
      E_intr_heavy = E_intr_heavy - LyPair(I_Z_heavy_sci,I_A_heavy_sci);
    }
    //Print "6331: ";I_mode,E_intr_light_mean,E_intr_light
    //Print "6332: ";I_mode,E_intr_heavy_mean,E_intr_heavy
    //Print "6333: ";EexcSIGrel,E_diss_Scission,SIGENECK

    //    E_intr_light = E_intr_light - Lypair(I_Z_light,I_A_light)
    //    E_intr_heavy = E_intr_heavy - Lypair(I_Z_heavy,I_A_heavy)  
    /* Staggering of BE of fragments by pairing */   
    /* Assumption: pairing only felt in the lowest nuclear levels, */
    /* at the end of the evaporation cascade */
    /* (This should be a good assumption for higher excitation energies. */
    /* Some deviations occur due to the even-odd effect at low exc. energies. */

  } // If I_Mode <> 0

  Qvalue = AME2012(I_Z_sci,I_A_sci) - (AME2012(I_Z_heavy_sci,I_A_heavy_sci) + AME2012(I_Z_light_sci,I_A_light_sci));

  Eexc_light = Eexc_light + E_intr_light;
  Eexc_heavy = Eexc_heavy + E_intr_heavy;
  /* Now: deformation + intrinsic excitation energy */
  /* Collective energy of fragments */

  Ecoll_mean = ECOLLFRAC * (De_Saddle_Scission(pow(I_Z_CN,2) / pow(I_A_CN,0.33333),ESHIFTSASCI_coll) - E_tunn);

  if(Ecoll_mean < 0.)
    Ecoll_mean = 0.;

  /* Experimental data of prompt neutron yields show an enhancement of the */  
  /* neutron yield for odd-Z CN, corresponding to an enhanced E* by about 1.6 MeV */
  /* The enhancement is equally shared to the light and the heavy fragment. */
  /* The neutron number of the CN has no influence. */ 
  /* The origin of this effect is not clear. */
  /* By technical reasons, this additional energy is introduced here into the */
  /* collective energy at scission, because this energy is divided equally between both fragments. */
  /*  KHS, 31. Jan. 2012 */
  //  If  (I_Z_CN Mod 2) = 1 Then Ecoll_mean = Ecoll_mean + 1.6
  // This is not used any more. It seems to be wrong or to be replaced by something else.

  //  Assuming positive correlation of Ecoll_heavy and Ecoll_light.
  //  This is probably not realistic!
  //  See Nix & Swiatecki, Nucl. Phys. 71 (1965) 1    
  /*  Ecoll = -1.
      While Ecoll < 0.
      Ecoll = PGauss(Ecoll_mean,EexcSIGrel* Ecoll_mean)
      Wend  
      Ecoll = Ecoll + E_coll_saddle(I_Mode)
      Ecoll_heavy = 0.5E0 * Ecoll
      Ecoll_light = 0.5E0 * Ecoll */

  //   Assuming no correlation (partly positive, partly negative correlation)
  //   See Nix & Swiatecki, Nucl. Phys. 71 (1965) 1    

  float Rmean,Rsigma; 
  Rmean = 0.5 * Ecoll_mean;
  Rsigma = sqrt(0.5) * EexcSIGrel * Ecoll_mean;
  Ecoll_light = -1;    
  while(Ecoll_light < 0.) 
    Ecoll_light = PGauss(Rmean,Rsigma);
  Ecoll_light = Ecoll_light + 0.5 * E_coll_saddle[I_Mode];
  Ecoll_heavy = -1.;
  while(Ecoll_heavy < 0.) 
    Ecoll_heavy = PGauss(Rmean,Rsigma);
  Ecoll_heavy = Ecoll_heavy + 0.5E0 * E_coll_saddle[I_Mode];
  Ecoll = Ecoll_heavy + Ecoll_light;

  Eexc_light = Eexc_light + Ecoll_light;
  Eexc_heavy = Eexc_heavy + Ecoll_heavy;
  /* Now: also collective excitation energy added */


  /* Excitation energy not including the rotational energy at scission: */

  /*** Angular momentum of fragments ***/           

  float J_Frag_light,J_Frag_heavy;
  float J_Frag_rot_light,J_Frag_rot_heavy; // collective spin
  float N_J_attempt;
  N_J_attempt = 0;

  bool J_attempt=true;
  while(J_attempt)
  {
    J_attempt=false;
    N_J_attempt = N_J_attempt + 1;
    if(I_Mode <= 4)
    {
      J_Frag_light = PLinGauss(SpinRMSNZ[I_Mode][1][I_N_light_sci][I_Z_light_sci]/sqrt(2.0)) - 0.5;
      J_Frag_heavy = PLinGauss(SpinRMSNZ[I_Mode][2][I_N_heavy_sci][I_Z_heavy_sci]/sqrt(2.0)) - 0.5;
    }
    if(I_Mode == 5)
    {
      J_Frag_light = PLinGauss(SpinRMSNZ[1][2][I_N_light_sci][I_Z_light_sci]/sqrt(2.0)) - 0.5;
      J_Frag_heavy = PLinGauss(SpinRMSNZ[1][2][I_N_heavy_sci][I_Z_heavy_sci]/sqrt(2.0)) - 0.5;
    }
    if(I_Mode == 6)
    {
      J_Frag_light = PLinGauss(SpinRMSNZ[2][2][I_N_light_sci][I_Z_light_sci]/sqrt(2.0)) - 0.5;
      J_Frag_heavy = PLinGauss(SpinRMSNZ[2][2][I_N_heavy_sci][I_Z_heavy_sci]/sqrt(2.0)) - 0.5;
    }

    if(J_Frag_light < 0)
      J_Frag_light = 0;

    if(J_Frag_heavy < 0)
      J_Frag_heavy = 0;

    IfragEff_light = U_Ired(I_Z_light_sci,I_A_light_sci);  
    Erotlight =  J_Frag_light*(J_Frag_light+1)/(2*IfragEff_light);
    // rotational energy at scission

    IfragEff_heavy = U_Ired(I_Z_heavy_sci,I_A_heavy_sci);
    Erotheavy =  J_Frag_heavy*(J_Frag_heavy+1)/(2*IfragEff_heavy);
    // rotational energy at scission


    float TXE, E_total;
    /* Kinetic energies of fragments */

    // TXE includes E*_CN 

    TXE = Eexc_light + Eexc_heavy + Erotlight + Erotheavy; // provisional value   

    E_total = Qvalue + R_E_exc_GS;
    TKE = E_total - TXE;

    static int Ntimes = 0;     
    if(TKE < 0)
    {
      //    If TKE < 0 Or I_Mode =0 Then
      if(N_J_attempt <= 3)
        J_attempt=true;
      float TXEcorr;
      cout<< "<E> Event with excessive excitation energy found."<<endl;;
      cout<< "I_A_light,I_Z_light,I_A_heavy,I_Z_heavy: "<<I_A_light_sci<<" "<<I_Z_light_sci<<" "<<I_A_heavy_sci<<" "<<I_Z_heavy_sci<<endl;
      cout<< "Q value,TKE; "<<Qvalue<<" "<<TKE<<endl;
      cout<< "Eexc_light,Eexc_heavy,Erot_light,Erot_heavy: "<<Eexc_light<<" "<<Eexc_heavy<<" "<<Erotlight<<" "<<Erotheavy<<endl;
      cout<< "J_rms_light,J_rms_heavy: "<<SpinRMSNZ[I_Mode][1][I_N_light_sci][I_Z_light_sci]<<" "<<SpinRMSNZ[I_Mode][2][I_N_heavy_sci][I_Z_heavy_sci]<<endl;
      cout<< "J_light,J_heavy "<<J_Frag_light<<" "<<J_Frag_heavy<<endl;
      cout<< "Ecoll_light,Eintr_light: "<<Ecoll_light<<" "<<E_intr_light<<endl;
      cout<< "Ecoll_heavy,Eintr_heavy: "<<Ecoll_heavy<<" "<<E_intr_heavy<<endl;
      if(I_Mode == 0)
      {
        cout<<"I_Mode,Edefo(light),Edefo_mic(heavy),Edefo_mac(heavy): "<<I_Mode<<" "<< Edefo[I_Mode+1][1][I_Z_light_sci]<<" "<<Edefo[0][2][I_Z_heavy_sci]<<" "<<Edefo[I_Mode+1][2][I_Z_heavy_sci]<<endl;
        cout<<"E_intr_light_S0_mic,..mac,..heavy_S0_mic,..heavy_S0_mac: "<<E_intr_light_S0_mic<<" "<<E_intr_light_S0_mac<<" "<<E_intr_heavy_S0_mic<<" "<<E_intr_heavy_S0_mac<<endl;        
      }
      else
        cout<< "I_Mode,Edefo(light),Edefo_(heavy): "<<I_Mode<<" "<<Edefo[I_Mode+1][1][I_Z_light_sci]<<" "<<Edefo[I_Mode+1][2][I_Z_heavy_sci]<<endl;
      TXEcorr = TXE + TKE - 1;
      TKE = 1;
      Eexc_light = Eexc_light * TXEcorr/TXE;
      Eexc_heavy = Eexc_heavy * TXEcorr/TXE;
      TXE = TXEcorr;
      Ntimes = Ntimes + 1;
      if(Ntimes > 99)
      {
        cout<< "<S> Calculation stopped due to numerical problem."<<endl;
        exit(-1);
      }
    }

  }
  TXElight = Eexc_light  + Erotlight;    // final values
  TXEheavy = Eexc_heavy  + Erotheavy;    //  "      " 

  // Print Erotlight,Erotlight_FF,Erotheavy,Erotheavy_FF      


  TXE = TXElight + TXEheavy;            //  "      "  


  Ekinlight_sci = TKE * I_A_heavy_sci / I_A_sci;
  Ekinheavy_sci = TKE * I_A_light_sci / I_A_sci;
  // Print 7058,TKE,I_A_light_sci,I_A_heavy_sci      
  // Print "7057 ";I_A_light_sci,Erotlight+Erotheavy,TXE,TKE,Qvalue, _    
  //               AME2012(I_Z_sci,I_A_sci) - _
  //              (AME2012(I_Z_heavy_sci,I_A_heavy_sci)  _
  //             + AME2012(I_Z_light_sci,I_A_light_sci))      

  // cout<< Ekinlight_sci<<" "<<Ekinheavy_sci<<endl;
}

void GEF::FromCMtoLAB(void)
{
  Thcm_light = acos(rn->Rndm()*2.-1.);
  Thcm_heavy = pi-Thcm_light;
  Phcm_light = rn->Rndm()*2*pi-pi;
  if(Phcm_light == 0)
    Phcm_heavy = pi;
  else
    Phcm_heavy = Phcm_light - Phcm_light/abs(Phcm_light)*pi;

  float Gfis,Bfis;
  Gfis = 1.;
  if(Afis>0){
    //Gfis = FIS.GetEnergy()/(FIS.Mass()+Ex)+1.;
    Gfis = Efis/(FIS.Mass()+Ex)+1.;
  }
  Bfis = sqrt(1.-1./pow(Gfis,2));

  float M_light,M_heavy;
  M_light=M_heavy=0;
  M_light = AME2012(I_Z_light_sci,I_A_light_sci)+I_Z_light_sci*(931.494+7.28897061)+(I_A_light_sci-I_Z_light_sci)*(931.494+8.0713171);
  M_heavy = AME2012(I_Z_heavy_sci,I_A_heavy_sci)+I_Z_heavy_sci*(931.494+7.28897061)+(I_A_heavy_sci-I_Z_heavy_sci)*(931.494+8.0713171);
  float Gcm_light, Gcm_heavy;
  float Bcm_light, Bcm_heavy;
  Gcm_light=Gcm_heavy=1.;
  if(M_light>0)
    Gcm_light = Ekinlight_sci/M_light+1.;
  if(M_heavy>0)
    Gcm_heavy = Ekinheavy_sci/M_heavy+1.;
  Bcm_light = sqrt(1.-1./pow(Gcm_light,2));
  Bcm_heavy = sqrt(1.-1./pow(Gcm_heavy,2));

 //cout << Ekinlight_sci << " " << Ekinheavy_sci << " " << Ekinlight_sci+ Ekinheavy_sci << endl;

  KElab_light = KElab_heavy= 0;
  KElab_light = (Gfis*Gcm_light*M_light +Gfis*Bfis*Gcm_light*Bcm_light*M_light*cos(Thcm_light))-M_light;
  KElab_heavy = (Gfis*Gcm_heavy*M_heavy +Gfis*Bfis*Gcm_heavy*Bcm_heavy*M_heavy*cos(Thcm_heavy))-M_heavy;

  Glab_light=Glab_heavy=1.;
  if(M_light>0)
    Glab_light = KElab_light/M_light+1.;
  if(M_heavy>0)
    Glab_heavy = KElab_heavy/M_heavy+1.;
  Blab_light = sqrt(1.-1./pow(Glab_light,2));
  Blab_heavy = sqrt(1.-1./pow(Glab_heavy,2));
  vlab_light = Blab_light*29.9792458;
  vlab_heavy=Blab_heavy*29.9792458;
 
  //Note that I am not taking account the angle of the fissioning system
  Thlab_light = Thlab_heavy=0;
  if(Blab_light>0)
    Thlab_light = acos(1./Glab_light/Blab_light*(Gfis*Bfis*Gcm_light+Gfis*Gcm_light*Bcm_light*cos(Thcm_light)));
  if(Blab_heavy>0)
    Thlab_heavy = acos(1./Glab_heavy/Blab_heavy*(Gfis*Bfis*Gcm_heavy+Gfis*Gcm_heavy*Bcm_heavy*cos(Thcm_heavy)));
  Phlab_light = Phcm_light;
  Phlab_heavy = Phcm_heavy;

  //Thlablab_light = acos(-1.*sin(Thlab_light)*cos(Phlab_light)*sin(FIS->GetThfis())-sin(Thlab_light)*sin(Phlab_light)*sin(FIS->GetThfis())+cos(Thlab_light)*cos(FIS->GetThfis()));
  //cout<<Thlab_light*180./pi<<" "<<Thlablab_light*180./pi<<endl;
}

///////////////////////////////////CHARGE STATES = f(Z,vlab) from E. Baron NIM A 328,177 (1993)//////////////////
float Qmean_Baron_lowZ(float Z, float v)
{
  float par[2]={-83.275,0.447};
  float beta = v/29.9792458;
  return Z*(1.-1.*exp(par[0]*beta/pow(Z,par[1])));
}
float Qmean_Baron(float Z, float v)
{
  float par[3]={-12.905,0.2124,-0.00122};
  float qmean1 = Qmean_Baron_lowZ(Z,v);
  float qmean2 = qmean1*(1.-exp(par[0]+par[1]*Z+par[2]*Z*Z));
  if(Z<54)
    return qmean1;
  else
    return qmean2;
}

float Qsigma_Baron(float Z, float v)
{
  float par[5]={0.5,1.67,0.07535,0.19,-0.2654};
  float qmean1 = Qmean_Baron_lowZ(Z,v);
  float sigma1 = par[0]*sqrt(qmean1*(1.-pow(qmean1/Z,par[1])));
  float sigma2 = sqrt(qmean1*(par[2]+par[3]*qmean1/Z+par[4]*pow(qmean1/Z,2)));
  if(Z<54)
    return sigma1;
  else
    return sigma2;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GEF::FissionProbability(void)
{
  // This section is for calculating the first-chance fission probablity for a large
  // range of excitation energy with an analytical function (4<Ex<80)MeV.
  // In contrast to the Monte-Carlo method, used for the "normal" calculation,
  // this works also for very low fission probabilities.
  float E_left;
  int Z_left,A_left;
  float SNTF,SNTFnopair,SNmean;
  float SP,BP,SPTF,SPTFnopair,SPmean,BPmean,BPnopair;
  float BF,BFNP,BFMNP,BFA,BFANP,BFB,BFBNP,BFM;
  float af_an;
  float Ztrans = 86.0;
  float dZtrans = 2;
  float DUN,DU;
  double Rho_cn;
  float G,GN,GP,GF;
  float T_n,T_f,T_p,Tm;
  float FtunnA,FtunnB;
  float hbom = 1.;  // effective curvature of fission barrier
  float Tequi = hbom /(2 * pi);  // equivalent temperature parameter
  float F1_coll = 10;
  float F1_diss = 20; 
  float GA = (0.14 / sqrt(pi/2));  // triaxiality (compared to quadrupole deformation)
  float GB = 0.5;                 // mass asymmetry (compared to quadrupole deformation)
  float GAmod,GBmod;
  float DUF = 0;   // Shell effect at barrier
  float ypsilon,l_lim_2;
  float Frot;
  float Red_diss;
  float Ggamma;
  if(Ex<4)
    PF=0;
  else
  {
    E_left = Ex;
    Z_left = P_Z_CN;
    A_left = P_A_CN;
    SN = AME2012(Z_left,A_left-1) - AME2012(Z_left,A_left);             
    SNTF = U_MASS(Z_left,A_left-1) + LyPair(Z_left,A_left-1) - (U_MASS(Z_left,A_left) + LyPair(Z_left,A_left));
    SNTFnopair = U_MASS(Z_left,A_left-1) - U_MASS(Z_left,A_left); 
    SNmean = 0.5 * (U_MASS(Z_left,A_left-2) - U_MASS(Z_left,A_left));

    SP = AME2012(Z_left-1,A_left-1) - AME2012(Z_left,A_left);             
    BP = SP + 1.44 * (Z_left-1) / (2.0 * (1 + pow(A_left-1,0.333333)));  
    SPTF =  U_MASS(Z_left-1,A_left-1) + LyPair(Z_left-1,A_left-1) - (U_MASS(Z_left,A_left) + LyPair(Z_left,A_left));
    SPTFnopair = U_MASS(Z_left-1,A_left-1) - U_MASS(Z_left,A_left);
    BPnopair = SPTFnopair + 1.44 * (Z_left-1) / (2.0 * (1 + pow(A_left-1,0.333333)));  
    SPmean = 0.5 * (U_MASS(Z_left-2,A_left-2) - U_MASS(Z_left,A_left));      
    BPmean = SPmean + 1.44 * (Z_left-1) / (2.0 * (1 + pow(A_left-1,0.333333)));                  

    BF = BFTF(P_Z_CN,A_left,1);     // with shell and enhanced pairing at barrier
    BFNP = BFTF(P_Z_CN,A_left,3);   // no enhanced pairing at barrier
    //     BFLD = BFTF(P_Z_CN,A_left,0)   // liquid drop
    //     BFA = BFTFA(P_Z_CN,A_left,1)   // with shell and enhanced pairing at barrier
    BFA = BFTFA(P_Z_CN,A_left,1) + 10 / sqrt(A_left) * Mod(A_left-P_Z_CN + 1,2);    
    // e-o in levdens at Bf
    BFANP = BFTFA(P_Z_CN,A_left,3); // no enhanced pairing at barrier
    BFB = BFTFB(P_Z_CN,A_left,1) + 10 / sqrt(A_left) * Mod(A_left-P_Z_CN + 1,2);
    // e-o-in levdens at Bf
    BFBNP = BFTFB(P_Z_CN,A_left,3); // no enhanced pairing at barrier
    BFM = Max(BFA,BFB);             // highest saddle
    BFMNP = Max(BFANP,BFBNP);       // highest saddle without enhanced pairing
    af_an = 1.2 - 0.2 * Gaussintegral(P_Z_CN - Ztrans, dZtrans);
    DUN = U_SHELL(P_Z_CN,A_left-1);
    DU = U_SHELL(P_Z_CN,A_left);
    Rho_cn = U_levdens(P_Z_CN,A_left, E_left, 1, 1, Tscale, Econd,1);
    //   Smean = 0.5 * (SNmean + BFMNP)
    if(E_left > SN)
    {
      T_n = U_Temp(P_Z_CN,A_left-1,E_left-SN,1,1,Tscale,Econd);
      GN = pow(A_left-1,0.66667) * 0.13 * pow(T_n,2) * U_levdens(P_Z_CN,A_left-1,E_left-SN,1,1,Tscale,Econd,1) / Rho_cn; 
    }
    else
      GN = 0; 
    if(E_left > BP)
    {
      T_p = U_Temp(Z_left-1,A_left-1,E_left-BP,1,1,Tscale,Econd);
      GP = pow(A_left-1,0.666667) * 0.13 * pow(T_p,2) * U_levdens(Z_left-1,A_left-1,E_left-BP,1,1,Tscale,Econd,1) / Rho_cn; 
    }
    else
      GP = 0;     

    float FcollA, EcollA, FcollB, EcollB;
    EcollA = E_left - BFANP;
    if(EcollA < 0)
      EcollA = 0;
    FcollA = GA / (1. - (1. - GA) * exp(-EcollA/F1_coll)); 
    FtunnA = 1./(1.+exp(-(E_left-BFA)/Tequi));   
    GAmod = GA / (FtunnA * FcollA);
    EcollB = E_left - BFBNP;
    if(EcollB < 0)
      EcollB = 0;
    FcollB = GB / (1. - (1. - GB) * exp(-EcollB/F1_coll)); 
    FtunnB = 1./(1.+exp(-(E_left-BFB)/Tequi));  
    GBmod = GB / (FtunnB * FcollB);
    T_f = U_Temp2(P_Z_CN,A_left,E_left-BFMNP,DUF,0,Tscale,Econd);
    G = GAmod * exp((BFANP-BFMNP)/T_f) + GBmod * exp((BFBNP-BFMNP)/T_f);
    ypsilon = 1.0 - pow(P_Z_CN,2) / (P_A_CN * 50.0);  // 1 - fissility
    ypsilon = Max(ypsilon,0.1);   // limitation for very heavy nuclei
    l_lim_2 = pow(15.0,2) / (ypsilon/0.28) * T_f / TEgidy(P_A_CN,DUF,1.0);  // Hasse & Myers
    Frot = exp( pow(P_I_rms_CN,2)/l_lim_2 );  // enhancement of GF by angular momentum
    GF = T_f / G * U_levdens(P_Z_CN,A_left-1,E_left-BFMNP,0,0,Tscale,Econd,af_an) / Rho_cn * Frot;
    Red_diss = exp(-(E_left - BFMNP)/(F1_diss)); 
    GF = GF * Red_diss; 
    // Low level density below the pairing gap          
    /*    If P_Z_CN Mod 2 = 0 And A_left Mod 2 = 0 Then
          If E_left - BFM < 2 * 12 / sqr(A_left) Then
          GF = GF * Max(0.01 + 0.99*(Abs( (E_left-BFM) / (2*12/sqr(A_left) )) ),0.1  ) 
          End If
          Else
          If P_Z_CN Mod 2 = 0 Or (A_left - P_Z_CN) Mod 2 = 0 Then
          If E_left - BFM < 12 / sqr(A_left) Then
          GF = GF * Max(0.01 + 0.99*(Abs( (E_left-BFM) / (12/sqr(A_left) )) ),0.1  ) 
          End If
          End If 
          End If  */

    Tm = U_Temp(P_Z_CN,A_left,E_left,1,0,Tscale,Econd);       // emitting nucleus 
    Ggamma = 0.624 * pow(A_left,1.6) * pow(Tm,5);    // in meV (Ignatyuk, Bologna)
    Ggamma = Ggamma * 1.e-9;        // in MeV

    PF = GF / (GF + GN + GP + Ggamma);  // Pfis does not include multi-chance fission  
  }
  // cout<<"Ex: "<<Ex<<" PF: "<<PF<<endl;
}


void GEF::InitCompound(double vEx, double vEfis, double vLfis, double vThfis, double vPhfis)
{
  Ex = vEx;
  Efis = vEfis;
  Lfis = vLfis;
  Thfis = vThfis;
  Phfis = vPhfis;
}

void GEF::Treat(void)
{
  /*Ex = FIS.GetExcitationEnergy();
  Lfis = 0;
  Efis = FIS.GetEnergy();
  Thfis = 0;
  Phfis = 0;*/

  P_E_exc = Ex;
  P_I_rms_CN = Lfis;
  Spin_pre_fission = P_I_rms_CN; //???????
  PotCurv_FM();
  R_E_exc_used = P_E_exc;
  EnergyTrans();
  BarriersEx_FM();
  TCollective();
  MeanValues_FM();
  EnergyDependence();
  RelationsZ_A_FM();//second interaccion after energy dependence
  Yields_FM();
  MassWidths_FM();
  ShellEff_FM();

  IntrinsicT();
  IntrinsicEx();
#ifndef GEF_SIMPLE
  FF_angularMomentum();
#endif

  FissionModeChoise();
  ZAChoise();
  N_Saddle_Scission();
  FragmentsEnergy();
  FissionProbability();

  FromCMtoLAB();

  Qffl = Qffh = 0;
  Qffl = int(rn->Gaus(Qmean_Baron(1.*I_Z_light_sci, vlab_light),Qsigma_Baron(1.*I_Z_light_sci, vlab_light)) +.5);
  Qffh = int(rn->Gaus(Qmean_Baron(1.*I_Z_heavy_sci, vlab_heavy),Qsigma_Baron(1.*I_Z_heavy_sci, vlab_heavy))+0.5);
  Brhoffl = Brhoffh=0.;
  if(Qffl>0)
    Brhoffl = 3.105*I_A_light_sci/Qffl*Blab_light*Glab_light;
  if(Qffh>0)
    Brhoffh = 3.105*I_A_heavy_sci/Qffh*Blab_heavy*Glab_heavy;
}

/*void GEF::outAttach(TTree *outT)
  {
  outT->Branch("Zffl",&I_Z_light_sci,"Zffl/I");
  outT->Branch("Zffh",&I_Z_heavy_sci,"Zffh/I");
  outT->Branch("Affl",&I_A_light_sci,"Affl/I");
  outT->Branch("Affh",&I_A_heavy_sci,"Affh/I");
  outT->Branch("Qffl",&Qffl,"Qffl/I");
  outT->Branch("Qffh",&Qffh,"Qffh/I");
  outT->Branch("KECMffl",&Ekinlight_sci,"KECMffl/F");
  outT->Branch("KECMffh",&Ekinheavy_sci,"KECMffh/F");
  outT->Branch("ThCMffl",&Thcm_light,"ThCMffl/F");
  outT->Branch("ThCMffh",&Thcm_heavy,"ThCMffh/F");
  outT->Branch("PhCMffl",&Phcm_light,"PhCMffl/F");
  outT->Branch("PhCMffh",&Phcm_heavy,"PhCMffh/F");
  outT->Branch("KELabffl",&KElab_light,"KELabffl/F");
  outT->Branch("KELabffh",&KElab_heavy,"KELabffh/F");
  outT->Branch("VLabffl",&vlab_light,"VLabffl/F");
  outT->Branch("VLabffh",&vlab_heavy,"VLabffh/F");
  outT->Branch("ThLabffl",&Thlab_light,"Thlabffl/F");
  outT->Branch("ThLabffh",&Thlab_heavy,"Thlabffh/F");
  outT->Branch("PhLabffl",&Phlab_light,"Phlabffl/F");
  outT->Branch("PhLabffh",&Phlab_heavy,"Phlabffh/F");
  outT->Branch("Brhoffl",&Brhoffl,"Brhoffl/F");
  outT->Branch("Brhoffh",&Brhoffh,"Brhoffh/F");
//  outT->Branch("ThLabLabffl",&Thlablab_light,"ThLabLabffl/F");
outT->Branch("FissionProb",&PF,"FissionProb/F");


}*/






