/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author : Pierre Morfouace  contact: pierre.morfouace2@cea.fr     *
 *                                                                           *
 * Creation Date   : Octobre 2020                                            *
 * Last update     :                                                         *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *   This Class hold data for fission decay of a given nuclei                *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "NPFissionDecay.h"
#include <iostream>
#include "NPOptionManager.h"
#include "NPFunction.h"
#include "NPCore.h"
#include "TF1.h"
#include "TRandom.h"
////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Fission Decay //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

NPL::FissionDecay::FissionDecay(std::string compound, std::string fission_model){

  m_FissionModelName = fission_model;
  m_CompoundName = compound;
}

////////////////////////////////////////////////////////////////////////////////
void NPL::FissionDecay::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("FissionDecay");

  if(blocks.size()>0 && NPOptionManager::getInstance()->GetVerboseLevel()){
    cout << endl << "\033[1;35m//// Fission decay found" << std::endl; 
  }
  if(blocks.size()==0) HasFissionToken=0;
  else HasFissionToken=1;

  std::vector<std::string> token = 
  {"CompoundNucleus","FissionModel","VamosChargeStates","Shoot_FF","Shoot_neutron","Shoot_gamma"};

  unsigned int size = blocks.size();
  for(unsigned int i = 0 ; i < size ; i++){
    if(blocks[i]->HasTokenList(token)){
      m_CompoundName = blocks[i]->GetString("CompoundNucleus");
      m_FissionModelName = blocks[i]->GetString("FissionModel");
      m_VamosChargeStates = blocks[i]->GetInt("VamosChargeStates");
      m_shoot_FF = blocks[i]->GetInt("Shoot_FF");
      m_shoot_neutron = blocks[i]->GetInt("Shoot_neutron");
      m_shoot_gamma = blocks[i]->GetInt("Shoot_gamma");

      m_Compound = NPL::Particle(m_CompoundName);
      if(m_FissionModelName=="GEF") m_FissionModel = new GEF(m_Compound);
    }
    else{
      cout << "ERROR: check your input file formatting \033[0m" << endl;
      exit(1);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
bool NPL::FissionDecay::GenerateEvent(string CompoundName, double MEx,double MEK,double MPx,double MPy,double MPz, 
    std::vector<NPL::Particle>& FissionFragments, std::vector<double>& Ex,
    std::vector<double>& DEK,
    std::vector<double>& DPx,std::vector<double>& DPy,std::vector<double>& DPz,
    double& TKE, double &KE1, double& KE2){

  bool worked=false;

  FissionFragments.clear();
  Ex.clear();
  DEK.clear();
  DPx.clear();
  DPy.clear();
  DPz.clear();

  TVector3 Momentum(MPx,MPy,MPz);
  Momentum.Unit();
  double Theta = Momentum.Theta();
  double Phi = Momentum.Phi();
  double Lfis = 0;

  m_Compound = NPL::Particle(CompoundName);
  m_Compound.SetExcitationEnergy(MEx);
  m_Compound.SetKineticEnergy(MEK);
  if(m_FissionModelName=="GEF"){
    if(m_FissionModel->IsValid(m_Compound.GetZ(), m_Compound.GetA())){
      worked=true;
      m_FissionModel->InitCompound(MEx,MEK,Lfis,Theta,Phi);
      m_FissionModel->Treat();

      int Ah = m_FissionModel->GetAffh();
      int Zh = m_FissionModel->GetZffh();
      int Al = m_FissionModel->GetAffl();
      int Zl = m_FissionModel->GetZffl();

      double KEl = m_FissionModel->GetKEffl();
      double KEh = m_FissionModel->GetKEffh();
      double Brhol = m_FissionModel->GetBrhoffl();
      double Brhoh = m_FissionModel->GetBrhoffh();

      NPL::Particle FFl = NPL::Particle(Zl,Al);
      NPL::Particle FFh = NPL::Particle(Zh,Ah);
      if(m_VamosChargeStates==1){
        // Include Charge states distribtuion from Baron et al. NIM 328 (1993) 177-182
        FFl.SetBrho(Brhol);
        FFh.SetBrho(Brhoh);
      }
      else{
        FFl.SetKineticEnergy(KEl);
        FFh.SetKineticEnergy(KEh);
      }
      FissionFragments.push_back(FFl);
      FissionFragments.push_back(FFh);
      Ex.push_back(0);
      Ex.push_back(0);
      DEK.push_back(KEl);
      DEK.push_back(KEh);

      double Massl = FFl.Mass();
      double Massh = FFh.Mass();
      double El = KEl+Massl;
      double Eh = KEh+Massh;

      double Pl = sqrt(El*El-Massl*Massl);
      double Ph = sqrt(Eh*Eh-Massh*Massh);

      double Thetal = m_FissionModel->GetThffl();
      double Thetah = m_FissionModel->GetThffh();
      double Phil   = m_FissionModel->GetPhffl();
      double Phih   = m_FissionModel->GetPhffh();

      TVector3 Momentuml = Pl * TVector3(sin(Thetal)*cos(Phil),
          sin(Thetal)*sin(Phil),
          cos(Thetal));

      TVector3 Momentumh = Ph * TVector3(sin(Thetah)*cos(Phih),
          sin(Thetah)*sin(Phih),
          cos(Thetah));

      DPx.push_back(Momentuml.X());
      DPx.push_back(Momentumh.X());
      DPy.push_back(Momentuml.Y());
      DPy.push_back(Momentumh.Y());
      DPz.push_back(Momentuml.Z());
      DPz.push_back(Momentumh.Z());

      TKE = m_FissionModel->GetTKE();
      KE1 = m_FissionModel->GetKE1();
      KE2 = m_FissionModel->GetKE2();

      // Neutron and gamma emission
      float* En1;
      float* Eg1;
      En1 = m_FissionModel->GetNeutronEnergyFrag1();
      //cout << "----- Neutron energy: " << endl;
      //for(int i=0; i<51; i++) {
      //  cout << "En= " << En1[i] << endl;
      //}

      Eg1 = m_FissionModel->GetGammaEnergyFrag1();
      //cout << "----- Gamma energy: " << endl;
      //for(int i=0; i<101; i++){
      //  cout << "Eg= " << Eg1[i] << endl;
      //}
    }
  }
  return worked;
}

























