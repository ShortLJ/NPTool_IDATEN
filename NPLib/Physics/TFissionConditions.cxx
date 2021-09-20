/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  address: pierre.morfouace2@cea.fr      *
 *                                                                           *
 * Creation Date  : 01/10/20                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: This class records all the information concerning the fission *
 *             fragment and the Ex of the fissioning system                  *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
using namespace std;

#include "TFissionConditions.h"

ClassImp(TFissionConditions)
  ////////////////////////////////////////////////////////////////////////////////
  TFissionConditions::TFissionConditions(){
  }
////////////////////////////////////////////////////////////////////////////////
TFissionConditions::~TFissionConditions(){
}
////////////////////////////////////////////////////////////////////////////////
void TFissionConditions::Clear(){
  // Fissionning system
  fFC_A_CN = -10;
  fFC_Z_CN = -10;
  fFC_Ex_CN = -10;
  fFC_ELab_CN = -10;
  fFC_ThetaLab_CN = -10;

  // Fission process
  fFC_TKE = -10;
  fFC_KE1 = -10;
  fFC_KE2 = -10;
  fFC_Neutron_Multiplicity = -10;

  // Fission Fragments
  fFC_Fragment_Name.clear();
  fFC_Fragment_A.clear();
  fFC_Fragment_Z.clear();
  fFC_Fragment_Theta.clear();
  fFC_Fragment_Phi.clear();
  fFC_Fragment_Kinetic_Energy.clear();
  fFC_Fragment_Brho.clear();
  fFC_Fragment_Momentum_Direction_X.clear();
  fFC_Fragment_Momentum_Direction_Y.clear();
  fFC_Fragment_Momentum_Direction_Z.clear();
}
////////////////////////////////////////////////////////////////////////////////
void TFissionConditions::Dump() const{
  cout << "--------- Fission Condition Dump ---------" << endl ;

  // Fissionning system
  cout << "\t ---- Fissionning system ---- " << endl;
  cout << "\t Z:  " << fFC_Z_CN << endl;
  cout << "\t A: " << fFC_A_CN << endl;
  cout << "\t Ex: " << fFC_Ex_CN << endl;
  cout << "\t TKE: " << fFC_TKE << endl;
  cout << "\t Neutron Mult: " << fFC_Neutron_Multiplicity << endl;


  // Fission fragments
  unsigned int size = fFC_Fragment_A.size();
  for(unsigned int i = 0 ; i < size; i ++){
    cout << "\t ---- Particle " << i << " ---- " << endl;
    cout << "\t Fragment Z: " <<   fFC_Fragment_Z[i] << endl;
    cout << "\t Fragment A: " <<   fFC_Fragment_A[i] << endl;
  }


}

