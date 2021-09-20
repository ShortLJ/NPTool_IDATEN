#ifndef TFISSIONCONDITIONS_H
#define TFISSIONCONDITIONS_H
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
 *             fragments and the Ex of the fissionning system                *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/



// STL Header
#include <vector>
#include <string>
#include <cmath>
using namespace std ;

// Root Header
#include "TObject.h"
#include "TVector3.h"

// NPTOOL headers
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

class TFissionConditions : public TObject{
  private:

    // Fissionning system
    int fFC_A_CN;
    int fFC_Z_CN;
    double fFC_Ex_CN;
    double fFC_ELab_CN;
    double fFC_ThetaLab_CN;

    // Fission Process
    double fFC_TKE;
    double fFC_KE1;
    double fFC_KE2;
    int fFC_Neutron_Multiplicity;

    // Fission Fragments
    vector<string> fFC_Fragment_Name;
    vector<int> fFC_Fragment_Z;
    vector<int> fFC_Fragment_A;
    vector<double> fFC_Fragment_Theta;
    vector<double> fFC_Fragment_Phi;
    vector<double> fFC_Fragment_Kinetic_Energy;
    vector<double> fFC_Fragment_Brho;
    vector<double> fFC_Fragment_Momentum_Direction_X;
    vector<double> fFC_Fragment_Momentum_Direction_Y;
    vector<double> fFC_Fragment_Momentum_Direction_Z;

  public:
    TFissionConditions();
    virtual ~TFissionConditions();

    void  Clear();
    void  Clear(const Option_t*) {Clear();};
    void  Dump() const;

    /////////////////////           SETTERS           ////////////////////////
    // Fissionning system parameter
    void SetA_CN        (const int A_CN)        {fFC_A_CN  = A_CN;}//!
    void SetZ_CN        (const int Z_CN)        {fFC_Z_CN  = Z_CN;}//!
    void SetEx_CN       (const double Ex_CN)    {fFC_Ex_CN  = Ex_CN;}//!
    void SetELab_CN     (const double ELab_CN)  {fFC_ELab_CN  = ELab_CN;}//!
    void SetThetaLab_CN (const double Theta_CN) {fFC_ThetaLab_CN  = Theta_CN;}//!

    // Fission process
    void Set_TKE (const double E) {fFC_TKE = E;}//!
    void Set_KE1 (const double E) {fFC_KE1 = E;}//!
    void Set_KE2 (const double E) {fFC_KE2 = E;}//!
    void SetNeutronMultiplicity (const int mult) {fFC_Neutron_Multiplicity = mult;}//!

    // Fission Fragments
    void SetFragmentName          (const string name) {fFC_Fragment_Name.push_back(name);}//!
    void SetFragmentZ             (const int Z) {fFC_Fragment_Z.push_back(Z);}//!
    void SetFragmentA             (const int A) {fFC_Fragment_A.push_back(A);}//!
    void SetFragmentTheta         (const double Theta) {fFC_Fragment_Theta.push_back(Theta);}//!
    void SetFragmentPhi           (const double Phi) {fFC_Fragment_Phi.push_back(Phi);}//!
    void SetFragmentKineticEnergy (const double E) {fFC_Fragment_Kinetic_Energy.push_back(E);}//!
    void SetFragmentBrho          (const double Brho) {fFC_Fragment_Brho.push_back(Brho);}//!
    void SetFragmentMomentumX     (const double P) {fFC_Fragment_Momentum_Direction_X.push_back(P);}//!
    void SetFragmentMomentumY     (const double P) {fFC_Fragment_Momentum_Direction_Y.push_back(P);}//!
    void SetFragmentMomentumZ     (const double P) {fFC_Fragment_Momentum_Direction_Z.push_back(P);}//!

    /////////////////////           GETTERS           ////////////////////////
    // Fissionning system parameter
    int GetA_CN() const {return fFC_A_CN;}//!
    int GetZ_CN() const {return fFC_Z_CN;}//!
    double GetEx_CN() const {return fFC_Ex_CN;}//!
    double GetELab_CN() const {return fFC_ELab_CN;}//!
    double GetThetaLab_CN() const {return fFC_ThetaLab_CN;}//!

    // Fission Process
    double GetTKE() const {return fFC_TKE;}//!
    double GetKE1() const {return fFC_KE1;}//!
    double GetKE2() const {return fFC_KE2;}//!
    int GetNeutronMultiplicity() const {return fFC_Neutron_Multiplicity;}//!

    // emmitted particles
    string GetFragmentName (const int &i) const {return fFC_Fragment_Name[i];}//!
    int GetFragmentZ (const int &i) const {return fFC_Fragment_Z[i];}//!
    int GetFragmentA (const int &i) const {return fFC_Fragment_A[i];}//!
    double GetFragmentTheta (const int &i) const {return fFC_Fragment_Theta[i];}//!
    double GetFragmentPhi (const int &i) const {return fFC_Fragment_Phi[i];}//!
    double GetFragmentKineticEnergy (const int &i) const {return fFC_Fragment_Kinetic_Energy[i];}//!
    double GetFragmentBrho (const int &i) const {return fFC_Fragment_Brho[i];}//!
    double GetFragmentMomentumX (const int &i) const {return fFC_Fragment_Momentum_Direction_X[i];}//!
    double GetFragmentMomentumY (const int &i) const {return fFC_Fragment_Momentum_Direction_Y[i];}//!
    double GetFragmentMomentumZ (const int &i) const {return fFC_Fragment_Momentum_Direction_Z[i];}//!

    ClassDef(TFissionConditions, 1) // TFissionConditions structure
};

#endif
