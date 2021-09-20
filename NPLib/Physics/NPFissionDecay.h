#ifndef NPFISSIONDECAY_H
#define NPFISSIONDECAY_H
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
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>
#include <vector>
#include <set>

// NPL Header
#include "NPInputParser.h"
#include "NPParticle.h"
#include "GEF.h"

namespace NPL{
  class FissionDecay{
    public: 
        FissionDecay(){};
        FissionDecay(std::string compound, std::string fission_model);
        ~FissionDecay(){};

    public:
        void ReadConfiguration(std::string Path);
        void ReadConfiguration(NPL::InputParser parser);
        int GetFissionToken() {return HasFissionToken;}
    private:
        GEF* m_FissionModel;
        std::string m_FissionModelName;
        std::string m_CompoundName; 
        NPL::Particle m_Compound;
        double m_MotherMass;
        std::vector<std::string> m_FissionFragmentName;
        std::vector<NPL::Particle> m_FissionFragment;
        std::vector<double> m_FissionFragmentMasses;
        int HasFissionToken; 
        bool m_VamosChargeStates;
        bool m_shoot_FF;
        bool m_shoot_neutron;
        bool m_shoot_gamma;

    public:
        // Given Energy and Momentum direction of the compound,
        // Send back Momemtum and Energy of fission fragments
        // Return false if the fission is not possible
        bool GenerateEvent(string CompoundName, double MEx,double MEK,double MPx, double MPy,double MPz,
        std::vector<NPL::Particle>& FissionFragments, std::vector<double>& Ex,
        std::vector<double>& DEK,
        std::vector<double>& DPx,std::vector<double>& DPy,std::vector<double>& DPz,
        double& TKE, double& KE1, double& KE2);
   
   public:// Getter
        inline std::vector<std::string> GetFissionFragmentName() {return m_FissionFragmentName;};
        inline std::vector<NPL::Particle> GetFissionFragment() {return m_FissionFragment;};
        inline NPL::Particle GetCompound() {return m_Compound;};
        inline std::string GetCompoundName() {return m_CompoundName;};
    };
  }
#endif
