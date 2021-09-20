#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : march 2015                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "NPVAnalysis.h"
#include "TMust2Physics.h"
#include "TSSSDPhysics.h"
#include "TInitialConditions.h"
#include "TReactionConditions.h"
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "TRandom3.h"

#include"NPDetectorManager.h"
class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();
    void InitOutputBranch();
    void InitInputBranch();
    void ReInitValue();
    static NPL::VAnalysis* Construct();

  private:
    double Ex;
    double ELab;
    double ThetaLab;
    double ThetaCM;
    double BeamEnergy;
    double OriginalELab;
    double OriginalThetaLab;
    double OriginalBeamEnergy;
    double ReactionVertexX;
    double ReactionVertexY;
    double ReactionVertexZ;
    NPL::Reaction* He10Reaction;

    // intermediate variable
    TRandom3 Rand;
    int DetectorNumber;
    double ThetaNormalTarget;
    double ThetaM2Surface; 
    double X_M2;
    double Y_M2;
    double Z_M2;
    double Si_E_M2;
    double CsI_E_M2; 
    double E_SSSD;
    double Energy ;
    double E_M2 ;
    double Si_X_M2;
    double Si_Y_M2;
    double TargetThickness;

    NPL::EnergyLoss He3CD2  ;
    NPL::EnergyLoss He3Al   ;
    NPL::EnergyLoss He3Si   ;
    NPL::EnergyLoss Li11CD2 ;

    TMust2Physics* M2;
    TSSSDPhysics* SSSD;
    TInitialConditions* Initial;
    TReactionConditions* ReactionConditions;
};
#endif
