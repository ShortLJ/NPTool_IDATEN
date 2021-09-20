/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 10/06/09                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: This class is an Abstract Base Class (ABC) from which should  *
 *             derive all different modules from the Gaspard tracker.        *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "GaspardTrackerModule.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "RootOutput.h"


TGaspardTrackerData *GaspardTrackerModule::ms_Event = 0;



GaspardTrackerModule::GaspardTrackerModule()
{
   if (ms_Event == 0) ms_Event = new TGaspardTrackerData();

   InitializeRootOutput();
   InitializeIndex();
   InitializeMaterial();
}



GaspardTrackerModule::~GaspardTrackerModule()
{
}



void GaspardTrackerModule::InitializeRootOutput()
{
   RootOutput *pAnalysis = RootOutput::getInstance();
   TTree *pTree = pAnalysis->GetTree();
   // if the branch does not exist yet, create it
   if (!pTree->FindBranch("GASPARD"))
      pTree->Branch("GASPARD", "TGaspardTrackerData", &ms_Event);
    pTree->SetBranchAddress("GASPARD", &ms_Event);

}



void GaspardTrackerModule::InitializeIndex()
{
   m_index["Square"]     =    0;
   m_index["Rectangle"]  =    0;
   m_index["Trapezoid"]  =  100;
   m_index["Annular"]    =  200;
   m_index["DummyShape"] = 1000;
}



void GaspardTrackerModule::InitializeMaterial()
{
   m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
   m_MaterialVacuum  = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}
