/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  All detector added in the project should derive from this virtual class  *
 *  A vector ofNPS::VDetector object is manage in the DetectorConstruction   *
 *  class and call the virtual method of this class implemented in the       *
 *  daughter class object.                                                   *
 *  This inheritance insure homogeneity and modularity of the code           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "NPSVDetector.hh"
#include "RootOutput.h"
#include "G4SDManager.hh"
#include <iostream>
 using namespace std;
TInteractionCoordinates* NPS::VDetector::ms_InterCoord = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Constructor
NPS::VDetector::VDetector(){
    if (ms_InterCoord == 0) ms_InterCoord = new TInteractionCoordinates();
    InitializeRootOutput();
}


// Destructor
NPS::VDetector::~VDetector(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NPS::VDetector::InitializeRootOutput(){
    TTree* Tree = RootOutput::getInstance()->GetTree();
    if(!Tree->FindBranch("InteractionCoordinates")){
        Tree->Branch("InteractionCoordinates", "TInteractionCoordinates", &ms_InterCoord);
    }
    Tree->SetBranchAddress("InteractionCoordinates", &ms_InterCoord);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4MultiFunctionalDetector* NPS::VDetector::CheckScorer(string name,bool &exist){
    exist = true;
    G4MultiFunctionalDetector* ptr = NULL;
    
    ptr =
    (G4MultiFunctionalDetector*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(name,false);

if(!ptr){
      ptr = new G4MultiFunctionalDetector(name.c_str());
      exist = false;
    }
    return ptr;
}
