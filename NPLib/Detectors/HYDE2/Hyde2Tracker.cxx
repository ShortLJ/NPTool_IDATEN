/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 31/08/09                                                 *
 * Last update    : 10/09/09                                                 *
 *---------------------------------------------------------------------------*
 * Decription: This class is mainly an interface to the                      *
 *             THyde2TrackerPhysics class and it deals with the geometrical*
 *             correspondance between strip number and absolute coordinates  *
 *             (X, Y, Z) of interaction.                                     *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *   + 10/09/09: Add support for Square and DummyShape shapes                *
 *               (N. de Sereville)                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "Hyde2Tracker.h"

// C++ headers
#include <iostream>   
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>

// NPL headers
#include "RootInput.h"
#include "RootOutput.h"

// ROOT headers
#include "TChain.h"

// Hyde2 headers
#include "Hyde2TrackerTrapezoid1.h"
#include "Hyde2TrackerTrapezoid2.h"
#include "Hyde2TrackerSquare1.h"

using namespace std ;   




Hyde2Tracker::Hyde2Tracker()
   : m_EventData (new THyde2TrackerData()),
     m_EventPhysics (new THyde2TrackerPhysics())
{
}   



Hyde2Tracker::~Hyde2Tracker()
{
   delete m_EventData;
   delete m_EventPhysics;
}   



void Hyde2Tracker::Clear()
{
}



// Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
void Hyde2Tracker::ReadConfiguration(string Path)    
{
   // open configuration file
   ifstream ConfigFile;
   ConfigFile.open(Path.c_str());

   bool HYD2TrkTrapezoid1  = false;
   bool HYD2TrkTrapezoid2  = false;
   bool HYD2TrkSquare1     = false;

   string LineBuffer;
   while (!ConfigFile.eof()) {
      getline(ConfigFile, LineBuffer);
      if (LineBuffer.compare(0, 11, "HYD2Square1") == 0  &&  HYD2TrkSquare1 == false) {
         HYD2TrkSquare1 = true;

         // instantiate a new "detector" corresponding to the Shape elements
         Hyde2TrackerModule* myDetector = new Hyde2TrackerSquare1(m_ModulesMap, m_EventPhysics);

         // Pass the data object to the Hyde2Tracker*** object
         myDetector->SetHyde2DataPointer(m_EventData);

         // read part of the configuration file corresponding to shape elements
         ConfigFile.close();
         myDetector->ReadConfiguration(Path);
         ConfigFile.open(Path.c_str());
      }
      else if (LineBuffer.compare(0, 14, "HYD2Trapezoid1") == 0  &&  HYD2TrkTrapezoid1 == false) {
         HYD2TrkTrapezoid1 = true;

         // instantiate a new "detector" corresponding to the Trapezoid elements
         Hyde2TrackerModule* myDetector = new Hyde2TrackerTrapezoid1(m_ModulesMap, m_EventPhysics);

         // Pass the data object to the Hyde2Tracker*** object
         myDetector->SetHyde2DataPointer(m_EventData);

         // read part of the configuration file corresponding to trapezoid elements
         ConfigFile.close();
         myDetector->ReadConfiguration(Path);
         ConfigFile.open(Path.c_str());
      }
      else if (LineBuffer.compare(0, 14, "HYD2Trapezoid2") == 0  &&  HYD2TrkTrapezoid2 == false) {
         HYD2TrkTrapezoid2 = true;

         // instantiate a new "detector" corresponding to the Trapezoid elements
         Hyde2TrackerModule* myDetector = new Hyde2TrackerTrapezoid2(m_ModulesMap, m_EventPhysics);

         // Pass the data object to the Hyde2Tracker*** object
         myDetector->SetHyde2DataPointer(m_EventData);

         // read part of the configuration file corresponding to trapezoid elements
         ConfigFile.close();
         myDetector->ReadConfiguration(Path);
         ConfigFile.open(Path.c_str());
      }
   }
}


// Read stream at Path and pick-up calibration parameter using Token
// If argument is "Simulation" no change calibration is loaded
void Hyde2Tracker::ReadCalibrationFile(string Path)   
{
   Path.c_str();   // remove warning at compilation time

/*   // Order of Polynom function used for calibration
   int Calibration_Si_E_Order;
   int Calibration_Si_T_Order;
   int Calibration_SiLi_E_Order;
   int Calibration_CsI_E_Order;

   // Calibration_Si_X_E[DetectorNumber][StripNumber][Order of Coeff]
   vector< vector< vector< double > > >   Calibration_Si_X_E   ;
   vector< vector< vector< double > > >   Calibration_Si_X_T   ;
   vector< vector< vector< double > > >   Calibration_Si_Y_E   ;
   vector< vector< vector< double > > >   Calibration_Si_Y_T   ;

   // Calibration_SiLi_E[DetectorNumber][PadNumber][Order of Coeff]
   vector< vector< vector< double > > >   Calibration_SiLi_E   ;

   // Calibration_SiLi_E[DetectorNumber][CrystalNumber][Order of Coeff]
   vector< vector< vector< double > > >   Calibration_CsI_E   ;

   if (Path == "Simulation") {   // Simulation case: data already calibrated
      Calibration_Si_E_Order   = 1;
      Calibration_Si_T_Order   = 1;
      Calibration_SiLi_E_Order = 1;
      Calibration_CsI_E_Order  = 1;

      vector<double> Coef;
      // Order 0            Order 1
      Coef.push_back(0) ; Coef.push_back(1)    ;

      vector< vector<double> > StripLine       ;
      StripLine.resize( 128 , Coef)         ;

      Calibration_Si_X_E.resize( m_NumberOfModule , StripLine)   ;
      Calibration_Si_X_T.resize( m_NumberOfModule , StripLine)   ;
      Calibration_Si_Y_E.resize( m_NumberOfModule , StripLine)   ;
      Calibration_Si_Y_T.resize( m_NumberOfModule , StripLine)   ;
            
      Calibration_SiLi_E.resize( m_NumberOfModule , StripLine)   ;
      Calibration_CsI_E .resize( m_NumberOfModule , StripLine)   ;
   }
   else {
   }
*/
}      

   

// Activated associated Branches and link it to the private member DetectorData address
// In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
void Hyde2Tracker::InitializeRootInput()       
{
   TChain* inputChain = RootInput::getInstance()->GetChain();
   inputChain->SetBranchStatus("HYDE2", true);
   inputChain->SetBranchStatus("fHYD2*", true);
   inputChain->SetBranchAddress("HYDE2", &m_EventData);
}



// Create associated branches and associated private member DetectorPhysics address
void Hyde2Tracker::InitializeRootOutput()    
{
   TTree* outputTree = RootOutput::getInstance()->GetTree();
   outputTree->Branch("HYDE2", "THyde2TrackerPhysics", &m_EventPhysics);
}



// This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
void Hyde2Tracker::BuildPhysicalEvent()      
{
   if (m_EventData->GetHYD2TrkFirstStageFrontEMult() > 0) {
      for (UInt_t i = 0; i < m_EventData->GetHYD2TrkFirstStageFrontEMult(); i++) {
         UShort_t detecNbr = m_EventData->GetHYD2TrkFirstStageFrontEDetectorNbr(i);
         m_ModulesMap[detecNbr]->BuildPhysicalEvent();
      }
   }
}



// Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
// This method aimed to be used for analysis performed during experiment, when speed is requiered.
// NB: This method can eventually be the same as BuildPhysicalEvent.
void Hyde2Tracker::BuildSimplePhysicalEvent()   
{
   BuildPhysicalEvent();
}



void Hyde2Tracker::DumpModulesMap()
{
   cout << "Hyde2Tracker::DumpMap()" << endl;

   map<int, Hyde2TrackerModule*>::iterator it;

   for (it = m_ModulesMap.begin(); it != m_ModulesMap.end(); ++it) {
      cout << "Module " << it->first << " -> Address " << it->second << endl;
   }
}



void Hyde2Tracker::DumpStrippingScheme(int moduleNumber)
{
   cout << "Hyde2Tracker::DumpStrippingScheme()" << endl;

   for (int i = 1; i < 128; i++) {   // front part
      for (int j = 1; j < 128; j++) {   // back part
         cout << "strips X, Y: " << i << "  " << j << "\t--->\t (X,Y,Z) mm: "
              << m_ModulesMap[moduleNumber]->GetStripPositionX(moduleNumber, i, j) << "\t"
              << m_ModulesMap[moduleNumber]->GetStripPositionY(moduleNumber, i, j) << "\t"
              << m_ModulesMap[moduleNumber]->GetStripPositionZ(moduleNumber, i, j) << endl;
      }
   }
}



double Hyde2Tracker::GetEnergyDeposit()      
{ 
   if (m_EventPhysics->GetEventMultiplicity() > 0) {
      return m_EventPhysics->GetTotalEnergy(0); 
   }
   else
      return -1000;
}



TVector3 Hyde2Tracker::GetPositionOfInteraction()
{
   TVector3 Position = TVector3(-1000,-1000,-1000);

   if (m_EventPhysics->GetEventMultiplicity() == 1) {
      Int_t ModuleNbr = m_EventPhysics->GetModuleNumber(0);
      Int_t FrontPos  = m_EventPhysics->GetFirstStageFrontPosition(0);
      Int_t BackPos   = m_EventPhysics->GetFirstStageBackPosition(0);
      Position = TVector3(m_ModulesMap[ModuleNbr]->GetStripPositionX(ModuleNbr, FrontPos, BackPos),
                          m_ModulesMap[ModuleNbr]->GetStripPositionY(ModuleNbr, FrontPos, BackPos),
                          m_ModulesMap[ModuleNbr]->GetStripPositionZ(ModuleNbr, FrontPos, BackPos));
   }

   return Position;
}
