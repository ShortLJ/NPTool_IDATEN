#ifndef ROOTOUTPUT_HH
#define ROOTOUTPUT_HH
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 21/07/09                                                 *
 * Last update    : 10/05/21                                                 *
 *---------------------------------------------------------------------------*
 * Decription: This class is a singleton class which deals with the ROOT     *
 *             output file and tree both for NPSimulation and NPAnalysis.    *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *   + 03/02/11: Add support for TAsciiFile objects (N. de Sereville)        *
 *   + 10/05/21: Add support for split tree output (A. Matta)                *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <map>

// NPL headers
#include "TAsciiFile.h"

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include <string>

class RootOutput{
public:
  // The analysis class is designed to be a singleton (i.e. only one instance
  // can exist). A member function called Instance is defined, which allows
  // the user to get a pointer to the existing instance or to create it if
  // it does not yet exist:
  // (see the constructor for an explanation of the arguments)
  static RootOutput* getInstance(std::string fileNameBase = "Simulation",
                                 std::string treeNameBase = "SimulatedTree",
                                 bool split = false);
  
  // The analysis class instance can be deleted by calling the Destroy
  // method (NOTE: The class destructor is protected, and can thus not be
  // called directly):
  static void Destroy();
  
protected:
  // Constructor (protected)
  RootOutput(std::string fileNameBase, std::string treeNameBase, bool split=false);
 
  // Destructor (protected)
  virtual ~RootOutput();
  
  // Prevent copying
  RootOutput(const RootOutput& only);
  const RootOutput& operator=(const RootOutput& only);
  
private:
  // The static instance of the RootOutput class:
  static RootOutput* instance;
  
private:
  void InitAsciiFiles();
  void CreateTreeAndFile(std::string name);
public:
  TFile*      GetFile(std::string name="global") ; 
  TTree*      GetTree(std::string name="global") ;
  TList*      GetList(std::string name="global") ;
  TAsciiFile* GetAsciiFileEventGenerator()        {return pEventGenerator;}
  TAsciiFile* GetAsciiFileDetectorConfiguration() {return pDetectorConfiguration;}
  TAsciiFile* GetAsciiFileCalibration()           {return pCalibrationFile;}
  TAsciiFile* GetAsciiFileRunToTreat()            {return pRunToTreatFile;}
  TAsciiFile* GetAsciiFileAnalysisConfig()        {return pAnalysisConfigFile;}
  void        Fill();
  
private:
  std::string                   pBaseName;
  std::string                   pTreeName;
  std::string                   pMasterFile;
  TDirectory*                   pCurrentDirectory;
  bool                          pSplit;
  std::map<std::string, TFile*> pRootFiles;
  std::map<std::string, TTree*> pRootTrees;
  std::map<std::string, TList*> pRootLists;
  TAsciiFile* pEventGenerator;
  TAsciiFile* pDetectorConfiguration;
  TAsciiFile* pCalibrationFile;
  TAsciiFile* pRunToTreatFile;
  TAsciiFile* pAnalysisConfigFile;
  
};

#endif // ROOTOUTPUT_HH
