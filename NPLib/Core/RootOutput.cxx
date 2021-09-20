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

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "RootOutput.h"
#include "NPOptionManager.h"

using namespace std;

RootOutput* RootOutput::instance = 0;
////////////////////////////////////////////////////////////////////////////////
RootOutput* RootOutput::getInstance(std::string fileNameBase, std::string treeNameBase,bool split){
  // A new instance of RootOutput is created if it does not exist:
  if (instance == 0) {
    instance = new RootOutput(fileNameBase.c_str(), treeNameBase.c_str(),split);
  }

  // The instance of RootOutput is returned:
  return instance;
}

////////////////////////////////////////////////////////////////////////////////
void RootOutput::Destroy(){ 
  if (instance != 0) {
    delete instance;
    instance = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
RootOutput::RootOutput(std::string fileNameBase, std::string treeNameBase,bool split){

  pSplit = (split || NPOptionManager::getInstance()->IsSplit());
  cout << endl << "/////////// ROOT Output files ///////////" << endl;
  cout << "Initializing ouput trees and files ";
  if(pSplit)
    cout << "in split mode (one tree per detector)" << endl;
  else
    cout << endl;
  pTreeName=treeNameBase;
  pCurrentDirectory= gDirectory;
  bool analysis=false;
  bool simulation=false;

  if(NPOptionManager::getInstance()->IsAnalysis()){
    analysis = true;
  }
  else if(NPOptionManager::getInstance()->IsSimulation()){
    simulation= true;
  }

  // Setup the base name
  if(analysis)
    pBaseName = NPOptionManager::getInstance()->GetAnalysisOutputPath();
  else if(simulation)
    pBaseName = NPOptionManager::getInstance()->GetSimulationOutputPath();
  else
    pBaseName="./";

  pBaseName += "/"+fileNameBase;

  if (fileNameBase.find("root")==std::string::npos) 
    pBaseName += ".root";

  // removing "//" from path
//  while(size_t pos = pBaseName.find("//")!=std::string::npos){
//    pBaseName.erase(pos);
//  }

  if(pSplit){
    // Create a folder for all the trees
    string stripname = pBaseName;
    stripname.erase(stripname.find(".root"),5);
    string path = stripname.substr(0,stripname.rfind("/"));
    string filebase = stripname.substr(stripname.rfind("/")+1);
    string command = "mkdir -p "+stripname;
    int res = system(command.c_str());
    if(res!=0){
      std::cout << "Error creating folder " << stripname << std::endl;
      exit(1);
    }
    // create the master file 
    pMasterFile=stripname+"/"+filebase+".tree";
    ofstream master(pMasterFile.c_str(),std::ofstream::trunc);
    master.close();
  }
  else{
    CreateTreeAndFile("global");
  }
  /////
  // Create the last file 
  if(treeNameBase=="SimulatedTree"){
    string path = getenv("NPTOOL");
    path+="/.last_sim_file";
    ofstream last_sim_file(path.c_str());
    last_sim_file << "Tree "<< pTreeName <<endl
      << " " << pBaseName <<endl;
    last_sim_file.close();
  }

  else if(treeNameBase=="PhysicsTree"){
    string path = getenv("NPTOOL");
    path+="/.last_phy_file";
    ofstream last_phy_file(path.c_str());
    last_phy_file << "Tree "<< pTreeName<<endl
      << " " << pBaseName <<endl;
    last_phy_file.close();
  }

  else if(treeNameBase=="ResultTree"){
    string path = getenv("NPTOOL");
    path+="/.last_res_file";
    ofstream last_res_file(path.c_str());
    last_res_file << "Tree " << pTreeName << endl 
      << " " << pBaseName<<endl;
    last_res_file.close();
  }

  else{
    string path = getenv("NPTOOL");
    path+="/.last_any_file";
    ofstream last_any_file(path.c_str());
    last_any_file << "Tree " << pTreeName <<endl
      << " " << pBaseName<< endl;
    last_any_file.close();
  }


  InitAsciiFiles();
}

////////////////////////////////////////////////////////////////////////////////
void RootOutput::CreateTreeAndFile(std::string name){
  // Create the tree only if does not exist already
  string file_name=pBaseName;

  if(pRootFiles.find(name)==pRootFiles.end()){
    if(pSplit){
      string  strip= pBaseName.substr(pBaseName.rfind("/"));
      strip = strip.substr(0,strip.rfind(".root"));
      string  insertion= "_"+name;
      file_name.insert(file_name.rfind(".root"),insertion);
      file_name.insert(file_name.rfind("/"),strip);
    }
    cout << " - Creating output file " << file_name.c_str() << endl;
    pRootFiles[name] = new TFile(file_name.c_str(), "RECREATE");
    pRootTrees[name] = new TTree(pTreeName.c_str(), "Data created / analysed with the nptool package");
    pRootFiles[name]->SetCompressionLevel(1);

    // Init TAsciiFile objects
    gDirectory->cd(pCurrentDirectory->GetPath()); 

    if(NPOptionManager::getInstance()->GetCircularTree()){
      cout << "Information: Output tree is set to circular mode" << endl;
      pRootTrees[name]->SetCircular(1000); 
    }

    if(pSplit){// Add the tree to the .tree file
      ofstream master(pMasterFile.c_str(),std::ofstream::app);
      file_name = file_name.substr(file_name.rfind("/")+1);
      master << file_name.c_str() << endl;
      master.close();
    }
  } 
  //
  /*
     for(auto it = m_DetectorMap.begin() ; it!=m_DetectorMap.end() ;++it){
     string insertion = "_"+it->first;
     master << filebase << insertion << ".root" << std::endl;
     string filename=path+"/"+filebase+"/"+filebase+insertion+".root";
     auto file = new TFile(filename.c_str(),"RECREATE");
     string treename = "RawTree_"+it->first;
     auto tree = new TTree("RawTree",treename.c_str());
     m_TreeMap[it->first]=tree;
     m_FileMap[it->first]=file;
     tree->SetDirectory(file);
     std::cout << "Splitting tree: " << filename << std::endl;
     it->second->InitBranch(tree);
     }
     master.close();

*/
}
////////////////////////////////////////////////////////////////////////////////
void RootOutput::Fill(){
 for(auto it = pRootTrees.begin();it!=pRootTrees.end();it++){
  it->second->Fill();
 } 
}

////////////////////////////////////////////////////////////////////////////////
void RootOutput::InitAsciiFiles(){
  // get NPOptionManager pointer
  NPOptionManager* OptionManager  = NPOptionManager::getInstance();

  // Event generator
  // Get file name from NPOptionManager
  std::string fileNameEG = OptionManager->GetReactionFile();
  pEventGenerator = new TAsciiFile();
  pEventGenerator->SetNameTitle("EventGenerator", fileNameEG.c_str());
  pEventGenerator->Append(fileNameEG.c_str());
  //pEventGenerator->Write(0,TAsciiFile::kOverwrite);

  // Detector configuration 
  // Get file name from NPOptionManager
  std::string fileNameDC = OptionManager->GetDetectorFile();
  pDetectorConfiguration = new TAsciiFile();
  pDetectorConfiguration->SetNameTitle("DetectorConfiguration", fileNameDC.c_str());
  pDetectorConfiguration->Append(fileNameDC.c_str());
  //pDetectorConfiguration->Write(0,TAsciiFile::kOverwrite);

  // Run to treat file
  // Get file name from NPOptionManager
  pRunToTreatFile = new TAsciiFile();
  if (!OptionManager->IsDefault("RunToTreat")) {
    std::string fileNameRT = OptionManager->GetRunToReadFile();
    pRunToTreatFile->SetNameTitle("RunToTreat", fileNameRT.c_str());
    pRunToTreatFile->Append(fileNameRT.c_str());
    //pRunToTreatFile->Write(0,TAsciiFile::kOverwrite);
  }

  // Calibration files
  pCalibrationFile = new TAsciiFile();
  if (!OptionManager->IsDefault("Calibration")) {
    std::string fileNameCal = OptionManager->GetCalibrationFile();
    pCalibrationFile->SetNameTitle("Calibration", fileNameCal.c_str());
    //pCalibrationFile->Write(0,TAsciiFile::kOverwrite);
  }

  // Analysis configuration files
  pAnalysisConfigFile = new TAsciiFile();
  pAnalysisConfigFile->SetNameTitle("AnalysisConfig", "AnalysisConfig");
  //pAnalysisConfigFile->Write(0,TAsciiFile::kOverwrite);
}

////////////////////////////////////////////////////////////////////////////////
RootOutput::~RootOutput(){ 
  // The data is written to the file and the tree is closed:
  if (pRootFiles.size()>0) {
    cout << endl << endl << "Root Output summary" << endl;
    TDirectory* pCurrentDirectory= gDirectory;
    for(auto it = pRootFiles.begin(); it!=pRootFiles.end();it++){
      cout << " - " <<it->first << " tree and file " << endl;
      gDirectory->cd(it->second->GetPath());
      // write TAsciiFile if used
      // EventGenerator
      if (!pEventGenerator->IsEmpty()) pEventGenerator->Write(0,TAsciiFile::kOverwrite);
      // DetectorConfiguration
      if (!pDetectorConfiguration->IsEmpty()) pDetectorConfiguration->Write(0,TAsciiFile::kOverwrite);
      // CalibrationFile
      if (!pCalibrationFile->IsEmpty()) pCalibrationFile->Write(0,TAsciiFile::kOverwrite);
      // RunToTreatFile
      if (!pRunToTreatFile->IsEmpty()) pRunToTreatFile->Write(0,TAsciiFile::kOverwrite);
      // Analysis ConfigFile
      if (!pAnalysisConfigFile->IsEmpty()) pAnalysisConfigFile->Write(0,TAsciiFile::kOverwrite);

      cout << "  -> Number of entries in the " << it->first << " Tree: " << pRootTrees[it->first]->GetEntries() << endl;
      cout << "  -> Number of bites written to file: " << pRootTrees[it->first]->Write(0, TObject::kOverwrite) << endl;
      it->second->Flush();
      it->second->Purge(1);

      gDirectory->cd(pCurrentDirectory->GetPath());
      it->second->Close();
    }
  }

  else {
    cout << "\033[1;31mNo histograms and Tree !\033[0m" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
TFile* RootOutput::GetFile(std::string name)  {
  if(!pSplit)
    name="global";

  if(pRootFiles.find(name)!=pRootFiles.end())
    return pRootFiles[name];
  else{
    std::cout << "Error: Requested file for detector " << name << " does not exist" << std::endl;
    exit(1);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
TTree* RootOutput::GetTree(std::string name) {
  if(!pSplit)
    name="global";

  if(pRootTrees.find(name)==pRootTrees.end())
    CreateTreeAndFile(name);
  
  return pRootTrees[name];
}

