#ifndef NPOPTIONMANAGER_HH
#define NPOPTIONMANAGER_HH

/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: A. MATTA  contact address: a.matta@surrey.ac.uk*
 *                                                                           *
 * Creation Date  : 21/07/09                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: This class is a singleton class which deals with input        *
 * arguments of the different NPTool programm (NPS and NPA)                  *
 *---------------------------------------------------------------------------*
 * Comment: The singleton form allow users to call the object from anywhere  *
 * in the code                                                               *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <iostream>
#include <string>
#include <set>

class NPOptionManager{
   public:
      // The analysis class is designed to be a singleton (i.e. only one instance
      // can exist). A member function called Instance is defined, which allows
      // the user to get a pointer to the existing instance or to create it if
      // it does not yet exist:
      // (see the constructor for an explanation of the arguments)
      static NPOptionManager* getInstance(int argc = 0, char** argv = NULL);
  
      // Added for compatibility with pROOF 
      static NPOptionManager* getInstance(std::string arg);

      // The analysis class instance can be deleted by calling the Destroy
      // method (NOTE: The class destructor is protected, and can thus not be
      // called directly):
      static void Destroy();

   protected:
      // Constructor (protected)
      NPOptionManager(int argc, char** argv);
      NPOptionManager(std::string arg);

      // Destructor (protected)
      ~NPOptionManager() {};

      // Prevent copying
      NPOptionManager(const NPOptionManager& only);
      const NPOptionManager& operator=(const NPOptionManager& only);

  private:
    // Read the input argument
    void ReadTheInputArgument(int argc = 0, char** argv = NULL);
    // Look for and Read the project config file if existing
    void ReadProjectConfigFile();

   private:
      // The static instance of the NPOptionManager class:
      static NPOptionManager* instance;

   private:
      void DisplayHelp();
      void DisplayVersion();
      void CheckArguments();
      void CheckEventGenerator();
      void CheckDetectorConfiguration();
      void CheckG4Macro();
      void CreateRunToTreatFile(std::string file, std::string tree );
   public:
      bool IsDefault(const char* type) const;
      void SendErrorAndExit(const char* type) const;

   public:
      // Getters
      // default values
      std::string GetDefaultReactionFile()     {return fDefaultReactionFileName;}
      std::string GetDefaultDetectorFile()     {return fDefaultDetectorFileName;}
      std::string GetDefaultRunToReadFile()    {return fDefaultRunToReadFileName;}
      std::string GetDefaultCalibrationFile()  {return fDefaultCalibrationFileName;}
      std::string GetDefaultOutputFile()       {return fDefaultOutputFileName;}
      std::string GetDefaultG4MacroPath()      {return fDefaultG4MacroPath;}

      // assigned values
      std::string GetReactionFile()             {return fReactionFileName;}
      std::string GetDetectorFile()             {return fDetectorFileName;}
      std::string GetRunToReadFile()            {return fRunToReadFileName;}
      std::string GetCalibrationFile()          {return fCalibrationFileName;}
      std::string GetOutputFile()               {return fOutputFileName;}
      std::string GetOutputTreeName()           {return fOutputTreeName;}
      std::string GetG4MacroPath()              {return fG4MacroPath;}     

      bool   GetDisableAllBranchOption()   {return fDisableAllBranchOption;}
      bool   GetInputPhysicalTreeOption()  {return fInputPhysicalTreeOption;}
      bool   GetGenerateHistoOption()      {return fGenerateHistoOption;}
      bool   GetCheckHistoOption()         {return fCheckHistoOption;}
      bool   GetOnline()                   {return fOnline;}
      bool   GetG4BatchMode()              {return fG4BatchMode;}
      bool   GetCircularTree()             {return fCircularTree;}
      bool   IsAnalysis()                  {return fIsAnalysis;};
      bool   IsSimulation()                {return fIsSimulation;}
      bool   IsSplit()                     {return fIsSplit;}

      int    GetVerboseLevel()             {return fVerboseLevel;}
      int    GetNumberOfEntryToAnalyse()   {return fNumberOfEntryToAnalyse;} 
      int    GetFirstEntryToAnalyse()      {return fFirstEntryToAnalyse;} 
      int    GetSpectraServerPort()        {return fSpectraServerPort;}
      int    GetRandomSeed()               {return fRandomSeed;}
      int    GetRecordTrack()              {return fRecordTrack;}
      std::string GetSharedLibExtension()       {return fSharedLibExtension;}     
      std::string GetLastFile();                 
      std::string GetAnalysisOutputPath(){return m_AnalysisOutputPath;};
      std::string GetSimulationOutputPath(){return m_SimulationOutputPath;};
      std::string GetEnergyLossPath(){return m_EnergyLossPath;};
      // Setters
      void SetReactionFile(const std::string& name)  {fReactionFileName = name;CheckEventGenerator();}
      void SetDetectorFile(const std::string& name)  {fDetectorFileName = name;CheckDetectorConfiguration();}
      void SetRunToReadFile(const std::string& name) {fRunToReadFileName = name;}
      void SetVerboseLevel(int VerboseLevel)         {fVerboseLevel = VerboseLevel;}
      void SetIsAnalysis(bool val=true){fIsAnalysis=val;};
      void SetIsSimulation(bool val=true){fIsSimulation=val;}

   public: // user definition
      bool HasDefinition(std::string def) {return(fDefinition.find(def)!=fDefinition.end());}

   private:
      // default values
      std::string fDefaultReactionFileName;
      std::string fDefaultDetectorFileName;
      std::string fDefaultRunToReadFileName;
      std::string fDefaultCalibrationFileName;
      std::string fDefaultOutputFileName;
      std::string fDefaultOutputTreeName;
      std::string fDefaultG4MacroPath;
      
      // assigned values
      std::string fReactionFileName;
      std::string fDetectorFileName;
      std::string fRunToReadFileName;
      std::string fCalibrationFileName;
      std::string fOutputFileName;
      std::string fOutputTreeName;
      bool   fIsSplit; // One tree per detector
      bool   fDisableAllBranchOption;
      bool   fInputPhysicalTreeOption;
      bool   fGenerateHistoOption;
      bool   fCheckHistoOption;
      bool   fOnline; // true if spectra server is started
      bool   fLastSimFile;
      bool   fLastPhyFile;
      bool   fLastResFile;
      bool   fLastAnyFile;
      bool   fCircularTree;
      bool   fIsAnalysis;
      bool   fIsSimulation;
      int    fVerboseLevel; // 0 for not talk, 1 for talking
      int    fNumberOfEntryToAnalyse; // use to limit the number of analysed in NPA
      int    fFirstEntryToAnalyse; // use to set the first event analysed in NPA (total: fFirstEntryToAnalyse -> fFirstEntryToAnalyse + fNumberOfEntryToAnalyse)
      int    fSpectraServerPort;
      int    fRandomSeed;
      int    fRecordTrack;
      std::string fSharedLibExtension; // lib extension is platform dependent
      std::string fG4MacroPath; // Path to a geant4 macro to execute at start of nps
      bool fG4BatchMode; // Execute geant4 in batch mode, running the given macro
      std::set<std::string> fDefinition; // a set of user defined definition 
      std::string m_AnalysisOutputPath;// output path of analysed tree
      std::string m_SimulationOutputPath;// output path of simulated tree
      std::string m_EnergyLossPath;// input/output path of energy loss table

};

#endif
