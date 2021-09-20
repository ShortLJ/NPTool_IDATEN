#ifndef ROOTINPUT_HH
#define ROOTINPUT_HH
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
 * Last update    : 10/05/21 Add support for Split tree input (A. Matta)     *
 *---------------------------------------------------------------------------*
 * Decription: This class is a singleton class which deals with the ROOT     *
 *             input file and tree both for NPSimulation and NPAnalysis.     *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/



// C++ header
#include <string>
#include <map>
#include <vector>
// ROOT headers
#include "TFile.h"
#include "TChain.h"

// NPL
#include "NPInputParser.h"

class RootInput
{
public:
   // The analysis class is designed to be a singleton (i.e. only one instance
   // can exist). A member function called Instance is defined, which allows
   // the user to get a pointer to the existing instance or to create it if
   // it does not yet exist:
   // (see the constructor for an explanation of the arguments)
   static RootInput* getInstance(std::string configFileName = "configFile");

   // The analysis class instance can be deleted by calling the Destroy
   // method (NOTE: The class destructor is protected, and can thus not be
   // called directly):
   static void Destroy();

protected:
   // Constructor (protected)
   RootInput(std::string configFileName);

   // Destructor (protected)
   virtual ~RootInput();

   // Prevent copying
   RootInput(const RootInput& only);
   const RootInput& operator=(const RootInput& only);

private:
   // The static instance of the RootInput class:
   static RootInput* instance;

public:
   std::string DumpAsciiFile(const char* type, const char* folder = "./.tmp");

public:
   // Return the private chain and file
   TChain*  GetChain()  {return pRootChain;}
   TFile*   GetFile()   {return pRootFile;}
   void     SetChain(TChain* c)  {pRootChain = c;} 

   // Read Input file
   void     ReadOldStyleInputFile(NPL::InputParser& parser);
   void     ReadInputFile(NPL::InputParser& parser);
   void     ReadTreeFile(std::string path);

private:
   TChain   *pRootChain;
   TFile    *pRootFile;
   std::string pTreeName;// the main tree name
   std::vector<std::string> pTreePath;// the main tree path
   // Used for user made tree friends
   std::multimap<std::string,std::vector<std::string>> pFriendsPath;// list of Friends tree indexed by their tree name
   // Used for .tree file input
   std::map<unsigned int, std::vector<std::string>> pFriendsTreePath;// list of Friends tree path indexed by number
   int NumberOfFriend;
  
};

#endif // ROOTINPUT_HH
