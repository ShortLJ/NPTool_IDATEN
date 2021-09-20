/*****************************************************************************
 * Copyright (C) 2009-2017   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 07/01/11                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: This script loads automatically the NPLib include path and    *
 *             shared libraries.                                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: This script should be called in your .rootlogon.C file           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include<string>

// ROOT headers
#include"TSystem.h"
#include"TROOT.h"
#include"TList.h"
#include"TSystemDirectory.h"

bool contains(std::string path,std::string search){
  if(path.find(search)!=std::string::npos) 
    return true;
  else
    return false;
}

/////////////////////////////////////////////////////
void NPToolLogon(){ 

#ifdef __APPLE__
  std::string Lib_Extension = ".dylib";
#endif
#ifdef __linux__
  std::string Lib_Extension = ".so";
#endif
#ifdef __FreeBSD__
  std::string Lib_Extension = ".so";
#endif

  // Create the NPTool Stype
  std::string NPLPath = gSystem->Getenv("NPTOOL");  
  gROOT->ProcessLine(Form(".x %s/NPLib/scripts/Style_nptool.C",NPLPath.c_str()));
  gROOT->ProcessLine(Form(".x %s/NPLib/scripts/Style_nponline.C",NPLPath.c_str()));

  // Change the standard random generator to TRandom2
  //gRandom = new TRandom2();

  std::string currentpath = gSystem->Getenv("PWD");
  std::string path = gSystem->Getenv("NPTOOL");

  // Add include path
  gROOT->ProcessLine(Form(".include %s/NPLib/include", path.c_str()));

  // Test if the root map exist, 
  // if yes exit
  // if no load the nptool lib
  std::string command = "ls "+path+"/NPLib/lib/*.rootmap > /dev/null 2> /dev/null";
  int return_value = system(command.c_str());
  bool check = false;

  if(return_value==0)
    check=true;

  if(!check){
    // Add shared libraries
    std::string libpath = Form("%s/NPLib/lib", path.c_str());
    TSystemDirectory libdir("libdir", libpath.c_str());
    TList* listfile = libdir.GetListOfFiles();

    // Since the list is ordered alphabetically and that the 
    // libVDetector.dylib library should be loaded before the 
    // lib*Physics.dylib libraries, it is then loaded manually 
    // first.
    // Test if the lib directory is empty or not
    std::string load_path = libpath+"/libNPCore"+Lib_Extension;
    if (listfile->GetEntries() > 2) gSystem->Load(load_path.c_str());

    gSystem->Load("libPhysics.so"); // Needed by Must2 and Sharc
    gSystem->Load("libHist.so"); // Needed by TSpectra Class
    gSystem->Load("libCore.so"); // Need by Maya
    // Loop on Data libraries
    Int_t i = 0;
    while (listfile->At(i)) {
      std::string libname = listfile->At(i++)->GetName();
      if (contains(libname,Lib_Extension) && contains(libname,"Data") && !contains(libname,"libVDetector"+Lib_Extension)) {
        std::string lib = libpath + "/" + libname;
        gSystem->Load(lib.c_str());
      }
    }

    // Loop on Physics Library
    i = 0;
    while (listfile->At(i)) {
      std::string libname = listfile->At(i++)->GetName();
      if (contains(libname,Lib_Extension) && contains(libname,"Physics") &&!contains(libname,"libVDetector"+Lib_Extension)) {
        std::string lib = libpath + "/" + libname;
        gSystem->Load(lib.c_str());
      }
    }

    // Loop on the Reset of the Library
    i = 0;
    while (listfile->At(i)) {
      std::string libname = listfile->At(i++)->GetName();
      if (contains(libname,Lib_Extension) && !contains(libname,"Physics") && !contains(libname,"Data")  &&!contains(libname,"libVDetector"+Lib_Extension)) {
        std::string lib = libpath + "/" + libname;
        gSystem->Load(lib.c_str());
      }
    }

    // gROOT->ProcessLine(".L $NPTOOL/NPLib/include/RootInput.h+");   
    // Since the libdir.GetListOfFiles() commands cds to the
    // libidr directory, one has to return to the initial
    // directory
    //gSystem->cd(currentpath.c_str());
  }
}


