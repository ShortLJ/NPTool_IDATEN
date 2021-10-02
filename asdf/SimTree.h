//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Oct  2 11:06:24 2021 by ROOT version 6.24/02
// from TTree SimulatedTree/Data created / analysed with the nptool package
// found on file: ../Outputs/Simulation/IDATEN07.root
//////////////////////////////////////////////////////////

#ifndef SimTree_h
#define SimTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TInteractionCoordinates.h"
#include "TFatimaData.h"
#include "TKhalaData.h"
#include "TTigressData.h"
#include "TInitialConditions.h"
#include "TTrackInfo.h"

class SimTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   TInteractionCoordinates *InteractionCoordinates;
   TFatimaData     *Fatima;
   TKhalaData      *Khala;
   TTigressData    *Tigress;
   TInitialConditions *InitialConditions;
   TTrackInfo      *TrackInfo;
   Int_t           Run;

   // List of branches
   TBranch        *b_InteractionCoordinates;   //!
   TBranch        *b_Fatima;   //!
   TBranch        *b_Khala;   //!
   TBranch        *b_Tigress;   //!
   TBranch        *b_InitialConditions;   //!
   TBranch        *b_TrackInfo;   //!
   TBranch        *b_Run;   //!

   SimTree(TTree *tree=0, int ifile=0);
   virtual ~SimTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
//   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef SimTree_cxx
SimTree::SimTree(TTree *tree, int ifile=0) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("/home/jlee/NPTools/test/NPTool_IDATEN/Outputs/Simulation/IDATEN%02d.root", ifile));
      if (!f || !f->IsOpen()) {
         f = new TFile(Form("/home/jlee/NPTools/test/NPTool_IDATEN/Outputs/Simulation/IDATEN%02d.root",ifile));
      }
      f->GetObject("SimulatedTree",tree);

   }
   Init(tree);
}

SimTree::~SimTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SimTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SimTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SimTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   InteractionCoordinates = 0;
   Fatima = 0;
   Khala = 0;
   Tigress = 0;
   InitialConditions = 0;
   TrackInfo = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("InteractionCoordinates", &InteractionCoordinates, &b_InteractionCoordinates);
   fChain->SetBranchAddress("Fatima", &Fatima, &b_Fatima);
   fChain->SetBranchAddress("Khala", &Khala, &b_Khala);
   fChain->SetBranchAddress("Tigress", &Tigress, &b_Tigress);
   fChain->SetBranchAddress("InitialConditions", &InitialConditions, &b_InitialConditions);
   fChain->SetBranchAddress("TrackInfo", &TrackInfo, &b_TrackInfo);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   Notify();
}

Bool_t SimTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SimTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SimTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef SimTree_cxx
