/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta   contact address: matta@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// NPL
#include "NPOnlineGUI.h"
#include "NPOptionManager.h"
#include "NPInputParser.h"
#include "NPCore.h"
// STL
#include <iostream>
#include <dirent.h>

// Root
#include "TROOT.h"
#include "TColor.h"
#include "TSystem.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMessage.h"
#include "TGSplitter.h"
#include "TG3DLine.h"
#include "TGDoubleSlider.h"
#include "TGTextEdit.h"
#include "TGLabel.h"
#include "TGComboBox.h"
#include "TASImage.h"
#include "TH2.h"
#include "NPCore.h"
ClassImp(NPL::OnlineGUI);
////////////////////////////////////////////////////////////////////////////////
void NPL::ExecuteMacro(const std::string& name){
  static DIR *dir;
  static struct dirent *ent;
  static std::string path; 
  path = "./online_macros/";
  static std::string filename;
  filename = name+".cxx";
  if ((dir = opendir (path.c_str())) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      if(ent->d_name==filename)
        gROOT->ProcessLine(Form(".x online_macros/%s",filename.c_str()));
    }
    closedir (dir);
  }
}
////////////////////////////////////////////////////////////////////////////////
NPL::OnlineGUI::OnlineGUI(NPL::SpectraClient* client):TGMainFrame(gClient->GetRoot(),1500,500,kMainFrame | kVerticalFrame){
  NPOptionManager::getInstance()->SetVerboseLevel(0);

  m_Client = client;
  m_Sock = 0;
  TString NPLPath = gSystem->Getenv("NPTOOL");
  gROOT->ProcessLine(Form(".x %s/NPLib/scripts/NPToolLogon.C+", NPLPath.Data()));
  gROOT->SetStyle("nponline");

  // Check Elog config 
  m_Elog.ReadConfiguration("elog.txt");

  // Build the interface
  MakeGui();

 
  // Link the button slot to the function
  m_Quit->SetCommand("gApplication->Terminate()");
  m_Connect->Connect("Clicked()", "NPL::OnlineGUI", this, "Connect()");
  m_Update->Connect("Clicked()", "NPL::OnlineGUI", this, "Update()");
  m_Clock->Connect("Clicked()","NPL::OnlineGUI",this,"AutoUpdate()");
  m_FitAll->Connect("Clicked()", "NPL::OnlineGUI", this, "FitAll()");
  m_FitCurrent->Connect("Clicked()", "NPL::OnlineGUI", this, "FitCurrent()");
  m_ResetCurrent->Connect("Clicked()","NPL::OnlineGUI",this,"ResetCurrent()");
  m_ResetAll->Connect("Clicked()","NPL::OnlineGUI",this,"ResetAll()");
  m_ApplyRangeCurrent->Connect("Clicked()","NPL::OnlineGUI",this,"ApplyRangeCurrent()");
  m_ApplyRangeAll->Connect("Clicked()","NPL::OnlineGUI",this,"ApplyRangeAll()");
  m_SaveAs->Connect("Clicked()","NPL::OnlineGUI",this,"SaveAs()");
  m_Eloging->Connect("Clicked()","NPL::OnlineGUI",this,"Eloging()");

  // Connect to server 
  Connect();



}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::SaveAs(){
  TCanvas* c = m_EmbeddedCanvas->GetCanvas();
 
  if(!c)
   return;

  if(m_SaveAsColor->IsOn()){
   gROOT->SetStyle("Modern");
   c->UseCurrentStyle();
   c->SaveAs(m_SaveAsFileName->GetText());
   gROOT->SetStyle("nponline");
   c->UseCurrentStyle();
  }

  else
    c->SaveAs(m_SaveAsFileName->GetText());

  c->Update();
}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::ResetAll(){
  TCanvas* c = m_EmbeddedCanvas->GetCanvas();
  if (!c)
    return;

  int size= ((TList*)c->GetListOfPrimitives())->GetSize();
  for(unsigned int i =  1 ; i < size ;i++){
    c->cd(i);
    ResetCurrent();
  }
  c->cd(1);
}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::ResetCurrent(){
  TCanvas* c = m_EmbeddedCanvas->GetCanvas();
  if (!c)
    return;

  // reset log scale attribute
  gPad->SetLogx(false);
  gPad->SetLogy(false);
  gPad->SetLogz(false);

  // loop on histograms, reset axis and content
  TList* list = gPad->GetListOfPrimitives();
  int Hsize = list->GetSize();
  for(int h = 0 ; h < Hsize ; h++){
     TObject* obj = list->At(h);
     if(obj->InheritsFrom(TH1::Class())){
        TH1* h = (TH1*) obj;
        h->GetXaxis()->UnZoom();
        h->GetYaxis()->UnZoom();
        h->GetZaxis()->UnZoom();
        h->Reset("ICESM");
     }
  }
  c->Update();
}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::Eloging(){

  std::vector<std::string> attributes;
  std::vector<std::string> val;

  std::map<std::string,TGTextEntry*>::iterator it ;
  for(it = m_ElogAttributes.begin(); it != m_ElogAttributes.end() ; it++){
    attributes.push_back(it->first);
    val.push_back(it->second->GetText());
  }

  // Create the file to be attached
  std::vector<std::string> attachement; 

  if(m_EmbeddedCanvas->GetCanvas()){
    m_EmbeddedCanvas->GetCanvas()->SaveAs("elog.pdf");
    attachement.push_back("elog.pdf");
  }


  m_Elog.CreateEntry(attributes,val,m_ElogEntry->GetText()->AsString().Data(),attachement);


}

////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::ApplyRangeAll(){
  TCanvas* c = m_EmbeddedCanvas->GetCanvas();
  if (!c)
    return;

  int size= ((TList*)c->GetListOfPrimitives())->GetSize();
  for(unsigned int i =  1 ; i < size ;i++){
    c->cd(i);
    ApplyRangeCurrent();
  }
  c->cd(1);
}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::ApplyRangeCurrent(){
  TCanvas* c = m_EmbeddedCanvas->GetCanvas();
  if (!c)
    return;

  // Log Scale
  if(m_CheckLogX->IsOn())
    gPad->SetLogx();
  else
    gPad->SetLogx(false);
  if(m_CheckLogY->IsOn())
    gPad->SetLogy();
  else
    gPad->SetLogy(false);
  if(m_CheckLogZ->IsOn())
    gPad->SetLogz();
  else
    gPad->SetLogz(false);

  if(m_Xmin->GetNumber() != m_Xmax->GetNumber()){
    TList* list = gPad->GetListOfPrimitives();
    int Hsize = list->GetSize();
    for(int h = 0 ; h < Hsize ; h++){
      TObject* obj = list->At(h);
      if(obj->InheritsFrom(TH1::Class())){
        TH1* h = (TH1*) obj;
        h->GetXaxis()->SetRangeUser(m_Xmin->GetNumber(),m_Xmax->GetNumber());
      }
    }
  }

  if(m_Ymin->GetNumber() != m_Ymax->GetNumber()){
    TList* list = gPad->GetListOfPrimitives();
    int Hsize = list->GetSize();
    for(int h = 0 ; h < Hsize ; h++){
      TObject* obj = list->At(h);
      if(obj->InheritsFrom(TH1::Class())){
        TH1* h = (TH1*) obj;
        h->GetYaxis()->SetRangeUser(m_Ymin->GetNumber(),m_Ymax->GetNumber());
      }
    }
  }

  gPad->Update(); 
}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::FitCurrent(){
  static std::string gauss_formula = "([0]*[3]/([2]*sqrt(2*3.14159265359)))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))";

  static std::string full_formula = gauss_formula +"+[P0]+[P1]*x";
  TList* list = gPad->GetListOfPrimitives();
    int Hsize = list->GetSize();
    for(int h = 0 ; h < Hsize ; h++){
      TObject* obj = list->At(h);
      if(obj->InheritsFrom(TH1::Class()) && !obj->InheritsFrom(TH2::Class())){
        TF1* fit1 = new TF1("Gaussian",gauss_formula.c_str(),-1000,1000);
        TF1* fit2 = new TF1("Gaussian+Background",full_formula.c_str(),-1000,1000);
        fit1->SetParNames("Integral","mean","#sigma","binning");
        fit2->SetParNames("Integral","mean","#sigma","binning","base line","slope");
        fit1->SetNpx(1000);
        fit2->SetNpx(1000);
        TH1* hh = (TH1*)obj;

        if(m_BackgroundFit->IsOn()){
          fit2->SetParameter(0,hh->GetMaximum()*hh->GetBinWidth(hh->GetMaximumBin()));
          fit2->SetParameter(1,hh->GetBinCenter(hh->GetMaximumBin()));
          fit2->SetParameter(2,10*hh->GetBinWidth(hh->GetMaximumBin()));
          fit2->FixParameter(3,hh->GetBinWidth(hh->GetMaximumBin()));
          fit2->SetParameter(4,hh->GetMinimum());
          fit2->SetParameter(5,0);
          hh->Fit(fit2,"IQ");
        }
        else{
          fit1->SetParameter(0,hh->GetMaximum()/hh->GetBinWidth(hh->GetMaximumBin()));
          fit1->SetParameter(1,hh->GetBinCenter(hh->GetMaximumBin()));
          fit1->SetParameter(2,10*hh->GetBinWidth(hh->GetMaximumBin()));
          fit1->FixParameter(3,hh->GetBinWidth(hh->GetMaximumBin()));
          hh->Fit(fit1,"IQ"); 
          
        }
      }
    }
  gPad ->Update();

}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::FitAll(){
  TCanvas* c = m_EmbeddedCanvas->GetCanvas();

    int size= ((TList*)c->GetListOfPrimitives())->GetSize();
    for(unsigned int i =  1 ; i < size ;i++){
      c->cd(i);
      FitCurrent(); 
    }
     
   c->cd(1);
}
////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::MakeGui(){
  m_BgColor = gROOT->GetColor(kGray+3)->GetPixel();
  m_FgColor = gROOT->GetColor(kAzure+7)->GetPixel();
  m_TabBgColor = gROOT->GetColor(kGray+3)->GetPixel();
  m_TabFgColor = gROOT->GetColor(kAzure+7)->GetPixel();
  m_Timer = 0;

  // main frame
  m_Main = this;
  m_Main->SetName("nponline");
  m_Main->SetBackgroundColor(m_BgColor);
  m_Main->SetForegroundColor(m_FgColor);
  m_Main->SetWindowName("nponline");
  // Button bar to hold the button
  m_ButtonBar= new TGVerticalFrame(m_Main,10000,42,kFixedSize);
  m_ButtonBar->SetBackgroundColor(m_BgColor);
  m_ButtonBar->SetForegroundColor(m_BgColor);
  m_ButtonBar->SetLayoutBroken(kTRUE);
  m_Main->AddFrame(m_ButtonBar,new TGLayoutHints(kLHintsLeft|kLHintsTop));

  std::string NPLPath = gSystem->Getenv("NPTOOL");  
  std::string path_quit = NPLPath+"/NPLib/Core/icons/power.xpm";
  m_Quit = new TGPictureButton(m_ButtonBar,gClient->GetPicture(path_quit.c_str()),-1,TGPictureButton::GetDefaultGC()(),kChildFrame);
  m_Quit->SetBackgroundColor(m_BgColor);
  m_Quit->SetToolTipText("Quit");

  m_ButtonBar->AddFrame(m_Quit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  m_Quit->MoveResize(10,5,32,32);

  std::string path_connect = NPLPath+"/NPLib/Core/icons/plugin.xpm";
  m_Connect = new TGPictureButton(m_ButtonBar,gClient->GetPicture(path_connect.c_str()),-1,TGPictureButton::GetDefaultGC()(),kChildFrame);
  std::string path_connected = NPLPath+"/NPLib/Core/icons/brightness.xpm"; 
  m_Connect->SetDisabledPicture(gClient->GetPicture(path_connected.c_str()));
  m_Connect->SetBackgroundColor(m_BgColor);
  m_Connect->SetForegroundColor(m_BgColor);

  m_Connect->SetToolTipText("Connect to server");
  m_ButtonBar->AddFrame(m_Connect, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  m_Connect->MoveResize(52,5,32,32);

  std::string path_update = NPLPath+"/NPLib/Core/icons/download.xpm";
  m_Update = new TGPictureButton(m_ButtonBar,gClient->GetPicture(path_update.c_str()),-1,TGPictureButton::GetDefaultGC()(),kChildFrame);
  m_Update->SetBackgroundColor(m_BgColor);
  m_Update->SetForegroundColor(m_BgColor);
  m_Update->SetToolTipText("Update spectra");
  m_ButtonBar->AddFrame(m_Update, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  m_Update->MoveResize(400,5,32,32);

  std::string path_clock= NPLPath+"/NPLib/Core/icons/clock.xpm";
  m_Clock = new TGPictureButton(m_ButtonBar,gClient->GetPicture(path_clock.c_str()),-1,TGPictureButton::GetDefaultGC()(),kChildFrame);
  m_Clock->SetBackgroundColor(m_BgColor);
  m_Clock->SetForegroundColor(m_BgColor);

  m_Clock->SetToolTipText("AutoUpdate");
  m_ButtonBar->AddFrame(m_Clock, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  m_Clock->MoveResize(442,5,32,32);

  TGFont* ufont;// will reflect user font changes
  ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

  TGGC  * uGC;// will reflect user GC changes
  // graphics context changes
  GCValues_t valress;
  valress.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
  gClient->GetColorByName("#000000",valress.fForeground);
  gClient->GetColorByName("#e7e7e7",valress.fBackground);
  valress.fFillStyle = kFillSolid;
  valress.fFont = ufont->GetFontHandle();
  valress.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&valress, kTRUE);
  m_Address = new TGTextEntry(m_ButtonBar, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kChildFrame | kOwnBackground);
  m_Address->SetMaxLength(4096);
  m_Address->SetAlignment(kTextLeft);
  m_Address->SetText(m_Client->GetAddress().c_str());
  m_Address->Resize(200,m_Address->GetDefaultHeight());
  m_ButtonBar->AddFrame(m_Address, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  m_Address->MoveResize(90,10,200,20);

  m_Port = new TGNumberEntry(m_ButtonBar, (Double_t) m_Client->GetPort(),9,-1,(TGNumberFormat::EStyle) 5);
  m_Port->SetName("m_Port");
  m_Port->GetButtonUp()->SetStyle(1);
  m_Port->GetButtonDown()->SetStyle(1);
  m_ButtonBar->AddFrame(m_Port, new TGLayoutHints(kLHintsLeft));
  m_Port->MoveResize(300,10,80,20);

  m_TimerEntry = new TGNumberEntry(m_ButtonBar, (Double_t) 1,9,-1,(TGNumberFormat::EStyle) 5);
  m_TimerEntry->SetName("m_TimerEntry");
  m_TimerEntry->GetButtonUp()->SetStyle(1);
  m_TimerEntry->GetButtonDown()->SetStyle(1);
  m_ButtonBar->AddFrame(m_TimerEntry, new TGLayoutHints(kLHintsLeft));
  m_TimerEntry->MoveResize(484,10,40,20);

  // Create the splitted frame
  m_Split = new TGHorizontalFrame(m_Main, 50, 50);
  TGVerticalFrame* fV1 = new TGVerticalFrame(m_Split, 10, 10, kFixedWidth);
  TGVerticalFrame* fV2 = new TGVerticalFrame(m_Split, 10, 10);
  TGVerticalFrame* fV3 = new TGVerticalFrame(m_Split, 10, 10, kFixedWidth);

  m_Left   = new TGCompositeFrame(fV1, 10, 10, kChildFrame);
  m_Center = new TGCompositeFrame(fV2, 10, 10, kChildFrame);
  m_Right  = new TGCompositeFrame(fV3, 10, 10, kChildFrame); 
  fV1->AddFrame(m_Left,  new TGLayoutHints(kLHintsLeft |    kLHintsExpandX | kLHintsExpandY,0, 0, 5, 10));
  fV2->AddFrame(m_Center,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,0, 0, 5, 10));
  fV3->AddFrame(m_Right, new TGLayoutHints(kLHintsRight |   kLHintsExpandX | kLHintsExpandY,0, 0, 5, 10));

  fV1->Resize(200, 400);
  fV2->Resize(400, 400);
  fV3->Resize(200, 400);


  // Change the size of the list tree
  m_Split->AddFrame(fV1, new TGLayoutHints(kLHintsLeft | kLHintsExpandY ));

  TGVSplitter* splitter1 = new TGVSplitter(m_Split,5,5);
  splitter1->SetFrame(fV1, kTRUE);
  m_Split->AddFrame(splitter1, new TGLayoutHints(kLHintsLeft| kLHintsTop |  kLHintsExpandY));
  m_Split->AddFrame(fV2,  new TGLayoutHints(kLHintsExpandX| kLHintsTop |  kLHintsExpandY));

  // Change the size of the tool bar
  TGVSplitter* splitter2 = new TGVSplitter(m_Split,5,5);
  splitter2->SetFrame(fV3, false);
  m_Split->AddFrame(splitter2, new TGLayoutHints(kLHintsLeft| kLHintsTop |  kLHintsExpandY));
  m_Split->AddFrame(fV3, new TGLayoutHints(kLHintsRight| kLHintsTop |  kLHintsExpandY));

  splitter1->SetBackgroundColor(m_BgColor);    
  splitter1->SetForegroundColor(m_BgColor);    

  splitter2->SetBackgroundColor(m_BgColor);    
  splitter2->SetForegroundColor(m_BgColor);    


  m_Split->SetBackgroundColor(m_BgColor); 
  m_Split->SetForegroundColor(m_BgColor); 

  m_Left->SetBackgroundColor(m_BgColor); 
  m_Right->SetBackgroundColor(m_FgColor);    
  m_Left->SetForegroundColor(m_BgColor); 
  m_Right->SetForegroundColor(m_BgColor);    
  fV1->SetBackgroundColor(m_BgColor); 
  fV2->SetBackgroundColor(m_BgColor);    
  fV1 ->SetForegroundColor(m_BgColor); 
  fV2->SetForegroundColor(m_BgColor);    
  fV3->SetBackgroundColor(m_BgColor);
  m_Main->AddFrame(m_Split, new TGLayoutHints(kLHintsRight | kLHintsExpandX |kLHintsExpandY));

  // Right
  // Create a fit tool bar
  //m_Right->SetLayoutBroken(kTRUE);

  // Navigation
  TGVerticalFrame* m_NavBar= new TGVerticalFrame(m_Right,10000,80);
  m_NavBar->SetBackgroundColor(m_FgColor);
  m_NavBar->SetForegroundColor(m_FgColor);
  m_Right->AddFrame(m_NavBar, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  TGVerticalFrame* m_LogBar= new TGVerticalFrame(m_NavBar,10000,80);
  m_LogBar->SetBackgroundColor(m_FgColor);
  m_LogBar->SetForegroundColor(m_FgColor);
  m_LogBar->SetLayoutManager(new TGHorizontalLayout(m_LogBar));
  m_NavBar->AddFrame(m_LogBar, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  m_CheckLogX = new TGCheckButton(m_LogBar,"Log&X");
  m_CheckLogX->SetBackgroundColor(m_FgColor);
  m_CheckLogX->SetForegroundColor(m_BgColor);        
  m_LogBar->AddFrame(m_CheckLogX, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 
  m_CheckLogX->Move(0,0); 

  m_CheckLogY = new TGCheckButton(m_LogBar,"Log&Y");
  m_CheckLogY->SetBackgroundColor(m_FgColor);
  m_CheckLogY->SetForegroundColor(m_BgColor);        
  m_LogBar->AddFrame(m_CheckLogY, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 
  m_CheckLogY->Move(60,0); 

  m_CheckLogZ = new TGCheckButton(m_LogBar,"Log&Z");
  m_CheckLogZ->SetBackgroundColor(m_FgColor);
  m_CheckLogZ->SetForegroundColor(m_BgColor);        
  m_LogBar->AddFrame(m_CheckLogZ, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 
  m_CheckLogZ->Move(120 ,0); 


  TGVerticalFrame* m_XRange= new TGVerticalFrame(m_NavBar,10000,80);
  m_XRange->SetBackgroundColor(m_FgColor);
  m_XRange->SetForegroundColor(m_FgColor);
  m_XRange->SetLayoutManager(new TGHorizontalLayout(m_XRange));
  m_NavBar->AddFrame(m_XRange, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  TGLabel* Xlabel = new TGLabel(m_XRange, "X");
  Xlabel->SetBackgroundColor(m_FgColor);
  m_XRange->AddFrame(Xlabel, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 

  m_Xmin = new TGNumberEntry(m_XRange,0,6);
  m_XRange->AddFrame(m_Xmin, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 
  m_Xmax = new TGNumberEntry(m_XRange,0,6);
  m_XRange->AddFrame(m_Xmax, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 

  TGVerticalFrame* m_YRange= new TGVerticalFrame(m_NavBar,10000,80);
  m_YRange->SetBackgroundColor(m_FgColor);
  m_YRange->SetForegroundColor(m_FgColor);
  m_YRange->SetLayoutManager(new TGHorizontalLayout(m_YRange));
  m_NavBar->AddFrame(m_YRange, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  TGLabel* Ylabel = new TGLabel(m_YRange, "Y");
  Ylabel->SetBackgroundColor(m_FgColor);
  m_YRange->AddFrame(Ylabel, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 

  m_Ymin = new TGNumberEntry(m_YRange,0,6);
  m_YRange->AddFrame(m_Ymin, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 
  m_Ymax = new TGNumberEntry(m_YRange,0,6);
  m_YRange->AddFrame(m_Ymax, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 

  TGVerticalFrame* m_RangeButton= new TGVerticalFrame(m_NavBar,10000,80);
  m_RangeButton->SetBackgroundColor(m_FgColor);
  m_RangeButton->SetForegroundColor(m_FgColor);
  m_RangeButton->SetLayoutManager(new TGHorizontalLayout(m_RangeButton));
  m_NavBar->AddFrame(m_RangeButton, new TGLayoutHints(kLHintsTop|kLHintsExpandX));


  m_ApplyRangeCurrent= new  TGTextButton(m_RangeButton, "C&urrent Pad",-1);
  m_RangeButton->AddFrame(m_ApplyRangeCurrent, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  m_ApplyRangeAll= new  TGTextButton(m_RangeButton, "All &Pad",-1);
  m_RangeButton->AddFrame(m_ApplyRangeAll, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  TGVerticalFrame* m_ResetButton= new TGVerticalFrame(m_NavBar,10000,80);
  m_ResetButton->SetBackgroundColor(m_FgColor);
  m_ResetButton->SetForegroundColor(m_FgColor);
  m_ResetButton->SetLayoutManager(new TGHorizontalLayout(m_ResetButton));
  m_NavBar->AddFrame(m_ResetButton, new TGLayoutHints(kLHintsTop|kLHintsExpandX));


  m_ResetCurrent= new  TGTextButton(m_ResetButton, "Reset Current",-1);
  m_ResetButton->AddFrame(m_ResetCurrent, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  m_ResetAll= new  TGTextButton(m_ResetButton, "Reset All",-1);
  m_ResetButton->AddFrame(m_ResetAll, new TGLayoutHints(kLHintsTop|kLHintsExpandX));



  // Horizontal separator
  TGHorizontal3DLine* m_NavLine = new TGHorizontal3DLine(m_Right,408,8);
  m_NavLine->SetBackgroundColor(m_FgColor);
  m_NavLine->SetForegroundColor(m_BgColor);
  m_Right->AddFrame(m_NavLine, new TGLayoutHints(kLHintsExpandX,2,2,2,2));  


  // Fit Bar
  TGVerticalFrame* m_FitBar = new TGVerticalFrame(m_Right,10000,80);
  m_FitBar->SetBackgroundColor(m_FgColor);
  m_FitBar->SetForegroundColor(m_FgColor);
  m_FitBar->SetLayoutManager(new TGVerticalLayout(m_FitBar));  
  m_Right->AddFrame(m_FitBar, new TGLayoutHints(kLHintsTop|kLHintsExpandX));


  // Check box container
  TGVerticalFrame* m_SelectBar = new TGVerticalFrame(m_FitBar,10000,80);
  m_SelectBar->SetBackgroundColor(m_FgColor);
  m_SelectBar->SetForegroundColor(m_FgColor);
  m_SelectBar->SetLayoutManager(new TGHorizontalLayout(m_SelectBar));
  m_FitBar->AddFrame(m_SelectBar, new TGLayoutHints(kLHintsLeft|kLHintsExpandX));

  // Fit check box
  m_BackgroundFit = new TGCheckButton(m_SelectBar,"&Background");
  m_BackgroundFit->SetBackgroundColor(m_FgColor);
  m_BackgroundFit->SetForegroundColor(m_BgColor);        
  m_SelectBar->AddFrame(m_BackgroundFit, new TGLayoutHints(kLHintsTop|kLHintsLeft,2,2,2,2)); 

  // Fit button
  TGVerticalFrame* m_FitButton= new TGVerticalFrame(m_FitBar,10000,80);
  m_FitButton->SetBackgroundColor(m_FgColor);
  m_FitButton->SetForegroundColor(m_FgColor);
  m_FitButton->SetLayoutManager(new TGHorizontalLayout(m_FitButton));
  m_FitBar->AddFrame(m_FitButton, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  m_FitCurrent= new  TGTextButton(m_FitButton, "Fit &Current",-1);
  m_FitButton->AddFrame(m_FitCurrent, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  m_FitAll= new  TGTextButton(m_FitButton, "Fit &All",-1);
  m_FitButton->AddFrame(m_FitAll, new TGLayoutHints(kLHintsTop|kLHintsExpandX));


  // Horizontal separator
  TGHorizontal3DLine* m_FitLine = new TGHorizontal3DLine(m_Right,408,8);
  m_FitLine->SetBackgroundColor(m_FgColor);
  m_FitLine->SetForegroundColor(m_BgColor);
  m_Right->AddFrame(m_FitLine, new TGLayoutHints(kLHintsExpandX,2,2,2,2));  

  // SaveAs Bar
  TGVerticalFrame* m_SaveAsBar = new TGVerticalFrame(m_Right,10000,80);
  m_SaveAsBar->SetBackgroundColor(m_FgColor);
  m_SaveAsBar->SetForegroundColor(m_FgColor);
  m_SaveAsBar->SetLayoutManager(new TGHorizontalLayout(m_SaveAsBar));  
  m_Right->AddFrame(m_SaveAsBar, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

  // SaveAs button
  std::string path_print= NPLPath + "/NPLib/Core/icons/print.xpm";
  m_SaveAs= new TGPictureButton(m_SaveAsBar,gClient->GetPicture(path_print.c_str()),-1,TGPictureButton::GetDefaultGC()(),kChildFrame);
  m_SaveAs->SetBackgroundColor(m_FgColor);
  m_SaveAs->SetToolTipText("SaveAs");
  m_SaveAsBar->AddFrame(m_SaveAs, new TGLayoutHints(kLHintsTop|kLHintsLeft,10,10,10,10));

  // Check box container
  TGVerticalFrame* m_OptionBar = new TGVerticalFrame(m_SaveAsBar,10000,80);
  m_OptionBar->SetBackgroundColor(m_FgColor);
  m_OptionBar->SetForegroundColor(m_FgColor);
  m_OptionBar->SetLayoutManager(new TGVerticalLayout(m_OptionBar));
  m_SaveAsBar->AddFrame(m_OptionBar, new TGLayoutHints(kLHintsLeft|kLHintsExpandX));

  // SaveAs color
  m_SaveAsColor = new TGCheckButton(m_OptionBar,"printer color");
  m_SaveAsColor->SetBackgroundColor(m_FgColor);
  m_SaveAsColor->SetForegroundColor(m_BgColor);        
  m_OptionBar->AddFrame(m_SaveAsColor, new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY,2,2,2,2)); 

  // Format
  m_SaveAsFileName = new TGTextEntry(m_OptionBar,"nponline.pdf");
  m_OptionBar->AddFrame(m_SaveAsFileName, new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY,1,1,1,1)); 


  // Horizontal separator
  TGHorizontal3DLine* m_SaveAsLine = new TGHorizontal3DLine(m_Right,408,8);
  m_SaveAsLine->SetBackgroundColor(m_FgColor);
  m_SaveAsLine->SetForegroundColor(m_BgColor);
  m_Right->AddFrame(m_SaveAsLine, new TGLayoutHints(kLHintsExpandX,2,2,2,2));  

  // Elog button
  std::string path_elog= NPLPath +  "/NPLib/Core/icons/booklet.xpm";
  m_Eloging= new TGPictureButton(m_Right,gClient->GetPicture(path_elog.c_str()),-1,TGPictureButton::GetDefaultGC()(),kChildFrame);
  m_Eloging->SetBackgroundColor(m_FgColor);
  m_Eloging->SetToolTipText("Elog");
  m_Right->AddFrame(m_Eloging, new TGLayoutHints(kLHintsTop|kLHintsLeft,10,10,10,10));

  // Elog attributes menu
  std::map<std::string , std::vector <std::string> > attributes  = m_Elog.GetAttributesValues();
  std::map<std::string , std::vector <std::string> >::iterator it;
  for(it = attributes.begin() ; it != attributes.end() ; it++){
    TGVerticalFrame* attframe= new TGVerticalFrame(m_Right,10000,80);
    attframe->SetBackgroundColor(m_FgColor);
    attframe->SetForegroundColor(m_FgColor);
    attframe->SetLayoutManager(new TGHorizontalLayout(attframe));
    m_Right->AddFrame(attframe, new TGLayoutHints(kLHintsTop|kLHintsExpandX));

    TGLabel* attlabel = new TGLabel(attframe, it->first.c_str());
    attlabel->SetBackgroundColor(m_FgColor);
    attframe->AddFrame(attlabel, new TGLayoutHints(kLHintsCenterX,1,1,1,1)); 

   // free type case 
   if(it->second.size()==1){
    TGTextEntry* attentry= new TGTextEntry(attframe,it->second[0].c_str());  
    attframe->AddFrame(attentry, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX,1,1,1,1)); 
    m_ElogAttributes[it->first] = attentry;
   }
   // selection menu case
   else{
    TGComboBox* attentry= new TGComboBox(attframe);  
    unsigned int size = it->second.size();
    for(unsigned int i = 0 ; i < size ; i++){
      attentry->AddEntry(it->second[i].c_str(),i);
    }
    attframe->AddFrame(attentry, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX|kLHintsExpandY,1,1,1,1)); 
    attentry->EnableTextInput(true);
    attentry->Select(0);
    m_ElogAttributes[it->first] = attentry->GetTextEntry();
   }
  }

  // Elog entry
  m_ElogEntry = new TGTextEdit(m_Right,m_Right->GetWidth(),100);
  m_Right->AddFrame(m_ElogEntry, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY |kLHintsCenterY, 1, 1, 1, 1));


  // Center //
  m_EmbeddedCanvas = new TRootEmbeddedCanvas("Display",m_Center,10000,10000,!kSunkenFrame);  
  m_EmbeddedCanvas->SetAutoFit(true);
  m_Center->AddFrame(m_EmbeddedCanvas,new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandX | kLHintsExpandY));
  TCanvas* c = NULL;
  m_EmbeddedCanvas->AdoptCanvas(c);

  // Left //
  // Create tabs for histo list and canvas list
  m_Tab = new TGTab(m_Left,700,500);
  m_Tab->Resize(m_Tab->GetDefaultSize());
  m_Tab->SetBackgroundColor(m_TabBgColor);
  m_Tab->SetForegroundColor(m_TabFgColor);
  m_Tab->ChangeSubframesBackground(m_BgColor);
  TGCompositeFrame* tfCanvas= m_Tab->AddTab("Canvas");
  TGCompositeFrame* tfHisto = m_Tab->AddTab("Histo");

  m_Left->AddFrame(m_Tab, new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandY | kLHintsExpandX));

  // canvas widget
  TGCanvas* m_ListCanvas = new TGCanvas(tfCanvas,120,500);
  m_ListCanvas->SetName("m_ListCanvas");

  m_ListCanvas->SetForegroundColor(m_BgColor); 
  m_ListCanvas->SetForegroundColor(m_BgColor);    

  // canvas viewport
  TGViewPort* CanvasViewPort = m_ListCanvas->GetViewPort();

  // list tree
  m_CanvasListTree = new CanvasList(m_Main,m_ListCanvas,m_EmbeddedCanvas,m_Client->GetSpectra());
  m_ListTree = m_CanvasListTree->GetListTree();

  CanvasViewPort->AddFrame(m_ListTree,new TGLayoutHints(kLHintsRight | kLHintsBottom | kLHintsExpandY | kLHintsExpandX));
  m_ListTree->SetLayoutManager(new TGHorizontalLayout(m_ListTree));
  m_ListTree->MapSubwindows();

  m_ListCanvas->SetContainer(m_ListTree);
  m_ListCanvas->MapSubwindows();
  tfCanvas->AddFrame(m_ListCanvas, new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandY | kLHintsExpandX));
  m_ListCanvas->MoveResize(10,50,120,500);

  m_Main->SetMWMHints(kMWMDecorAll,kMWMFuncAll,kMWMInputModeless);
  m_Main->MapSubwindows();

  m_Main->Resize(m_Main->GetDefaultSize());
  m_Main->MapWindow();
  m_Main->MoveResize(50,50,1250,800);

}

////////////////////////////////////////////////////////////////////////////////
NPL::OnlineGUI::~OnlineGUI(){
  delete m_Main; 
  delete m_ListCanvas;
  delete m_ListTree;
  delete m_Tab;
  delete m_Quit;
  delete m_Connect;
  delete m_Update;
  delete m_Sock;
  delete m_Port;
  delete m_Address;
  delete m_StatusBar;
}

////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::Connect(){
  m_Client->SetAddressAndPort((std::string) m_Address->GetDisplayText().Data(),(int) m_Port->GetNumber());
  m_Client->Connect();
  m_CanvasListTree->LoadCanvasList(m_Client->GetSpectra());
}

////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::Update(){
    TCanvas* c = m_EmbeddedCanvas->GetCanvas();
    int current = gPad->GetNumber();
    TList* first = c->GetListOfPrimitives();
    int size= first->GetSize();
    for(unsigned int i =  0 ; i < size ;i++){
      if(first->At(i)->InheritsFrom(TPad::Class())){
      dynamic_cast<TPad*>(first->At(i))->cd();
      TList* list = gPad->GetListOfPrimitives();
      int Hsize = list->GetSize();
      for(int h = 0 ; h < Hsize ; h++){
        TObject* obj = list->At(h);
        if(obj->InheritsFrom(TH1::Class())){
          m_Client->Update(obj->GetName());
          obj->Paint();
          ExecuteMacro(obj->GetName());
          gPad->Update();
        }
      }
     }
    }
    c->cd(current);
    c->Update();

}

////////////////////////////////////////////////////////////////////////////////
void NPL::OnlineGUI::AutoUpdate(){

  if(m_Timer){
    delete m_Timer;
    m_Timer = 0 ;
    return;
  }

  else if(m_TimerEntry->GetNumber()>0){
    m_Timer = new TTimer(m_TimerEntry->GetNumber()*1000);
    m_Timer->Connect("Timeout()", "NPL::OnlineGUI", this, "Update()");
    m_Timer->TurnOn();
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* CanvasList Class */

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

NPL::CanvasList::CanvasList(TGMainFrame* main, TGCanvas* parent,TRootEmbeddedCanvas* canvas,TList* Spectra){
  std::string NPLPath = gSystem->Getenv("NPTOOL");  
  std::string path_icon = NPLPath+"/NPLib/Core/icons/polaroid.xpm";
  std::string path_icon_folder = NPLPath+"/NPLib/Core/icons/folder.xpm";

  m_popen = gClient->GetPicture(path_icon.c_str());
  m_pclose = gClient->GetPicture(path_icon.c_str());
  m_pfolder = gClient->GetPicture(path_icon_folder.c_str()); 

  m_BgColor = gROOT->GetColor(kGray+3)->GetPixel();
  m_FgColor = gROOT->GetColor(kWhite)->GetPixel();

  m_ListTree = new TGListTree(parent,kHorizontalFrame);
  m_ListTree->Connect("DoubleClicked(TGListTreeItem*,Int_t)","NPL::CanvasList",this,"OnDoubleClick(TGListTreeItem*,Int_t)");
  m_Main = main;
  m_EmbeddedCanvas = canvas;
  m_OldCurrent     = "_void_";
}
////////////////////////////////////////////////////////////////////////////////
NPL::CanvasList::~CanvasList(){
}
////////////////////////////////////////////////////////////////////////////////
void NPL::CanvasList::OnDoubleClick(TGListTreeItem* item, Int_t btn){
  TCanvas* c = m_Canvas[item->GetText()];

  if(c){
    m_EmbeddedCanvas->AdoptCanvas(c);
    TGDimension size = m_EmbeddedCanvas->GetContainer()->GetSize();
    m_EmbeddedCanvas->GetContainer()->Resize(0,0);
    m_EmbeddedCanvas->GetContainer()->Resize(size);
    c->SetHighLightColor(kAzure+7); // color of active pad
    c->Draw();
    c->Update();
    c->cd(1);
  }
}
////////////////////////////////////////////////////////////////////////////////
void NPL::CanvasList::AddItem(TCanvas* c,TGListTreeItem* parent){
  TGListTreeItem* item  = m_ListTree->AddItem(parent,c->GetName());
  item->SetPictures(m_popen, m_pclose);
  if(parent)
    parent->SetPictures(m_pfolder,m_pfolder);
  m_Canvas[c->GetName()]=c;
}
////////////////////////////////////////////////////////////////////////////////
void NPL::CanvasList::Clear(){
  m_OldCurrent=m_EmbeddedCanvas->GetCanvas()->GetName();
  
  std::map<std::string,TCanvas*>::iterator it ;
  // delete all old canvas
  for(it=m_Canvas.begin(); it!=m_Canvas.end();it++){
    delete it->second;
    }
  m_Canvas.clear();

  // Clear the list tree
  TGListTreeItem* item =  m_ListTree->GetFirstItem() ;
  while(item){
    m_ListTree->DeleteItem(item);
    item = m_ListTree->GetFirstItem() ;
  }
}
////////////////////////////////////////////////////////////////////////////////
TGListTree* NPL::CanvasList::GetListTree(){
  return m_ListTree;
}
////////////////////////////////////////////////////////////////////////////////
void NPL::CanvasList::SetTab(TGTab* tab){
  m_Tab=tab;
}
////////////////////////////////////////////////////////////////////////////////
void NPL::CanvasList::SetStatusText(const char* txt, int pi){
  for(unsigned int i = 0 ; i < m_StatusBar.size(); i++)
    m_StatusBar[i]->SetText(txt,pi);
}
////////////////////////////////////////////////////////////////////////////////
void NPL::CanvasList::EventInfo(int event,int px,int py,TObject* selected){
  const char *text0, *text1, *text3;
  char text2[50];
  text0 = selected->GetTitle();
  SetStatusText(text0,0);
  text1 = selected->GetName();
  SetStatusText(text1,1);
  if (event == kKeyPress)
    sprintf(text2, "%c", (char) px);
  else
    sprintf(text2, "%d,%d", px, py);
  SetStatusText(text2,2);
  text3 = selected->GetObjectInfo(px,py);
  SetStatusText(text3,3);
}
////////////////////////////////////////////////////////////////////////////////
void NPL::CanvasList::LoadCanvasList(TList* Spectra){
  if(!Spectra)
    return;
  Clear();
  NPL::InputParser parser("CanvasList.txt",false);
  std::vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Canvas");
  std::vector<std::string> token = {"Path","Divide","Histo"};
  gROOT->ProcessLine("gROOT->SetBatch(kTRUE)");
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      std::vector<std::string> path = blocks[i]->GetVectorString("Path");
      std::vector<int> divide = blocks[i]->GetVectorInt("Divide");
      std::vector<std::string> histo = blocks[i]->GetVectorString("Histo");
      std::string name = path[path.size()-1];
      TCanvas* c = new TCanvas(name.c_str(), 5000,5000,0);
      c->Divide(divide[0],divide[1]);

      unsigned int size = histo.size();
      for(unsigned int h = 0 ; h < size ; h++){
        c->cd(h+1);
        std::string padname=name+"_"+NPL::itoa(h);
        gPad->SetName(padname.c_str());
        TH1* hist = (TH1*) Spectra->FindObject(histo[h].c_str());
        if(hist){
          hist->UseCurrentStyle();
          hist->Draw("colz");
          ExecuteMacro(hist->GetName());
        }
      }
      if(m_EmbeddedCanvas && m_OldCurrent == c->GetName()){
        m_EmbeddedCanvas->AdoptCanvas(c);
      }


      TGListTreeItem*  item  =  NULL;
      TGListTreeItem*  pitem =  NULL;

      std::string item_path="";
      for(unsigned int j = 0 ; j < path.size()-1 ; j++){
        item_path+="/"+path[j];
        item = m_ListTree->FindItemByPathname(item_path.c_str());
        if(!item){
          item= m_ListTree->AddItem(pitem,path[j].c_str());
          if(pitem)
            pitem->SetPictures(m_pfolder,m_pfolder); 
        }
        pitem = item;

      }

      if(item)
        AddItem(c,item);
      else
        AddItem(c);
    }
    else
      NPL::SendWarning("NPL::CanvasList","CanvasList.txt has incorrect formatting");
  }
  gROOT->ProcessLine("gROOT->SetBatch(kFALSE)");

}
