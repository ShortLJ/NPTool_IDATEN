/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiFDC2 treated data                                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TSamuraiFDC2Physics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPOptionManager.h"
#include "NPDetectorFactory.h"
#include "NPSystemOfUnits.h"
//   ROOT
using namespace NPUNITS;
///////////////////////////////////////////////////////////////////////////

ClassImp(TSamuraiFDC2Physics)
  ///////////////////////////////////////////////////////////////////////////
  TSamuraiFDC2Physics::TSamuraiFDC2Physics(){
    m_EventData         = new TSamuraiFDC2Data ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
 //   ToTThreshold_L = 180;
 //   ToTThreshold_H = 1000;
    ToTThreshold_L = 0;
    ToTThreshold_H = 10000; 
    DriftLowThreshold=0.1;
    DriftUpThreshold=9.9;
    //PowerThreshold=15;

    PowerThreshold=1e6;
  }

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::BuildPhysicalEvent(){
  PreTreat();
//  RemoveNoise();

  // Map[plane angle, vector of spatial information]
  static map<double, vector<double> > X ; 
  static map<double, vector<double> > Z ; 
  static map<double, vector<double> > R ; 
  static int det,layer,wire;
  X.clear();Z.clear();R.clear();

  unsigned int size = Detector.size();
  for(unsigned int i = 0 ; i < size ; i++){
    if(DriftLength[i] > DriftLowThreshold && DriftLength[i] < DriftUpThreshold){
      det = Detector[i];
      layer = Layer[i];
      wire = Wire[i]; 
      SamuraiDCIndex idx(det,layer,wire);
      X[Wire_Angle[idx]].push_back(Wire_X[idx]); 
      Z[Wire_Angle[idx]].push_back(Wire_Z[idx]); 
      R[Wire_Angle[idx]].push_back(DriftLength[i]); 
    }
  }

  // Reconstruct the vector for each of the plane of each of the detector
  static double X0,X100,a,b; // store the BuildTrack2D results
  static map<double, TVector3 > VX0 ;  
  static map<double, TVector3 > VX100 ;  
  static map<double, double > D ;// the minimum distance  
  static unsigned int uid; uid=0;

  VX0.clear();VX100.clear(),D.clear();

  for(auto it = X.begin();it!=X.end();++it){
#if __cplusplus > 199711L && NPMULTITHREADING
    m_reconstruction.AddPlan(uid++,X[it->first],Z[it->first],R[it->first]); 
#else
    D[it->first]=m_reconstruction.BuildTrack2D(X[it->first],Z[it->first],R[it->first],X0,X100,a,b); 
#endif 
    }

#if __cplusplus > 199711L && NPMULTITHREADING
  // do all plan at once in parallele, return when all plan are done
  m_reconstruction.BuildTrack2D();
  uid=0;
#endif
  
  for(auto it = X.begin();it!=X.end();++it){
#if __cplusplus > 199711L && NPMULTITHREADING
 
  D[it->first]=m_reconstruction.GetResults(uid++,X0,X100,a,b); 
#endif

/*   // for Debug, write a file of 
   { std::ofstream f("distance.txt", std::ios::app);
   f<< D[it->first] << endl;
   f.close();
   }
*/    
   // very large a means track perpendicular to the chamber, what happen when there is pile up
   if(abs(a)>5000)
      PileUp++;

    Mult+=X[it->first].size();
    // Position at z=0
    TVector3 P(X0,0,0);
    P.RotateZ(it->first);
    VX0[it->first]=P;
    // Position at z=100
    TVector3 P100= TVector3(X100,0,0);
    P100.RotateZ(it->first);
    VX100[it->first]=P100;
  }

  // Reconstruct the central position (z=0) for each detector
  static vector<TVector3>  C ;  
  static vector<double>    W ; // weight based on D  
  C.clear(),W.clear();
  TVector3 P;

  for(auto it1 = VX0.begin();it1!=VX0.end();++it1){
    for(auto it2 = it1;it2!=VX0.end();++it2){
      if(it1!=it2){// different plane, same detector
        m_reconstruction.ResolvePlane(it1->second,it1->first,it2->second,it2->first,P);
        if(P.X()!=-10000 /*&& D[it1->first]<PowerThreshold && D[it2->first]<PowerThreshold*/){
          C.push_back(P);
          // Mean pos are weighted based on the the sum of distance from track
          // to hit obtained during the minimisation
          W.push_back(1./sqrt(D[it1->first]*D[it2->first]));
          }
      }
    }
  }
  
  // Reconstruct the position at z=100 for each detector
  static vector<TVector3> C100 ;  
  C100.clear();
  for(auto it1 = VX100.begin();it1!=VX100.end();++it1){
    for(auto it2 = it1;it2!=VX100.end();++it2){
      if(it1!=it2){// different plane
        m_reconstruction.ResolvePlane(it1->second,it1->first,it2->second,it2->first,P);

        if(P.X()!=-10000/*&& D[it1->first]<PowerThreshold && D[it2->first]<PowerThreshold*/)
          C100.push_back(P);
      }
    }
  }

  // Build the Reference position by averaging all possible pair 
  size = C.size();
  static double PosX100,PosY100,norm;
  if(size){
    
    PosX=0;
    PosY=0;
    PosX100=0;
    PosY100=0;
    norm=0;
    for(unsigned int i = 0 ; i < size ; i++){
      PosX+= C[i].X()*W[i]; 
      PosY+= C[i].Y()*W[i]; 
      PosX100+= C100[i].X()*W[i]; 
      PosY100+= C100[i].Y()*W[i]; 
      norm+=W[i];
      //        cout << C[2][i].X() << " (" << C[2][i].Y() << ") ";
    } 
    // cout << endl;
    MultMean=size;
    // Mean position at Z=0
    PosX=PosX/norm; 
    PosY=PosY/norm; 
    // Mean position at Z=100
    PosX100=PosX100/norm; 
    PosY100=PosY100/norm; 
    
    devX=0;
    devY=0;
    for(unsigned int i = 0 ; i < size ; i++){
      devX+=W[i]*(C[i].X()-PosX)*(C[i].X()-PosX);
      devY+=W[i]*(C[i].Y()-PosY)*(C[i].Y()-PosY);
    }

    devX=sqrt(devX/((size-1)*norm));
    devY=sqrt(devY/((size-1)*norm));
   
    if(m_invertX){
      PosX*=-1;
      PosX100*=-1;
    }

    if(m_invertY){
      PosY*=-1;
      PosY100*=-1;
    }

    // Compute ThetaX, angle between the Direction vector projection in XZ with
    // the Z axis
    //ThetaX=atan((PosX100-PosX)/100.);
    ThetaX = (PosX100-PosX)/100.;
    // Compute PhiY, angle between the Direction vector projection in YZ with
    // the Z axis
    //PhiY=atan((PosY100-PosY)/100.);
    PhiY=(PosY100-PosY)/100.;
    Dir=TVector3(PosX100-PosX,PosY100-PosY,100).Unit();
    PosX+=m_offset.X();
    PosY+=m_offset.Y();
  }
/*  if(PosX==-10000)
    cout << " bad " <<  Detector.size()<< " " << size << endl;
  else
    cout << " okay" <<  Detector.size()<< " " << size << endl;*/
  return;
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::PreTreat(){
  static CalibrationManager* Cal = CalibrationManager::getInstance();
  static string channel;
  //m_EventData->Print();
  unsigned int size = m_EventData->Mult();
  for(unsigned int i = 0 ; i < size ; i++){
    // EDGE=1 is the leading edge, IE, the real time.
    // EDGE=0 is the trailing edge, so it helps build Tot
    if(m_EventData->GetEdge(i)==1){
      int det   = m_EventData->GetDetectorNbr(i); 
      int layer = m_EventData->GetLayerNbr(i); 
      int wire  = m_EventData->GetWireNbr(i); 
      double time = m_EventData->GetTime(i);
      double etime = 0;
      // look for matching trailing edge   
      for(unsigned int j = 0 ; j < size ; j++){
        if(m_EventData->GetEdge(j)==0){
          int edet   = m_EventData->GetDetectorNbr(j); 
          int elayer = m_EventData->GetLayerNbr(j); 
          int ewire  = m_EventData->GetWireNbr(j); 
          // same wire
          if(wire==ewire && layer==elayer && det==edet){
            etime = m_EventData->GetTime(j); 
          }    
        }
        if(etime && etime>time)
          break;
        else
          etime=0;
      }
      // a valid wire must have an edge
      if(etime && time && etime-time>ToTThreshold_L && etime-time<ToTThreshold_H){
        if(!(wire==93 && layer ==7)){// remove noisy wire
         Detector.push_back(det);
          Layer.push_back(layer);       
          Wire.push_back(wire);
          Time.push_back(time);
          ToT.push_back(etime-time);
          channel="SamuraiFDC2/L" + NPL::itoa(layer);
          // rescalling is needed because calib are bad.
          // to be fixed
          if(!m_invertD)
            DriftLength.push_back(Cal->ApplySigmoid(channel,etime));
          else
            DriftLength.push_back(10-Cal->ApplySigmoid(channel,etime));
        }
      }
    }
  }
  return;
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::RemoveNoise(){
  // Remove the noise by looking if a matching wire exist in the adjacent layer
  // this done by looking at the closest plane with the same orientation
  unsigned int size = Detector.size(); 
  for(unsigned int i = 0 ; i < size ; i++){
    bool match=false;
    int det = Detector[i];
    int layer = Layer[i];
    int wire = Wire[i];
    // look for matching adjacent wire   

    for(unsigned int j = 0 ; j < size ; j++){
      int adet = Detector[j];
      int alayer = Layer[j];
      int awire = Wire[j];
      bool blayer = false;
      if(layer%2==0){
        if(layer+1==alayer)
          blayer=true;
      }

      else{
        if(layer-1==alayer)
          blayer=true;
      }

      if(det==adet && blayer && abs(wire-awire)<=1){
        match=true;
        break;
      }
    }

    if(match)
      Matched.push_back(true);
    else
      Matched.push_back(false);
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////
TVector3 TSamuraiFDC2Physics::ProjectedPosition(double Z){
  TVector3 pos(-10000,-10000,-10000);
  if(PosX!=-10000){
  pos = TVector3(PosX,PosY,0)+(Z/Dir.Z())*Dir;
  //cout << pos.X() << " " << pos.Y() << " " << pos.Z() << endl;
  }
  return pos;
}
////////////////////////////////////////////////////////////////////////////////
double TSamuraiFDC2Physics::ProjectedPositionX(double Z){
  return ProjectedPosition(Z).X();
}
////////////////////////////////////////////////////////////////////////////////
double TSamuraiFDC2Physics::ProjectedPositionY(double Z){
  return ProjectedPosition(Z).Y();
}
   
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::Clear(){
  MultMean=0;
  PileUp=0;
  Mult=0;
  PosX=PosY=-10000;
  ThetaX=PhiY=-10000;
  devX=devY=-10000;
  DriftLength.clear();
  Detector.clear();
  Layer.clear();
  Wire.clear();
  Time.clear();
  ToT.clear();
  ParticleDirection.clear();
  MiddlePosition.clear();
  Matched.clear();
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SAMURAIFDC2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML","Offset","InvertX","InvertY","InvertD"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    cout << endl << "////  Samurai FDC2 (" << i+1 << ")" << endl;
    string xmlpath = blocks[i]->GetString("XML");
    NPL::XmlParser xml;
    xml.LoadFile(xmlpath);
    AddDC("SAMURAIFDC2",xml);
    m_offset = blocks[i]->GetTVector3("Offset","mm"); 
    m_invertX = blocks[i]->GetInt("InvertX"); 
    m_invertY = blocks[i]->GetInt("InvertY"); 
    m_invertD = blocks[i]->GetInt("InvertD"); 
  }


#if __cplusplus > 199711L && NPMULTITHREADING
 
  if(blocks.size()){
    // one thread for each plan X,U,V = 3
    // ! more than this will not help !
    m_reconstruction.SetNumberOfThread(3);
    m_reconstruction.InitThreadPool();
   }
#endif 

    GetOffset().Print();
    PosX=1;
    cout << m_invertY << endl;
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::AddDC(string name, NPL::XmlParser& xml){
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName(name);  
  // FDC2 case
  if(name=="SAMURAIFDC2"){
    unsigned int det=2;
    unsigned int size = b.size();
    for(unsigned int i = 0 ; i < size ; i++){
      unsigned int layer = b[i]->AsInt("layer"); 
      unsigned int wire  = b[i]->AsInt("wireid"); 
      double X = b[i]->AsDouble("wirepos");  
      double Z = b[i]->AsDouble("wirez");  
      string sDir = b[i]->AsString("anodedir");
      double T=0;
      if(sDir=="X")
        T=0*deg;
      else if(sDir=="Y")
        T=90*deg;
      else if(sDir=="U")
        T=-30*deg;
      else if(sDir=="V")
        T=+30*deg;
      else{
        cout << "ERROR: Unknown layer orientation for Samurai FDC2"<< endl;
        exit(1);
      }
      SamuraiDCIndex idx(det,layer,wire);
      Wire_X[idx]=X;
      Wire_Z[idx]=Z;
      Wire_Angle[idx]=T;
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitSpectra(){  
  //m_Spectra = new TSamuraiFDC2Spectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::FillSpectra(){  
  //  m_Spectra -> FillRawSpectra(m_EventData);
  //  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::CheckSpectra(){  
  //  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TSamuraiFDC2Physics::GetSpectra() {
  /*  if(m_Spectra)
      return m_Spectra->GetMapHisto();
      else{
      map< string , TH1*> empty;
      return empty;
      }*/
  map< string , TH1*> empty;
  return empty;

} 

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::WriteSpectra(){
  // m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::AddParameterToCalibrationManager(){
  CalibrationManager* Cal = CalibrationManager::getInstance();

  // each layer
  for( int l = 0 ; l < 14 ; ++l){
    Cal->AddParameter("SamuraiFDC2", "L"+ NPL::itoa(l),"FDC2_L"+ NPL::itoa(l));
  }

}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "SamuraiFDC2" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fDC_*"))
    inputChain->SetBranchStatus( "fDC_*",true);
  inputChain->SetBranchAddress( "SamuraiFDC2" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitializeRootInputPhysics(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "SamuraiFDC2" , true );
  inputChain->SetBranchAddress( "SamuraiFDC2" , &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiFDC2Physics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree("SamuraiFDC2");
  outputTree->Branch( "SamuraiFDC2" , "TSamuraiFDC2Physics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSamuraiFDC2Physics::Construct(){
  return (NPL::VDetector*) new TSamuraiFDC2Physics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_samuraiFDC2{
    public:
      proxy_samuraiFDC2(){
        NPL::DetectorFactory::getInstance()->AddToken("SAMURAIFDC2","Samurai");
        NPL::DetectorFactory::getInstance()->AddDetector("SAMURAIFDC2",TSamuraiFDC2Physics::Construct);
      }
  };

  proxy_samuraiFDC2 p_samuraiFDC2;
}

