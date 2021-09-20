/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiBDC treated data                                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TSamuraiBDCPhysics.h"

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
#include "TAsciiFile.h"
#include "NPOptionManager.h"
#include "NPDetectorFactory.h"
#include "NPSystemOfUnits.h"
//   ROOT
using namespace NPUNITS;
///////////////////////////////////////////////////////////////////////////

ClassImp(TSamuraiBDCPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TSamuraiBDCPhysics::TSamuraiBDCPhysics(){
    m_EventData         = new TSamuraiBDCData ;
    m_PreTreatedData    = new TSamuraiBDCData ;
    m_EventPhysics      = this ;
    //m_Spectra           = NULL;
    ToTThreshold_L = 0;
    ToTThreshold_H = 1000;
    DriftLowThreshold=0 ;
    DriftUpThreshold=2.5;
    PowerThreshold=5;
  }

///////////////////////////////////////////////////////////////////////////
TVector3 TSamuraiBDCPhysics::GetPos(unsigned int det){
  TVector3 res(-10000,-10000,-10000); 
  unsigned int size = PosX.size();
  for(unsigned int i = 0 ; i < size ; i++){
    if(Detector[i]==det){
      res = TVector3(PosX[i],PosY[i],PosZ[i]);
    }
  }
  return res;


}
///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::BuildPhysicalEvent(){
  PreTreat();

  static unsigned int det,layer,wire,size;
  static double dl;
  // Map[plane angle, vector of spatial information]
  static map<double, vector<double> > X ; 
  static map<double, vector<double> > Z ; 
  static map<double, vector<double> > R ; 
  static double X0,X100,a,b; // store the BuildTrack2D results
  static map<double, TVector3 > VX0 ;  
  static map<double, TVector3 > VX100 ;  
  static map<double, double > D ;// the minimum distance  
  static unsigned int uid; uid=0;
  static vector<TVector3> C ;  
  static vector<double  > W ; // weight based on D  
  static double PosX100,PosY100,norm;
  int count = 0 ;
  for(auto it = m_DCHit.begin(); it!=m_DCHit.end(); it++){
    // Each entry in the map is a detector 
    det = it->first;
    Detector.push_back(det);
    PosX.push_back(0);
    PosY.push_back(0);
    PosZ.push_back(0);
    ThetaX.push_back(0);
    PhiY.push_back(0);
    devX.push_back(0);
    devY.push_back(0);
    Dir.push_back(TVector3());
    PileUp.push_back(0);


    X.clear();Z.clear();R.clear();
    // Build the necessary X,Z,R vector
    for(auto itt = it->second.begin() ; itt!= it->second.end() ; itt++){
      dl = (*itt).DriftLength; 
      if(dl > DriftLowThreshold && dl < DriftUpThreshold){
        SamuraiDCIndex idx(det,(*itt).Layer,(*itt).Wire);
        X[Wire_Angle[idx]].push_back(Wire_X[idx]); 
        Z[Wire_Angle[idx]].push_back(Wire_Z[idx]); 
        R[Wire_Angle[idx]].push_back(dl);
      }
    }
    // Reconstruct the vector for each of the plane of each of the detector
    VX0.clear();VX100.clear(),D.clear();


    uid=0;
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

    // Loop over the results
    for(auto it = X.begin();it!=X.end();++it){
#if __cplusplus > 199711L && NPMULTITHREADING
      D[it->first]=m_reconstruction.GetResults(uid++,X0,X100,a,b); 
#endif
      // for Debug, write a file of 
      //  { std::ofstream f("distance.txt", std::ios::app);
      //  f<< D[it->first] << endl;
      //  f.close();
      //  }

      // very large "a" means track perpendicular to the chamber, what happen when there is pile up
      if(abs(a)>5000)
        PileUp[count]++;
//cout << a << " " << b << endl;
      // Position at z=0
      TVector3 P(X0,0,0);
      P.RotateZ(it->first);
      VX0[it->first]=P;
      // Position at z=100
      TVector3 P100= TVector3(X100,0,0);
      P100.RotateZ(it->first);
      VX100[it->first]=P100;

      // Reconstruct the central position (z=0) for each detector
      C.clear(),W.clear();
      for(auto it1 = VX0.begin();it1!=VX0.end();++it1){
        for(auto it2 = it1;it2!=VX0.end();++it2){
          if(it1!=it2){// different plane
            //cout << "BDC" << endl;
            m_reconstruction.ResolvePlane(it1->second,it1->first,it2->second,it2->first,P);
            // cout << "done " << D[it1->first] << " " << D[it2->first] << endl;
            if(P.X()!=-10000 /*&& D[it1->first]<PowerThreshold&& D[it2->first]<PowerThreshold*/){
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
          if(it1!=it2 ){// different plane, same detector
            m_reconstruction.ResolvePlane(it1->second,it1->first,it2->second,it2->first,P);

            if(P.X()!=-10000/*&& D[it1->first]<PowerThreshold && D[it2->first]<PowerThreshold*/)
              C100.push_back(P);
          }
        }
      }
      // Build the Reference position by averaging all possible pair 
      size = C.size();
      if(size){
        norm=0;PosX100=0;PosY100=0;
        for(unsigned int i = 0 ; i < size ; i++){
          PosX[count]+= C[i].X()*W[i]; 
          PosY[count]+= C[i].Y()*W[i]; 
          PosX100+= C100[i].X()*W[i]; 
          PosY100+= C100[i].Y()*W[i]; 
          norm+=W[i];
        } 
        //MultMean=size;
        // Mean position at Z=0
        PosX[count]/=norm; 
        PosY[count]/=norm; 
        // Mean position at Z=100
        PosX100/=norm; 
        PosY100/=norm; 
        
        for(unsigned int i = 0 ; i < size ; i++){
          devX[count]+=W[i]*(C[i].X()-PosX[count])*(C[i].X()-PosX[count]);
          devY[count]+=W[i]*(C[i].Y()-PosY[count])*(C[i].Y()-PosY[count]);
        }
        
        devX[count]=sqrt(devX[count]/((size-1)*norm));
        devY[count]=sqrt(devY[count]/((size-1)*norm));

        // Compute ThetaX, angle between the Direction vector projection in XZ with
        // the Z axis
        ThetaX[count]=(PosX100-PosX[count])/100.;
        // Compute PhiY, angle between the Direction vector projection in YZ with
        // the Z axis
        PhiY[count]=(PosY100-PosY[count])/100.;
        Dir[count]=TVector3(PosX100-PosX[count],PosY100-PosY[count],100).Unit();
        if(m_invertX[det])
          PosX[count]*=-1;
        if(m_invertY[det])
          PosY[count]*=-1;
        
        PosX[count]+=m_offset[det].X();
        PosY[count]+=m_offset[det].Y();
        PosZ[count]=m_offset[det].Z();
      }
    }

    if(PosX[count]==0&&PosY[count]==0){
      PosX.erase(PosX.begin()+count);
      PosY.erase(PosY.begin()+count);
      PosZ.erase(PosZ.begin()+count);
      ThetaX.erase(ThetaX.begin()+count);
      PhiY.erase(PhiY.begin()+count);
      Detector.erase(Detector.begin()+count);
      count--;
    }
    count++;

  }// detector loop
  return;
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::PreTreat(){
  ClearPreTreatedData();
  static CalibrationManager* Cal = CalibrationManager::getInstance();
  static string channel;

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
        channel="SamuraiBDC"+NPL::itoa(det)+"/L" + NPL::itoa(layer);
        SamuraiDCIndex idx(det,layer,wire);
        if(!m_invertD[det])
          m_DCHit[det].push_back(DCHit(det,layer,wire,time,etime-time,Cal->ApplySigmoid(channel,etime)));
        else
          m_DCHit[det].push_back(DCHit(det,layer,wire,time,etime-time,2.5-Cal->ApplySigmoid(channel,etime)));
      }
    }
  }
  return;
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::Clear(){
  m_DCHit.clear();
  // Computed variable
  PosX.clear();
  PosY.clear();
  PosZ.clear();
  ThetaX.clear();
  PhiY.clear();
  devX.clear();
  devY.clear();
  Dir.clear();
  PileUp.clear();
  Detector.clear();
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SAMURAIBDC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML","Offset","InvertX","InvertY","InvertD"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      cout << endl << "////  Samurai BDC (" << i+1 << ")" << endl;
      unsigned int det = std::atoi(blocks[i]->GetMainValue().c_str());
      string xmlpath = blocks[i]->GetString("XML");
      NPL::XmlParser xml;
      xml.LoadFile(xmlpath);
      AddDC(det,xml);
      TVector3 offset = blocks[i]->GetTVector3("Offset","mm"); 
      bool invertX = blocks[i]->GetInt("InvertX"); 
      bool invertY = blocks[i]->GetInt("InvertY"); 
      bool invertD = blocks[i]->GetInt("InvertD"); 
      m_offset[det] = offset;
      m_invertX[det] = invertX;
      m_invertY[det] = invertY;
      m_invertD[det] = invertD;
    }
    else{
      cout << " --- ERROR : BDC block wrongly formatted" << endl;
      exit(1);
    }
  }

#if __cplusplus > 199711L && NPMULTITHREADING 
  if(blocks.size()){ // if a detector is found, init the thread pool
    // one thread for each plan X,Y = 2
    // ! more than that this will not help !
    m_reconstruction.SetNumberOfThread(2);
    m_reconstruction.InitThreadPool();
  }
#endif 
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::AddDC(int det, NPL::XmlParser& xml){
  std::string name = "SAMURAIBDC"+NPL::itoa(det);
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName(name);  
  unsigned int sizeB = b.size();
  for(unsigned int i = 0 ; i < sizeB ; i++){
    unsigned int layer = b[i]->AsInt("layer"); 
    unsigned int wire  = b[i]->AsInt("wireid"); 
    double X = b[i]->AsDouble("wirepos");  
    double Z = b[i]->AsDouble("wirez");  
    string sDir = b[i]->AsString("anodedir");
    double T=0;
    if(sDir=="X")
      T= 0*deg;
    else if(sDir=="Y")
      T= 90*deg;
    else if(sDir=="U")
      T=-30*deg;
    else if(sDir=="V")
      T=+30*deg;
    else{
      cout << "ERROR: Unknown layer orientation for Samurai BDC"<< endl;
      exit(1);
    }
    SamuraiDCIndex idx(det,layer,wire);
    Wire_X[idx]=X;
    Wire_Z[idx]=Z;
    Wire_Angle[idx]=T;
  }
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::InitSpectra(){  
  //m_Spectra = new TSamuraiBDCSpectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::FillSpectra(){  
  //  m_Spectra -> FillRawSpectra(m_EventData);
  //  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::CheckSpectra(){  
  //  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TSamuraiBDCPhysics::GetSpectra() {
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
void TSamuraiBDCPhysics::WriteSpectra(){
  // m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::AddParameterToCalibrationManager(){
  CalibrationManager* Cal = CalibrationManager::getInstance();

  // for each det
  for( int d = 1 ; d < 3 ; ++d){
    // each layer
    for( int l = 0 ; l < 8 ; ++l){
      Cal->AddParameter("SamuraiBDC"+NPL::itoa(d), "L"+ NPL::itoa(l),"BDC"+NPL::itoa(d)+"_L"+ NPL::itoa(l));
    }
  }

}

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "SamuraiBDC" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fBDC_*"))
    inputChain->SetBranchStatus( "fBDC_*",true);
  inputChain->SetBranchAddress( "SamuraiBDC" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::InitializeRootInputPhysics(){
}

///////////////////////////////////////////////////////////////////////////
void TSamuraiBDCPhysics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree("SamuraiBDC");
  outputTree->Branch( "SamuraiBDC" , "TSamuraiBDCPhysics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSamuraiBDCPhysics::Construct(){
  return (NPL::VDetector*) new TSamuraiBDCPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_samuraiBDC{
    public:
      proxy_samuraiBDC(){
        NPL::DetectorFactory::getInstance()->AddToken("SAMURAIBDC","Samurai");
        NPL::DetectorFactory::getInstance()->AddToken("SAMURAIBDC1","Samurai");
        NPL::DetectorFactory::getInstance()->AddToken("SAMURAIBDC2","Samurai");
        NPL::DetectorFactory::getInstance()->AddDetector("SAMURAIBDC",TSamuraiBDCPhysics::Construct);
        NPL::DetectorFactory::getInstance()->AddDetector("SAMURAIBDC1",TSamuraiBDCPhysics::Construct);
        NPL::DetectorFactory::getInstance()->AddDetector("SAMURAIBDC2",TSamuraiBDCPhysics::Construct);
      }
  };

  proxy_samuraiBDC p_samuraiBDC;
}

