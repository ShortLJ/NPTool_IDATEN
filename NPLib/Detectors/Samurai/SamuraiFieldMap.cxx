/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : May  2021                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Samurai field map data                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "SamuraiFieldMap.h"
#include "NPPhysicalConstants.h"
#include "Math/Factory.h"
using namespace NPUNITS;

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(SamuraiFieldMap);

SamuraiFieldMap::SamuraiFieldMap(){
  m_BrhoScan=NULL;
  m_min=ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
  m_func=ROOT::Math::Functor(this,&SamuraiFieldMap::Delta,1); 
  m_min->SetFunction(m_func);
  m_min->SetPrintLevel(0);



}
////////////////////////////////////////////////////////////////////////////////
double SamuraiFieldMap::Delta(const double* parameter){
    static vector<TVector3>pos ;
    static TVector3 diff;
    pos =Propagate(parameter[0],m_FitPosFDC0,m_FitDirFDC0,false); 
    // Move the fdc2 pos from lab frame to fdc2 frame 
    pos.back().RotateY(-m_fdc2angle+m_angle); 

    //double d = (pos.back().X()-m_FitPosFDC2.X())*(pos.back().X()-m_FitPosFDC2.X());
   // return d;
    diff = pos.back()-m_FitPosFDC2;
    return diff.Mag2();
}

////////////////////////////////////////////////////////////////////////////////
double SamuraiFieldMap::FindBrho(TVector3 p_fdc0,TVector3 d_fdc0,TVector3 p_fdc2,TVector3 d_fdc2){
  m_FitPosFDC0=p_fdc0;
  m_FitDirFDC0=d_fdc0;
  m_FitPosFDC2=p_fdc2;
  m_FitDirFDC2=d_fdc2;

  if(!m_BrhoScan)
    BrhoScan(1,10,0.1);
  // do a first guess based on fdc2 pos
  double b0[1] ={m_BrhoScan->Eval(p_fdc2.X())}; 
  
  m_min->Clear(); 
  m_min->SetPrecision(1e-6);
  m_min->SetMaxFunctionCalls(1000);
  m_min->SetLimitedVariable(0,"B",b0[0],0.1,1,10);
  m_min->Minimize(); 
  return m_min->X()[0];
}

////////////////////////////////////////////////////////////////////////////////
TGraph* SamuraiFieldMap::BrhoScan(double min, double max,double step){
  if(m_BrhoScan)
    delete m_BrhoScan;
  m_BrhoScan=new TGraph;
  unsigned int size = (max-min)/step;
  m_BrhoScan->Set(size);
  unsigned int i=0;
  TVector3 p(0,0,-3500);
  TVector3 d(0,0,1);
  p.RotateY(m_angle);
  d.RotateY(m_angle);
  for(double b = min ; b < max ; b+=step){
    vector<TVector3> pos= Propagate(b,p,d,false);
    pos.back().RotateY(-m_fdc2angle);
    m_BrhoScan->SetPoint(i++,pos.back().X(),b); 
  }
  m_BrhoScan->Sort();
  return m_BrhoScan;
}

////////////////////////////////////////////////////////////////////////////////
TVector3 SamuraiFieldMap::PropagateToFDC2(TVector3 pos, TVector3 dir){
  // go to FDC2 frame reference
  pos.RotateY(-m_fdc2angle);
  dir.RotateY(-m_fdc2angle);

  double deltaZ=m_fdc2R-pos.Z();
  dir*=deltaZ/dir.Z();
  pos+=dir;
  pos.SetX(pos.X());
  pos.RotateY(m_fdc2angle);
  return pos;
}

////////////////////////////////////////////////////////////////////////////////
std::vector< TVector3 > SamuraiFieldMap::Propagate(double Brho, TVector3 pos, TVector3 dir,bool store){
  pos.RotateY(m_angle);
  dir.RotateY(m_angle);
  dir=dir.Unit();
  // Property of a particle with the correct Brho:
  // We assume a 4He to compute v
  // The choice of the particle is of no importance
  static NPL::Particle N("4He");
  N.SetBrho(Brho);

  // track result
  static std::vector< TVector3 > track;
  track.clear();

  // starting point of the track
  if(store){
    pos.RotateY(-m_angle);
    track.push_back(pos);
    pos.RotateY(m_angle);
  }
  dir=dir.Unit();
  static double r;
  r = sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
  // number of step taken
  static unsigned int count,limit;
  count = 0;
  // maximum number of state before giving up
  limit = 1000;

  // First propagate to r_max with one line
  while(r>m_Rmax && count<limit){
    pos+=(r-m_Rmax)/cos(dir.Theta())*dir.Unit();
    r= 1.01*sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
  }

  if(r<=m_Rmax){ // success
    if(store){
      pos.RotateY(-m_angle);
      track.push_back(pos);
      pos.RotateY(m_angle);
    }
  }
  else {// failure
    //cout << "Fail" << endl;
    return track;
  }

  static TVector3 xk1,xk2,xk3,xk4; // position
  static TVector3 pk1,pk2,pk3,pk4; // impulsion
  static TVector3 imp;
  static double K,m,P,px,py,pz;
  K = N.GetEnergy(); // kinetic energy
  m = N.Mass(); // mc2
  P = sqrt(K*K+2*K*m)/c_light; // P
  px = P*dir.X();//px
  py = P*dir.Y();//py
  pz = P*dir.Z();//pz
  imp = P*dir;
  static double h = 1*nanosecond;
  while(r<=m_Rmax && count < limit){
    func(N, pos           , imp            , xk1, pk1);
    func(N, pos+xk1*(h/2.), imp+pk1*(h/2.) , xk2, pk2);
    func(N, pos+xk2*(h/2.), imp+pk2*(h/2.) , xk3, pk3);
    func(N, pos+xk3*h     , imp+pk3*h      , xk4, pk4);
    pos +=(xk1+2*xk2+2*xk3+xk4)*(h/6.); 
    imp +=(pk1+2*pk2+2*pk3+pk4)*(h/6.); 
    if(store){
      pos.RotateY(-m_angle);
      track.push_back(pos);
      pos.RotateY(m_angle);
    }
    r = sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
    count++;
  }
  imp=imp.Unit();
  pos = PropagateToFDC2(pos, imp);
  pos.RotateY(-m_angle);
  track.push_back(pos);

  return track;

}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::func(NPL::Particle& N, TVector3 pos, TVector3 imp, TVector3& new_pos, TVector3& new_imp){
  static double px,py,pz,vx,vy,vz,Bx,By,Bz,q,P2,D,m2c4;
  static vector<double> B; 
  px=imp.X(); 
  py=imp.Y();
  pz=imp.Z();

  P2=imp.Mag2(); // P2
  m2c4 = N.Mass()*N.Mass();
  D=sqrt(m2c4+P2*c_squared); // sqrt(m2c4+P2c2)
  vx=px*c_squared/D;// pxc * c / D = pxc2/D
  vy=py*c_squared/D;
  vz=pz*c_squared/D;
  new_pos.SetX(vx);
  new_pos.SetY(vy);
  new_pos.SetZ(vz);
  B = InterpolateB(pos);
  Bx= B[0]; 
  By= B[1];
  Bz= B[2];
  q = N.GetZ()*eplus; // issue with the tesla/coulomb definition
  new_imp.SetX(q*(vy*Bz-vz*By));// q*pyc2*Bz/D -q*pzc2*By/D
  new_imp.SetY(q*(vz*Bx-vx*Bz));
  new_imp.SetZ(q*(vx*By-vy*Bx));
}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadMap(double angle,std::string file,unsigned int bin){
  m_bin=bin;
  m_angle=angle;
  if(file.find(".bin")!=std::string::npos)
    LoadBinary(file);
  else
    LoadAscii(file);
}
////////////////////////////////////////////////////////////////////////////////
std::vector<double> SamuraiFieldMap::GetB(std::vector<double>& pos){
  static vector<double> nullv ={0,0,0};
  // the map is only 1/4 of the detector so we apply symetrie:
  double x,y,z ;

  if(pos[0]<0)
    pos[0] = -pos[0];

  if(pos[2]<0)
    pos[2] = -pos[2];

  auto it=m_field.find(pos);
  if(it!=m_field.end()){
    return it->second;
  }
  else 
    return nullv;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<double> SamuraiFieldMap::InterpolateB(const std::vector<double>& pos){
  static vector<double> nullv ={0,0,0};
  // the map is only 1/4 of the detector so we apply symetrie:
  double x,y,z ;

  if(pos[0]>0)
    x = pos[0];
  else
    x = -pos[0];

  y = pos[1];

  if(pos[2]>0)
    z = pos[2];
  else
    z = -pos[2];

  // out of bound 
  if(x<m_x_min || x>m_x_max)
    return nullv;
  if(y<m_y_min || y>m_y_max)
    return nullv;
  if(z<m_z_min || z>m_z_max)
    return nullv;



  double xm = (double)((int)x/m_bin*m_bin);
  double ym = (double)((int)y/m_bin*m_bin);
  double zm = (double)((int)z/m_bin*m_bin);

  vector<double> p0={xm,ym,zm};
  vector<double> p1={xm+m_bin,ym,zm};
  vector<double> p2={xm,ym+m_bin,zm};
  vector<double> p3={xm,ym,zm+m_bin};
  vector<double> p4={xm-m_bin,ym,zm};
  vector<double> p5={xm,ym-m_bin,zm};
  vector<double> p6={xm,ym,zm-m_bin};

  vector<map<vector<double>,vector<double>>::iterator> it=
  { m_field.lower_bound(p0),
    m_field.lower_bound(p1),m_field.lower_bound(p2),m_field.lower_bound(p3),
    m_field.lower_bound(p4),m_field.lower_bound(p5),m_field.lower_bound(p6)};

  double Bx=0;
  double By=0;
  double Bz=0;
  double totalW=0;
  auto end=m_field.end();
  unsigned int size = it.size();
  for(unsigned int i = 0 ; i < size; i++){
    if(it[i]!=end){
      double d = 1e-6+sqrt( (x-it[i]->first[0])*(x-it[i]->first[0])+
          (y-it[i]->first[1])*(y-it[i]->first[1])+
          (z-it[i]->first[2])*(z-it[i]->first[2]));

      Bx+=it[i]->second[0]/(d*d);
      By+=it[i]->second[1]/(d*d);
      Bz+=it[i]->second[2]/(d*d);
      totalW+=1./(d*d);
    }
  }
  vector<double> res = {Bx/totalW,By/totalW,Bz/totalW};
  return res;
}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadAscii(std::string file){
  ifstream in(file.c_str());
  if(!in.is_open()){
    cout << "Error: failed to load samurai field map " << file << endl;
    exit(1);
  }

  cout << "//////// Loading Ascii Samurai field map " << file << endl; 
  double x,y,z,Bx,By,Bz;

  m_x_max=m_y_max=m_z_max=-1e32;
  m_x_min=m_y_min=m_z_min=1e32;
  unsigned int  count =0 ;

  // ignore 8 first line 
  string buffer;
  for(unsigned int i = 0 ; i < 8 ; i++){
    getline(in,buffer);
  }

  while(in >> x >> y >> z >> Bx >> By >> Bz){
    if(++count%50000==0)
      cout << "\r  - Loading " << count << " values " << flush; 
    vector<double> p = {x,y,z};
    Bx*=tesla;
    By*=tesla;
    Bz*=tesla;
    vector<double> B = {Bx,By,Bz};
    m_field[p]=B;
    if(x<m_x_min)
      m_x_min=x;
    if(x>m_x_max)
      m_x_max=x;  
    if(y<m_y_min)
      m_y_min=y;
    if(y>m_y_max)
      m_y_max=y;  
    if(z<m_z_min)
      m_z_min=z;
    if(z>m_z_max)
      m_z_max=z;  
  }

  m_Rmax=m_x_max;
  cout << "\r  - " << count << " values loaded" << endl; 
  cout << "  - min(" << m_x_min <<";"<< m_y_min <<";" << m_z_min<< ") max(" << m_x_max <<";"<< m_y_max <<";" << m_z_max<< ")" << endl; 
  cout << "  - Rmax = " << m_Rmax << endl;
  in.close();
}
////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadBinary(std::string file){
  ifstream in(file.c_str(),std::ifstream::binary);
  if(!in.is_open()){
    cout << "Error: failed to load samurai field map " << file << endl;
    exit(1);
  }

  cout << "//////// Loading Binary Samurai field map " << file << endl; 
  double x,y,z,Bx,By,Bz;

  m_x_max=m_y_max=m_z_max=-1e32;
  m_x_min=m_y_min=m_z_min=1e32;
  unsigned int  count =0 ;
  while(!in.eof()){

    if(++count%50000==0)
      cout << "\r  - Loading " << count << " values " << flush; 

    in.read((char*)&x,sizeof(x));
    in.read((char*)&y,sizeof(y));
    in.read((char*)&z,sizeof(z));
    in.read((char*)&Bx,sizeof(Bx));
    in.read((char*)&By,sizeof(By));
    in.read((char*)&Bz,sizeof(Bz));
    vector<double> p = {x,y,z};
    Bx*=tesla;
    By*=tesla;
    Bz*=tesla;
    vector<double> B = {Bx,By,Bz};
    m_field[p]=B;
    if(x<m_x_min)
      m_x_min=x;
    if(x>m_x_max)
      m_x_max=x;  
    if(y<m_y_min)
      m_y_min=y;
    if(y>m_y_max)
      m_y_max=y;  
    if(z<m_z_min)
      m_z_min=z;
    if(z>m_z_max)
      m_z_max=z;  
  }
  
  m_Rmax=m_x_max;
  cout << "\r  - " << count << " values loaded" << endl; 
  cout << "  - min(" << m_x_min <<";"<< m_y_min <<";" << m_z_min<< ") max(" << m_x_max <<";"<< m_y_max <<";" << m_z_max<< ")" << endl; 
  cout << "  - Rmax = " << m_Rmax << endl;
  in.close();
}
