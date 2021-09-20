#ifndef SamuraiFieldMap_h
#define SamuraiFieldMap_h
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

#include<string>
#include<vector>
#include<map>
#include"TObject.h"
#include"TGraph.h"
#include"TVector3.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"

#include "NPParticle.h"
class SamuraiFieldMap{

  public:
    SamuraiFieldMap();
    SamuraiFieldMap(std::string file);
    ~SamuraiFieldMap(){};
  
  public: // Map reading
    void LoadMap(double angle, std::string file, unsigned int bin);

  private:
    void LoadAscii(std::string file);
    void LoadBinary(std::string file);

  private:
    // map[Pos]=B;
    std::map<std::vector<double>,std::vector<double>> m_field;
    double m_x_max,m_y_max,m_z_max,m_x_min,m_y_min,m_z_min;
    int m_bin;
    double m_angle;
    double m_Rmax ;

  public: // getting the field at a point in space
    // return B at an existing point
    std::vector<double> GetB(std::vector<double>& pos); 
    inline std::vector<double> GetB(double x,double y ,double z){
      std::vector<double> pos = {x,y,z};
      return GetB(pos);
    };
    
    // interpolate B witin volume (0 outside volume)
    std::vector<double> InterpolateB(const std::vector<double>& pos);
    // interpolate B witin volume (0 outside volume)
    inline std::vector<double> InterpolateB(const TVector3& pos){
      std::vector<double> p={(double)pos.X(),(double)pos.Y(),(double)pos.Z()};
      return InterpolateB(p);
    };
 
  public: // Propagation of a particule in the field
    // return a 3D track of the particle in the field
    std::vector< TVector3 > Propagate(double Brho, TVector3 pos, TVector3 dir,bool store=true);
    void func(NPL::Particle& N, TVector3 pos, TVector3 imp, TVector3& new_pos, TVector3& new_dir);
  private:
    double m_fdc2angle;
    double m_fdc2R;
  public:
    void SetFDC2Angle(double angle){m_fdc2angle=angle;};
    void SetFDC2R(double R){m_fdc2R=R;};
    TVector3 PropagateToFDC2(TVector3 pos, TVector3 dir);

  public:
    TGraph* BrhoScan(double min,double max,double step);
    double  FindBrho(TVector3 p_fdc0,TVector3 d_fdc0,TVector3 p_fdc2,TVector3 d_fdc2);

 
  private:
    TGraph* m_BrhoScan;
    ROOT::Math::Minimizer* m_min;
    ROOT::Math::Functor    m_func;
    double Delta(const double* parameter);
    TVector3 m_FitPosFDC0,m_FitDirFDC0,m_FitPosFDC2,m_FitDirFDC2;
    

    //
    ClassDef(SamuraiFieldMap,1);
};

#endif
