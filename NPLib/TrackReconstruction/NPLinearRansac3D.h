#ifndef NPLinearLinearRansac3D_h
#define NPLinearLinearRansac3D_h
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta  contact address: matta@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : April 2021                                              *
 * Last update     : April 2021                                              *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Find linear tracks using the ransac method                               *
 *****************************************************************************/
//std
#include<vector>
#include<set>
#include<map>
//root
#include "TVector3.h"
#include "TRandom2.h"
// nptool
#include "NPTrackingUtility.h"
namespace NPL{
  /////////////////////////////////
  class LinearCluster3D{
    public:
      LinearCluster3D(){};
      LinearCluster3D(std::vector<double>* X, std::vector<double>* Y, std::vector<double>* Z){
        m_X=X; m_Y=Y; m_Z=Z;
      }
      ~LinearCluster3D(){};

    public:
      void AddIndex(const unsigned int& i) {m_index.push_back(i);};
      unsigned int size() const {return m_index.size();};
      unsigned int GetIndex (const unsigned int& i) const {return m_index[i];};

      // Fit the cluster
      void LinearFit();

    private:
      std::vector<double>* m_X;//!
      std::vector<double>* m_Y;//!
      std::vector<double>* m_Z;//!
      std::vector<unsigned int> m_index;//!
      TVector3 m_P0;// a point on the line
      TVector3 m_Dir;// direction of the line

    public:
      std::vector<unsigned int> GetIndex() const {return m_index;};
      TVector3 GetP0()const {return m_P0;};
      TVector3 GetDir()const {return m_Dir;};

    public: //operator
      // reversed logic so the cluster are ordered from largest to smallest
      bool operator<(const LinearCluster3D& B){
        return this->size()>B.size();
      }

      friend bool operator>(const LinearCluster3D& A, const LinearCluster3D& B){
        return (A.size()<B.size());
      } 
      friend bool operator<(const LinearCluster3D& A, const LinearCluster3D& B){
        return (A.size()>B.size());
      } 
  };
  /////////////////////////////////
  class LinearRansac3D{

    public:

      LinearRansac3D(){
        SetParameters(100,5,40,10);// reasonnable parameters for typical detector
      }
      LinearRansac3D(unsigned int max_iteration,double match_distance,double max_distance,unsigned int min_count){
       SetParameters(max_iteration,match_distance,max_distance,min_count);
      };

      ~LinearRansac3D(){};

    private:
      unsigned int m_max_iteration;// max iteration before giving up
      double m_match_distance;// max point to Line distance to belong to the cluster
      double m_max_distance;// max distance a test couple is allowed
      unsigned int m_min_count; // minimum point a cluster to be valid
      unsigned int m_iteration; //number of iteration within a TreatEvent
      unsigned int m_iteration_d; //number of iteration to find a valid pair
      TRandom2 rand; // Random generator
      std::multiset<LinearCluster3D> m_cluster;// all of the build cluster
      std::set<int>                  m_assigned;// list of unique point that are in any cluster
      std::vector<int> cluster_id;
      std::vector<LinearCluster3D> clusters;


    public:
      unsigned int Iterations(){return m_iteration;};
      void SetParameters(unsigned int max_iteration,double match_distance,double max_distance,unsigned int min_count){
        m_max_iteration = max_iteration;
        m_max_distance  = max_distance;
        m_match_distance  = match_distance;
        m_min_count  = min_count;
      };

    public:
      std::vector<LinearCluster3D> TreatEvent(std::vector<double>& X, std::vector<double>&Y, std::vector<double>&Z);
  };


}

#endif 
