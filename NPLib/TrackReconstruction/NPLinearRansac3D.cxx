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

#include"NPLinearRansac3D.h"
using namespace std;
using namespace NPL;
////////////////////////////////////////////////////////////////////////////////
void LinearCluster3D::LinearFit(){

  unsigned int sizeI = m_index.size();   
  double min_distance = 1e20;
  unsigned int min_i = 0 ; unsigned int min_j = 0; 

  double meanDX=0;
  double meanDY=0;
  double meanDZ=0;
  double totalW=0;
  double meanX=0;
  double meanY=0;
  double meanZ=0;
  double w;

  for(unsigned int i = 0 ; i < sizeI ; i++){
    TVector3 Pi((*m_X)[m_index[i]],(*m_Y)[m_index[i]],(*m_Z)[m_index[i]]);
    for(unsigned int j = i+1 ; j < sizeI ; j++){
      double dij=0;
      TVector3 Pj((*m_X)[m_index[j]],(*m_Y)[m_index[j]],(*m_Z)[m_index[j]]);
      // compute the average distance of all point to this line
      for(unsigned int p = 0 ; p < sizeI ; p++){
        TVector3 Pp((*m_X)[m_index[p]],(*m_Y)[m_index[p]],(*m_Z)[m_index[p]]);
        dij+=MinimumDistancePointLine(Pi,Pj,Pp);
      }
      
      w = 1./dij; 
      totalW+=w;

      meanX+=w*((*m_X)[m_index[i]]);
      meanY+=w*((*m_Y)[m_index[i]]);
      meanZ+=w*((*m_Z)[m_index[i]]);

      meanDX+=w*((*m_X)[m_index[i]]-(*m_X)[m_index[j]]);
      meanDY+=w*((*m_Y)[m_index[i]]-(*m_Y)[m_index[j]]);
      meanDZ+=w*((*m_Z)[m_index[i]]-(*m_Z)[m_index[j]]);
    }  
  }

  meanDX/=totalW;
  meanDY/=totalW;
  meanDZ/=totalW;
  meanX/=totalW;
  meanY/=totalW;
  meanZ/=totalW;

  m_P0  = TVector3(meanX,meanY,meanZ);
  m_Dir = (TVector3(meanDX,meanDY,meanDZ)).Unit();

};

////////////////////////////////////////////////////////////////////////////////
  vector<LinearCluster3D> LinearRansac3D::TreatEvent(vector<double>& X, vector<double>&Y, vector<double>&Z){
      cluster_id.clear();
      clusters.clear();
      unsigned int sizeX = X.size(); 
      cluster_id.resize(sizeX,0);
      m_cluster.clear();
      m_assigned.clear();

      if(sizeX<m_min_count)
        return clusters;

      unsigned int p1,p2;
      double d;
      TVector3 Vp1,Vp2,D,P;
      m_iteration_d=0;
      m_iteration=0;

      while(m_iteration++ < m_max_iteration && m_assigned.size()<sizeX){
        LinearCluster3D cluster(&X,&Y,&Z);
        // take 2 distant point randomly that has not been match before
        d=0 ; m_iteration_d=0; 
        while(d<3*m_match_distance  && m_iteration_d++<m_max_iteration ){
          p1 = rand.Integer(sizeX); 
          p2 = rand.Integer(sizeX); 

          Vp1.SetXYZ(X[p1],Y[p1],Z[p1]);
          Vp2.SetXYZ(X[p2],Y[p2],Z[p2]);
          D=Vp1-Vp2;
          d  = D.Mag();
          if(d>m_max_distance)
            d=0;
        }

        // loop over all points
        for(unsigned int i = 0 ; i < sizeX ; i++){
          P.SetXYZ(X[i],Y[i],Z[i]);
          if(MinimumDistancePointLine(Vp1,Vp2,P) < m_match_distance){
            cluster.AddIndex(i);
            m_assigned.insert(i);
          }
        }
        // insert the newly formed cluster
        if(cluster.size()>m_min_count){
          m_cluster.insert(cluster);
        }

      }//while
      
      // loop over the cluster starting with the biggest
      unsigned int current_cluster=0;
      unsigned int index;
      for(auto it = m_cluster.begin() ; it!=m_cluster.end() ; ++it){
        current_cluster++;
        //  cout << current_cluster << endl;
        unsigned int sizeC = (*it).size();
        unsigned int cluster_size = 0;
        for(unsigned int i = 0 ; i < sizeC ; i++){
          // Assigned cluster id to identified points
          unsigned int index = (*it).GetIndex(i);
          if(!cluster_id[index]){
            cluster_id[index]=current_cluster;
            cluster_size++;
          }
        }
        if(cluster_size<m_min_count){
          for(unsigned int i = 0 ; i < sizeC ; i++){
            unsigned int index = (*it).GetIndex(i);
            // remove the assigned point 
            if(cluster_id[index]==current_cluster){
              cluster_id[index]=0;
            }
          }
        }
        else{
          clusters.push_back(*it);
        }
      }
      return clusters;
    }

