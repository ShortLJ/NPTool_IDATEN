//----------! Developed by C. Santamaria / CEA Saclay !----------
//----------!      Version date :: 2014/11/12         !----------
// Tracking Functions
#ifndef TRACKING_H
#define TRACKING_H

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "TGraph.h"
#include "TVector3.h"

using namespace std;

namespace NPL{
class Tracking {
    public:
		Tracking();
		virtual ~Tracking();
        /* double conv_fit(double *x, double *p); */
        int Hough_modified(vector<double> *x,vector<double> *y,vector<double> *q, vector<double> *x_out,vector<double> *y_out, vector<double> *q_out, vector<int> *ringbool);
        double FitFunction(double *x, double *p);
        void FindStart(double pStart[4], double chi[2],  int fitStatus[2],TGraph *grxz, TGraph *gryz);
        double distance2(double x,double y,double z, double *p);
        void Hough_3D(vector<double> *x,vector<double> *y,vector<double> *z,vector<double> *q,vector<double> *x_out,vector<double> *y_out,vector<double> *z_out,vector<double> *q_out);
        void vertex(double *p, double *pp, double &xv,double &yv,double &zv, double &min_dist, double &Theta_tr1, double &Theta_tr2, double &Phi1, double &Phi2,double *VectorTrack1, double *VectorTrack2);
        void ParFor_Vertex(double *a, double *b, double *parFit);

	ClassDef(Tracking,1);

};
}
#endif
