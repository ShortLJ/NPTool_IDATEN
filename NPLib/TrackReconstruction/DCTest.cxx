#include"NPDCReconstruction.h"
#include"NPSystemOfUnits.h"
using namespace NPUNITS;
R__LOAD_LIBRARY(libNPTrackReconstruction.dylib)
  //
  void DrawTrack(double a, double b, int style);
  //
  void DrawCircle(double x, double z, double r);
  //
  void DCTest(){
    // load the lib
    auto c = new TCanvas();
    c->Divide(2,1);
    TRandom3 rand;
    // The ronconstruction utility object
    NPL::DCReconstruction dcr;
    vector<double> x,z,r;
    // maximum drift length
    double maxR = 10; 
    double stepZ = 50;

    // a drawing artefact
    // give axis and background to ResolvePlane
    auto bck = new TH2F("bck","bck",10,-150,150,10,-150,150);

    // infinite loop
    while(true){
      // Check 2D minimisation
      c->cd(1);
      x.clear();z.clear();r.clear();
      // create a random trajectory in 2D //
      // a and b are the parameter of the 2D track
      double a = rand.Uniform(-10,10);
      double b = rand.Uniform(-100,100); 
      DrawTrack(a,b,2); 
      // Generate n points from the line
      int n = (int) rand.Uniform(3,20);
      for(unsigned int i = 0 ; i < n ; i++){
        double Z = i*stepZ-100;
        z.push_back(Z); 
        double R = rand.Uniform(-maxR,maxR);
        r.push_back(R);
        double X = (Z-b)/a+R;
        x.push_back(X);
        DrawCircle(X,Z,R);
      }
/*      // add second track for pile up
      a = rand.Uniform(-10,10);
      b = rand.Uniform(-100,100); 
      DrawTrack(a,b,3); 
      // Generate n points from the line
      n = (int) rand.Uniform(3,20);
      for(unsigned int i = 0 ; i < n ; i++){
        double Z = i*stepZ-100;
        z.push_back(Z); 
        double R = rand.Uniform(-maxR,maxR);
        r.push_back(R);
        double X = (Z-b)/a+R;
        x.push_back(X);
        DrawCircle(X,Z,R);
      }
*/
      // add noise
      // Generate n points from the line
/*      n = (int) rand.Uniform(1,3);
      for(unsigned int i = 0 ; i < n ; i++){
        double Z = i*stepZ-100;
        z.push_back(Z); 
        double R = rand.Uniform(-maxR,maxR);
        r.push_back(R);
        double X = rand.Uniform(-100,100);
        x.push_back(X);
        DrawCircle(X,Z,R);
      }
*/
 
      double x0,x100,af,bf;
      double res = dcr.BuildTrack2D(x,z,r,x0,x100,af,bf);
      cout << res << endl;
      DrawTrack(af,bf,1);

      // Check 2D Minimize Func
      c->cd(2);
      TGraph* ga = dcr.Scan(af,bf,0,af*0.5,af*2);  
      TGraph* gb = dcr.Scan(af,bf,1,bf*0.5,bf*2);  
      ga->SetLineColor(kAzure+7);
      ga->Draw("ac");
      ga->GetYaxis()->SetRangeUser(0,1000);

      gb->Draw("c");
      
      c->Update();
      gPad->WaitPrimitive();
    }
  }

// draw a track based on z=ax+b
void DrawTrack(double a, double b, int style){
  double x[2]={-100,100};
  double z[2]={a*-100+b,a*100+b};
  auto l = new TGraph(2,x,z);
  l->SetLineColor(kBlack);
  if(style>1)
    l->SetLineWidth(2);
  l->SetLineStyle(style);
  if(style==2)
    l->Draw("al");
  else
    l->Draw("l");
}
// draw a circle
void DrawCircle(double x, double z, double r){
  auto e = new TEllipse(x,z,r,r);
  e->SetLineColor(kBlack);
  e->Draw("same"); 
}
