#include "Garfield/GeometrySimple.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewMedium.hh"
		 
double GARFDRIFT(double energy, double t_start, G4ThreeVector start, G4ThreeVector delta_pos, double R, int Pad_start, double ID,
		 double ScoLength, double SegLength, double event, double raw_energy_check, string shape) {

  ofstream file;
  /*        
  file.open("garf_E_check.dat", std::ios::app);
  file << raw_energy_check << "\t" << R << endl;
  file.close();
  */
  double rWire;
  if(shape == "Cylindrical") rWire = 1.2;
  if(shape == "Long_Chamber") rWire = -2.5;
  const double rTube = 5.;
  const double lTube = 25.19;
  //TH1F *hist = new TH1F("","",1000,0,10000); // 10 ns bins
  //For R6 = 2, R7 = 1
  //const double vWire = -375.; //For V_total = 750V
  //const double vWire = -645.3312; //For V_total = 1200V //SHOULD BE 583.741 !?
  //const double vWire = -535.0959; //For V_total = 1100V
  //const double vWire = -875.6116; //For V_total = 1800V 
  //const double vWire =  -681.0312; //For V_total = 1400V
  //const double vWire = -730.; 
  const double vWire = -1000.;
  //const double vWire = -584.7139; // ND run 4230 (2008)
  //const double vWire = -603.6855; // ND run 4279 (2008), HV = 1241, R1 = 1, R2 = 2 
  const double vTube = 0.;
  double x_start, y_start, z_start;
  //double production_bias = 0.01;
  double production_bias = 0.001;
  double t_end, x_end, y_end, z_end;
  int status;
  double Pad_end;
  //double w_value = 35.; //for HeCO2 --guess
  //double w_value = 26.31; //for argon
  double w_value = 41.1;// for He
  energy = energy*production_bias;
  //int electrons  = (int)(energy/w_value*production_bias);
  int electrons = (int)(energy/w_value);
  double energy_excess = 0.;
  int sector;
  vector<double> pad_vec, time_vec; 
  
  Garfield::MediumMagboltz* gas = new Garfield::MediumMagboltz();
  Garfield::ViewMedium* mediumView = new Garfield::ViewMedium();
  Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();
  Garfield::SolidTube* tube = new Garfield::SolidTube(0, 0, 0, 0, rTube, lTube/2);
  Garfield::SolidBox* box = new Garfield::SolidBox(0, 0, 0, 3.47, 5.0, lTube/2); 
  Garfield::ComponentAnalyticField* cmp = new Garfield::ComponentAnalyticField();
  Garfield::Sensor* sensor = new Garfield::Sensor();
  Garfield::AvalancheMC* drift = new Garfield::AvalancheMC();
  
  //cout << "\n" << electrons << "\n" << endl;
  
  energy_excess = (energy - electrons*w_value)/production_bias;
  
  //  if(R>rWire) {

    //    energy_excess = (energy/w_value*production_bias - electrons)*w_value/production_bias;

    //    energy_excess = (energy/w_value - electrons)*w_value;
  
  if(R>rWire && electrons>0) {
    
      //gas->LoadGasFile("he_90_co2_10_200mbar.gas");
      //gas->LoadGasFile("ar_90_ch4_10.gas"); //atomic fraction in this case (P10).
    //gas->LoadGasFile("he_90_co2_10_500mbar.gas"); //mass fraction in this case
      //gas->LoadGasFile("he_90_co2_10_1bar.gas");

    gas->LoadGasFile("he_90_co2_10_100mbar.gas");
    //gas->LoadGasFile("he_90_co2_10_700mbar.gas");
    
    gas->Initialise(true);
      
    mediumView->SetMedium(gas);
    if(shape == "Cylindrical") geo->AddSolid(tube, gas);
    if(shape == "Long_Chamber") geo->AddSolid(box, gas);
    cmp->SetGeometry(geo);
    if(shape == "Cylindrical") {
      cmp->AddWire(0, 0, 2*rWire, vWire, "c");
      cmp->AddTube(rTube, vTube, 0, "a");
    }
    if(shape == "Long_Chamber") {
      cmp->AddPlaneY(-2.5,-500.,"c"); //CHANGE CATHODE AND ANODE VOLTAGE HERE (FOR LONG CHAMBER)
      cmp->AddPlaneY(2.5,0.,"a");
    }
    sensor->AddComponent(cmp);
    drift->SetSensor(sensor);

      //drift->DisableAttachment();
    
    for(int e=0;e<electrons;e++) {
      
      double randomize = (double)std::rand() / (double)RAND_MAX;
	
      x_start = start.x() + delta_pos.x()*randomize;
      y_start = start.y() + delta_pos.y()*randomize;
      z_start = start.z() + delta_pos.z()*randomize;
      /*
	file.open("coord_check.dat", std::ios::app);
	file <<  start.x() << "\t" << start.y() << "\t" << start.z() << x_start  endl;	  
	file.close();
      */
      drift->DriftElectron(x_start, y_start, z_start, t_start);
      
      int np = drift->GetNumberOfDriftLinePoints();
      //drift->GetEndPoint(x_end, y_end, z_end, t_end, status);
      drift->GetDriftLinePoint(np-1, x_end, y_end, z_end, t_end); //np-1 is last DriftLineEndPoint
      //int np = drift->GetNuberOfElectronEndPoints();
      //drift->GetElectronEndPoint(np-1, x_end, y_end, z_end, t_end); 
      
      Pad_end = (int)((z_end + ScoLength / 2.) / SegLength ) + 1; //new Pad number 

      if(shape=="Cylindrical") {
      
	if(x_end > 0 and y_end > 0) {
	  if(abs(x_end) > abs(y_end)) sector = 0;
	  if(abs(x_end) < abs(y_end)) sector = 1;
	}
	if(x_end < 0 and y_end > 0) {
	  if(abs(x_end) < abs(y_end)) sector = 2;
	  if(abs(x_end) > abs(y_end)) sector = 3;
	}
	if(x_end < 0 and y_end < 0) {
	  if(abs(x_end) > abs(y_end)) sector = 4;
	  if(abs(x_end) < abs(y_end)) sector = 5;
	}
	if(x_end > 0 and y_end < 0) {
	  if(abs(x_end) < abs(y_end)) sector = 6;
	  if(abs(x_end) > abs(y_end)) sector = 7;
	}
	
      }

      if(shape=="Long_Chamber") sector = 0;
      
      //if(pow(pow(x_end,2)+pow(y_end,2),0.5)>4.9) { //reached the anode
      file.open("signal.dat", std::ios::app);
      if(shape == "Cylindrical") file << event << "\t" << ID << "\t" << sector << "\t" <<  Pad_end << "\t" << t_end << "\t" << pow(pow(x_start,2)+pow(y_start,2),0.5) << "\t" << pow(pow(x_end,2)+pow(y_end,2),0.5) << endl;
      if(shape == "Long_Chamber") file << event << "\t" << ID << "\t" << sector << "\t" <<  Pad_end << "\t" << z_end <<  "\t" << t_end << "\t" << y_start << "\t" << y_end << endl;
      //file << event << "\t" << ID << "\t" << sector << "\t" <<  Pad_end << "\t" << t_end
      file.close();
      //}
      
    }
  }
  //}

delete gas; delete mediumView; delete geo; delete tube; delete cmp; delete sensor; delete drift; //otherwise these are held in memory and cause a killed: 9 crash.

  return energy_excess;

}
