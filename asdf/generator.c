{
	ofstream submit_macro("/home/jlee/NPTools/NPTools_IDATEN/asdf/submit.sh");

	for (auto i = 1; i <= 20; i++) 
	{ 
		TString name = Form("/home/jlee/NPTools/NPTools_IDATEN/asdf/gamma%02d.source",i);
		ofstream submit_file(name);

		submit_file << "Isotropic" << endl;
		submit_file << " EnergyLow= "<<0.1*i<<" MeV" << endl;
		submit_file << " EnergyHigh= "<<0.1*i<<" MeV" << endl;
		submit_file << " HalfOpenAngleMin= 0 deg" << endl;
		submit_file << " HalfOpenAngleMax= 180 deg" << endl;
		submit_file << " x0= 0 mm" << endl;
		submit_file << " y0= 0 mm" << endl;
		submit_file << " z0= 0 mm" << endl;
		submit_file << " particle= gamma" << endl;

		submit_macro << "npsimulation -D IDATEN.detector -E "<<name<<" -O IDATEN"<<Form("%02d",i)<<" -B /home/jlee/NPTools/NPTools_IDATEN/asdf/batch.mac &" << endl;
	}
}
