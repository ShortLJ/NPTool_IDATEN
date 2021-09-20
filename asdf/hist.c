#include "SimTree.h"


void hist(){
	const int nfile=18; const double div=100;
	int Entries[nfile]; int permille[nfile];
	SimTree *tree[nfile];
	TCanvas *c1[nfile];
	TH1D *h1_KE[nfile];
	TH1D *h1_FE[nfile];

	TF1 *f1_KE[nfile];
	TF1 *f1_FE[nfile];

	TGraphErrors *gr_Keff = new TGraphErrors(); int igrK=0;
	gr_Keff->SetName("gr_Keff");
	gr_Keff->SetTitle("gr_Keff; Energy(keV); efficiency(%)");
	TGraphErrors *gr_Feff = new TGraphErrors(); int igrF=0;
	gr_Feff->SetName("gr_Feff");
	gr_Feff->SetTitle("gr_Feff; Energy(keV); efficiency(%)");


	//for(int ifile=5; ifile<=8; ifile++){
	for(int ifile=1; ifile<=nfile; ifile++){
		c1[ifile] = new TCanvas(Form("c_%.0fkeV", div*ifile),Form("c_%.0fkeV", div*ifile),2400,1000);
		h1_KE[ifile] = new TH1D(Form("h1_KE%.0fkeV", div*ifile),Form("h1_KE%.0fkeV;Energy(keV)", div*ifile),2000,0,2000);
		h1_FE[ifile] = new TH1D(Form("h1_FE%.0fkeV", div*ifile),Form("h1_FE%.0fkeV;Energy(keV)", div*ifile),2000,0,2000);

		tree[ifile] = new SimTree(0,ifile);
		Entries[ifile]=tree[ifile]->fChain->GetEntries();
		permille[ifile]=Entries[ifile]/1000;

		printf("ifile %d\n",ifile);
		for(int ient=0; ient<Entries[ifile]; ient++){
			if(ient%permille[ifile]) printf("%dient of %dent, %.1f%%\r", ient, Entries[ifile], (float)ient/Entries[ifile]*100);
			tree[ifile]->GetEntry(ient);
			for(int ibr=0; ibr<tree[ifile]->Khala->GetKhalaLaBr3EMult(); ibr++){
				h1_KE[ifile]->Fill(tree[ifile]->Khala->GetKhalaLaBr3EEnergy(ibr)*1000);
			}
			for(int ibr=0; ibr<tree[ifile]->Fatima->GetFatimaLaBr3EMult(); ibr++){
				h1_FE[ifile]->Fill(tree[ifile]->Fatima->GetFatimaLaBr3EEnergy(ibr)*1000);
			}
		}
		tree[ifile]->~SimTree();
		printf("\n");

		c1[ifile]->Divide(2,1);

		c1[ifile]->cd(1);
		h1_KE[ifile]->Draw();
		f1_KE[ifile] = new TF1(Form("f1_KE%.0fkeV", div*ifile), "gaus(0)+[3]*(x-[1])+[4]", div*ifile-3*10, div*ifile+3*10);
		f1_KE[ifile]->SetParameter(0,8000/ifile);
		f1_KE[ifile]->SetParameter(1,div*ifile);
		f1_KE[ifile]->SetParameter(2,10);
		f1_KE[ifile]->SetParameter(3,0);
		f1_KE[ifile]->SetParameter(4,0);
		f1_KE[ifile]->SetParLimits(0,0,1000000);
		f1_KE[ifile]->SetParLimits(1,div*ifile-3*f1_KE[ifile]->GetParameter(2), div*ifile+3*f1_KE[ifile]->GetParameter(2));
		f1_KE[ifile]->SetParLimits(2,0,50);
		f1_KE[ifile]->SetRange(div*ifile-3*f1_KE[ifile]->GetParameter(2), div*ifile+3*f1_KE[ifile]->GetParameter(2));
		h1_KE[ifile]->Fit(f1_KE[ifile], "RQ");
		f1_KE[ifile]->SetRange(div*ifile-3*f1_KE[ifile]->GetParameter(2), div*ifile+3*f1_KE[ifile]->GetParameter(2));
		h1_KE[ifile]->Fit(f1_KE[ifile], "R");

		gr_Keff->SetPoint(igrK++, f1_KE[ifile]->GetParameter(1), sqrt(2*TMath::Pi()) * f1_KE[ifile]->GetParameter(0)*f1_KE[ifile]->GetParameter(2)/10000000*100);


		c1[ifile]->cd(2);
		h1_FE[ifile]->Draw();
		f1_FE[ifile] = new TF1(Form("f1_FE%.0fkeV", div*ifile), "gaus(0)+[3]*(x-[1])+[4]", div*ifile-3*10, div*ifile+3*10);
		f1_FE[ifile]->SetParameter(0,8000/ifile);
		f1_FE[ifile]->SetParameter(1,div*ifile);
		f1_FE[ifile]->SetParameter(2,10);
		f1_FE[ifile]->SetParameter(3,0);
		f1_FE[ifile]->SetParameter(4,0);
		f1_FE[ifile]->SetParLimits(0,0,1000000);
		f1_FE[ifile]->SetParLimits(1,div*ifile-3*f1_FE[ifile]->GetParameter(2), div*ifile+3*f1_FE[ifile]->GetParameter(2));
		f1_FE[ifile]->SetParLimits(2,0,50);
		f1_FE[ifile]->SetRange(div*ifile-3*f1_FE[ifile]->GetParameter(2), div*ifile+3*f1_FE[ifile]->GetParameter(2));
		h1_FE[ifile]->Fit(f1_FE[ifile], "RQ");
		f1_FE[ifile]->SetRange(div*ifile-3*f1_FE[ifile]->GetParameter(2), div*ifile+3*f1_FE[ifile]->GetParameter(2));
		h1_FE[ifile]->Fit(f1_FE[ifile], "R");

		gr_Feff->SetPoint(igrF++, f1_FE[ifile]->GetParameter(1), sqrt(2*TMath::Pi()) * f1_FE[ifile]->GetParameter(0)*f1_FE[ifile]->GetParameter(2)/10000000*100);


	}

	TCanvas *c2 = new TCanvas("c_eff","c_eff", 1800,600);
	c2->Divide(2,1);
	c2->cd(1);
	gr_Keff->Draw("APL");
	c2->cd(2);
	gr_Feff->Draw("APL");



}


