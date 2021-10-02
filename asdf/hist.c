#include "SimTree.h"


void hist(){
	const int nfile=40; const double div=100;
	int Entries[nfile]; int permille[nfile];
	SimTree *tree[nfile];
	TCanvas *c1[nfile];
	TH1D *h1_KE[nfile];
	TH1D *h1_FE[nfile];
	TH1D *h1_TE[nfile];

	TF1 *f1_KE[nfile];
	TF1 *f1_FE[nfile];
	TF1 *f1_TE[nfile];

	TGraphErrors *gr_Keff = new TGraphErrors(); int igrK=0;
	gr_Keff->SetName("gr_Keff");
	gr_Keff->SetTitle("gr_Keff; Energy(keV); efficiency(%)");
	TGraphErrors *gr_Feff = new TGraphErrors(); int igrF=0;
	gr_Feff->SetName("gr_Feff");
	gr_Feff->SetTitle("gr_Feff; Energy(keV); efficiency(%)");
	TGraphErrors *gr_Teff = new TGraphErrors(); int igrT=0;
	gr_Teff->SetName("gr_Teff");
	gr_Teff->SetTitle("gr_Teff; Energy(keV); efficiency(%)");

	TFile *output = new TFile("Hist_Eff.root", "recreate");


	//for(int ifile=5; ifile<=8; ifile++){
	for(int ifile=1; ifile<=nfile; ifile++){
		c1[ifile] = new TCanvas(Form("c_%.0fkeV", div*ifile),Form("c_%.0fkeV", div*ifile),2400, 800);
		h1_KE[ifile] = new TH1D(Form("h1_KE%.0fkeV", div*ifile),Form("h1_KE%.0fkeV;Energy(keV)", div*ifile),4500,0,4500);
		h1_FE[ifile] = new TH1D(Form("h1_FE%.0fkeV", div*ifile),Form("h1_FE%.0fkeV;Energy(keV)", div*ifile),4500,0,4500);
		h1_TE[ifile] = new TH1D(Form("h1_TE%.0fkeV", div*ifile),Form("h1_TE%.0fkeV;Energy(keV)", div*ifile),4500,0,4500);

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
			for(int ibr=0; ibr<tree[ifile]->Tigress->GetMultiplicityGe(); ibr++){
				h1_TE[ifile]->Fill(tree[ifile]->Tigress->GetGeEnergy(ibr)*1000);
			}
		}
		tree[ifile]->~SimTree();
		printf("\n");

		c1[ifile]->Divide(3,1);

		c1[ifile]->cd(1);
		h1_KE[ifile]->Draw();
		f1_KE[ifile] = new TF1(Form("f1_KE%.0fkeV", div*ifile), "gaus(0) + [0]*[3]*(0.5-0.5*erf((x-[1])/[2]/TMath::Sqrt(2.))) + [4]*(x-[1])+[5]", div*ifile-3*10, div*ifile+3*10);
		f1_KE[ifile]->SetParameter(0,8000/ifile);
		f1_KE[ifile]->SetParameter(1,div*ifile);
		f1_KE[ifile]->SetParameter(2,10);
		f1_KE[ifile]->SetParameter(3,0);
		f1_KE[ifile]->SetParameter(4,0);
		f1_KE[ifile]->SetParameter(5,0);
		f1_KE[ifile]->SetParLimits(0,0,1000000);
		f1_KE[ifile]->SetParLimits(1,div*ifile-3*f1_KE[ifile]->GetParameter(2), div*ifile+3*f1_KE[ifile]->GetParameter(2));
		f1_KE[ifile]->SetParLimits(2,0,50);
		f1_KE[ifile]->SetParLimits(3,0,1);
		//f1_KE[ifile]->SetParLimits(4,0,1);
		f1_KE[ifile]->SetParLimits(5,0,1000000);
		f1_KE[ifile]->SetRange(div*ifile-3*f1_KE[ifile]->GetParameter(2), div*ifile+3*f1_KE[ifile]->GetParameter(2));
		h1_KE[ifile]->Fit(f1_KE[ifile], "RQ");
		f1_KE[ifile]->SetRange(div*ifile-3*f1_KE[ifile]->GetParameter(2), div*ifile+3*f1_KE[ifile]->GetParameter(2));
		h1_KE[ifile]->Fit(f1_KE[ifile], "R");

		gr_Keff->SetPoint(igrK++, f1_KE[ifile]->GetParameter(1), sqrt(2*TMath::Pi()) * f1_KE[ifile]->GetParameter(0)*f1_KE[ifile]->GetParameter(2)/10000000*100);


		c1[ifile]->cd(2);
		h1_FE[ifile]->Draw();
		f1_FE[ifile] = new TF1(Form("f1_FE%.0fkeV", div*ifile), "gaus(0) + [0]*[3]*(0.5-0.5*erf((x-[1])/[2]/TMath::Sqrt(2.))) + [4]*(x-[1])+[5]", div*ifile-3*10, div*ifile+3*10);
		f1_FE[ifile]->SetParameter(0,8000/ifile);
		f1_FE[ifile]->SetParameter(1,div*ifile);
		f1_FE[ifile]->SetParameter(2,10);
		f1_FE[ifile]->SetParameter(3,0);
		f1_FE[ifile]->SetParameter(4,0);
		f1_FE[ifile]->SetParameter(5,0);
		f1_FE[ifile]->SetParLimits(0,0,1000000);
		f1_FE[ifile]->SetParLimits(1,div*ifile-3*f1_FE[ifile]->GetParameter(2), div*ifile+3*f1_FE[ifile]->GetParameter(2));
		f1_FE[ifile]->SetParLimits(2,0,50);
		f1_FE[ifile]->SetParLimits(3,0,1);
		//f1_FE[ifile]->SetParLimits(4,0,1);
		f1_FE[ifile]->SetParLimits(5,0,1000000);
		f1_FE[ifile]->SetRange(div*ifile-3*f1_FE[ifile]->GetParameter(2), div*ifile+3*f1_FE[ifile]->GetParameter(2));
		h1_FE[ifile]->Fit(f1_FE[ifile], "RQ");
		f1_FE[ifile]->SetRange(div*ifile-3*f1_FE[ifile]->GetParameter(2), div*ifile+3*f1_FE[ifile]->GetParameter(2));
		h1_FE[ifile]->Fit(f1_FE[ifile], "R");

		gr_Feff->SetPoint(igrF++, f1_FE[ifile]->GetParameter(1), sqrt(2*TMath::Pi()) * f1_FE[ifile]->GetParameter(0)*f1_FE[ifile]->GetParameter(2)/10000000*100);

		c1[ifile]->cd(3);
		h1_TE[ifile]->Draw();
		f1_TE[ifile] = new TF1(Form("f1_TE%.0fkeV", div*ifile), "gaus(0) + [0]*[3]*(0.5-0.5*erf((x-[1])/[2]/TMath::Sqrt(2.))) + [4]*(x-[1])+[5]", div*ifile-3*3, div*ifile+3*3);
		f1_TE[ifile]->SetParameter(0,8000/ifile);
		f1_TE[ifile]->SetParameter(1,div*ifile);
		f1_TE[ifile]->SetParameter(2,3);
		f1_TE[ifile]->SetParameter(3,0);
		f1_TE[ifile]->SetParameter(4,0);
		f1_TE[ifile]->SetParameter(5,0);
		f1_TE[ifile]->SetParLimits(0,0,1000000);
		f1_TE[ifile]->SetParLimits(1,div*ifile-3*f1_TE[ifile]->GetParameter(2), div*ifile+3*f1_TE[ifile]->GetParameter(2));
		f1_TE[ifile]->SetParLimits(2,0,50);
		f1_TE[ifile]->SetParLimits(3,0,1);
		//f1_TE[ifile]->SetParLimits(4,0,1);
		f1_TE[ifile]->SetParLimits(5,0,1000000);
		f1_TE[ifile]->SetRange(div*ifile-3*f1_TE[ifile]->GetParameter(2), div*ifile+3*f1_TE[ifile]->GetParameter(2));
		h1_TE[ifile]->Fit(f1_TE[ifile], "RQ");
		f1_TE[ifile]->SetRange(div*ifile-3*f1_TE[ifile]->GetParameter(2), div*ifile+3*f1_TE[ifile]->GetParameter(2));
		h1_TE[ifile]->Fit(f1_TE[ifile], "R");

		gr_Teff->SetPoint(igrT++, f1_TE[ifile]->GetParameter(1), sqrt(2*TMath::Pi()) * f1_TE[ifile]->GetParameter(0)*f1_TE[ifile]->GetParameter(2)/10000000*100);


		output->cd();	c1[ifile]->Write();

	}

	TCanvas *c2 = new TCanvas("c_eff","c_eff", 1800,600);
	c2->Divide(3,1);
	c2->cd(1);
	gr_Keff->Draw("APL");
	c2->cd(2);
	gr_Feff->Draw("APL");
	c2->cd(3);
	gr_Teff->Draw("APL");

	c2->Write();
	output->Close();

}


