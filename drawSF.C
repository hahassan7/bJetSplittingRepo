
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <iostream>

using std::cout;
using std::endl;

void drawSF(){
	
	const int nBins = 1;
	const string workpts[nBins] = {"0.8"};
	const string legendTitleFull[nBins] = {"Fulljet CSV>0.8"};
	
	const int ptBins = 2;
	const int colorArr[nBins*ptBins] = {4,4};
	const string legendTitle[nBins*ptBins] = {"120 < p_{T} < 160 GeV","160 < p_{T} < 400 GeV"};
	
	TFile *ff[nBins][ptBins];
	TFile *ff_light[nBins][ptBins];
	TH1D *sfs[nBins][ptBins];
	TH1D *dataEffs[nBins][ptBins];
	TH1D *mcEffs[nBins][ptBins];
	
	TH1D *dataPurs[nBins][ptBins];
	TH1D *mcPurs[nBins][ptBins];
	TH1D *sfsLight[nBins][ptBins];
	
	TH1D *dataPurs_sub[nBins][ptBins];
	TH1D *mcPurs_sub[nBins][ptBins];
	
	TFile *ff_full[nBins][ptBins];
	TH1D *sfs_full[nBins][ptBins];
	TH1D *dataEff_full[nBins][ptBins];
	TH1D *mcEff_full[nBins][ptBins];
	for(int i=0; i<nBins; i++){
		for(int j=0; j<ptBins; j++){
			//if(i==1) ff[i][j] = new TFile("sf_subjets_zg__muFiltered_csv0.8_sv1p5_dataJPCalibOnMC.root");
			ff[i][j] = new TFile(Form("efficiencies/sf_subjets_zg__muFiltered_csv%s_sv1p2_modDataPurErr_fullJetPtBin%d.root",workpts[i].c_str(),j+1));
		//ff_light[i] = new TFile(Form("sf_lightJets_subjets_csv%s.root",workpts[i].c_str()));
			ff_full[i][j] = new TFile(Form("efficiencies/sf_fullJets_csv0.8.root"));//,workpts[i].c_str()));
			
			sfs[i][j] = (TH1D*)ff[i][j]->Get("scaleFactors_btagSF");
			dataEffs[i][j] = (TH1D*)ff[i][j]->Get("dataEff");
			dataEffs[i][j]->SetMarkerStyle(24); dataEffs[i][j]->SetLineColor(1); dataEffs[i][j]->SetMarkerColor(1);
			mcEffs[i][j] = (TH1D*)ff[i][j]->Get("mcEff");
			sfs_full[i][j] = (TH1D*)ff_full[i][j]->Get("scaleFactors_btagSF");
			dataEff_full[i][j] = (TH1D*)ff_full[i][j]->Get("dataEff");
			dataEff_full[i][j]->SetMarkerStyle(24); dataEff_full[i][j]->SetLineColor(1); dataEff_full[i][j]->SetMarkerColor(1);
			mcEff_full[i][j] = (TH1D*)ff_full[i][j]->Get("mcEff");
			
			mcPurs[i][j] = (TH1D*)ff_full[i][j]->Get("mcPur");
			sfsLight[i][j] = (TH1D*)ff_full[i][j]->Get("sfLight");
			dataPurs[i][j] = (TH1D*)ff_full[i][j]->Get("dataPur");
			
			mcPurs_sub[i][j] = (TH1D*)ff[i][j]->Get("mcPur");
			dataPurs_sub[i][j] = (TH1D*)ff[i][j]->Get("dataPur");
		}
	}

	
	TLegend *ll = new TLegend(0.5,0.5,0.9,0.9);
	TLegend *ll2 = new TLegend(0.1,0.1,0.5,0.5);
	TLatex *llat[ptBins];
	TCanvas *cc = new TCanvas("cc","",1200,800);
	cc->Divide(ptBins,2);
	for(int i=0; i<ptBins; i++){
		cc->cd(i+1); //1,2,3
		mcEffs[0][i]->SetMarkerColor(colorArr[i]);
		mcEffs[0][i]->SetLineColor(colorArr[i]);
		mcEff_full[0][i]->SetMarkerColor(colorArr[i]);
		mcEff_full[0][i]->SetLineColor(colorArr[i]);
		dataEffs[0][i]->SetMarkerStyle(24);
		dataEff_full[0][i]->SetMarkerStyle(25);
		mcEff_full[0][i]->SetMarkerStyle(21);
		mcEffs[0][i]->SetMarkerStyle(20);
		dataEffs[0][i]->SetYTitle("Tagging Efficiency");
		//dataEffs[0][i]->SetXTitle("Subjet p_{T}");
		dataEffs[0][i]->SetXTitle("Z_{g}");
		dataEffs[0][i]->Draw("e");
		dataEffs[0][i]->GetXaxis()->SetNdivisions(505);
		dataEffs[0][i]->GetYaxis()->SetRangeUser(0.,1.0);
		mcEffs[0][i]->Draw("same");
		//mcEff_full[i]->Draw("same");
		//dataEff_full[i]->Draw("e,same");
		llat[i] = new TLatex(0.1,0.4,legendTitle[i].c_str());
		llat[i]->Draw("same");
		if(i==1){
			ll->AddEntry(mcEffs[0][i], "MC Subjets");
			ll->AddEntry(dataEffs[0][i], "Data Subjets");
			//ll->AddEntry(mcEff_full[i], "MC Full jets");
			//ll->AddEntry(dataEff_full[i], "Data Full jets");
			ll->Draw();
		} 
		
		cc->cd(i+1+ptBins); //4,5,6
		sfs[0][i]->SetMarkerColor(colorArr[i]);
		sfs[0][i]->SetLineColor(colorArr[i]);
		sfs[0][i]->SetMarkerSize(1);
		sfs[0][i]->SetMarkerStyle(20);
		sfs[0][i]->SetXTitle("Z_{g}");
		//sfs[0][i]->SetXTitle("Subjet p_{T}");
		sfs[0][i]->SetYTitle("Tag Eff. Scale Factor (Data/MC)");
		//sfs[i]->Draw();
		sfs[0][i]->Draw();
		sfs[0][i]->GetXaxis()->SetNdivisions(505);
		sfs[0][i]->GetYaxis()->SetRangeUser(0.5,1.5);
		//ll->AddEntry(sfs[i],legendTitle[i].c_str(),"lp");
		
		sfs_full[0][i]->SetMarkerStyle(25);
		sfs_full[0][i]->SetMarkerSize(1);
		sfs_full[0][i]->SetMarkerColor(colorArr[i]);
		sfs_full[0][i]->SetLineColor(colorArr[i]);
		sfs_full[0][i]->GetYaxis()->SetRangeUser(0.5,1.5);
		//sfs_full[i]->Draw("same");
		if(i==nBins-1){
			ll2->AddEntry(sfs[0][i],"Subjet SF","lp");
			//ll2->AddEntry(sfs_full[0][i],"Full Jet SF","lp");
			//ll2->Draw();
		}
		
	}
	//ll->Draw("Same");
	//ll2->Draw("same");
	
	TCanvas *cc4 = new TCanvas("cc4","",1200,600);
	cc4->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cc4->cd(i+1);
		for(int ibin=1; ibin<=dataPurs_sub[0][i]->GetNbinsX(); ibin++){
			if(1.-dataPurs_sub[0][i]->GetBinContent(ibin) < dataPurs_sub[0][i]->GetBinError(ibin) ) dataPurs_sub[0][i]->SetBinError(ibin, (1.-dataPurs_sub[0][i]->GetBinContent(ibin)));
			if(dataPurs_sub[0][i]->GetBinError(ibin) != dataPurs_sub[0][i]->GetBinError(ibin)) dataPurs_sub[0][i]->SetBinError(ibin, 0.001);
		}
		
		dataPurs_sub[0][i]->SetXTitle("zG");
		dataPurs_sub[0][i]->SetYTitle("Tagging Purity");
		dataPurs[0][i]->SetLineColor(1);
		dataPurs_sub[0][i]->SetLineColor(1);
		dataPurs_sub[0][i]->SetMarkerStyle(24);
		dataPurs[0][i]->SetMarkerStyle(25);
		mcPurs[0][i]->SetMarkerStyle(21);
		mcPurs_sub[0][i]->SetMarkerStyle(20);
		dataPurs_sub[0][i]->GetYaxis()->SetRangeUser(0,1);
		mcPurs[0][i]->SetLineColor(colorArr[i]);
		mcPurs_sub[0][i]->SetMarkerColor(colorArr[i]);
		mcPurs_sub[0][i]->SetLineColor(colorArr[i]);
		dataPurs_sub[0][i]->Draw();
		mcPurs[0][i]->SetMarkerColor(colorArr[i]);
		//mcPurs[0][i]->Draw("same");
		mcPurs_sub[0][i]->Draw("same");
		//dataPurs[0][i]->Draw("same");
		llat[i]->Draw("same");
		if(i==1) ll->Draw();
	}
	
	TCanvas *cc5 = new TCanvas("cc5","",1200,600);
	cc5->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cc5->cd(i+1);
		//sfs[0][i]->SetXTitle("Subjet pT");
		sfs[0][i]->SetYTitle("Scale Factors");
		sfs[0][i]->SetLineColor(1);
		sfs[0][i]->SetMarkerStyle(20);
		sfs[0][i]->Draw();
		//sfs[i]->GetYaxis()->SetRangeUser(0,1);
		sfsLight[0][i]->SetMarkerColor(6);
		sfsLight[0][i]->SetLineColor(6);
		sfsLight[0][i]->Draw("same");
		llat[i]->Draw("same");
	}
	
	
}