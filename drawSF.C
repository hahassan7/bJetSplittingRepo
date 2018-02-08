
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <iostream>

using std::cout;
using std::endl;

void drawSF(){
	
	const int nBins = 3;
	const string workpts[nBins] = {"0.5","0.7","0.8"};
	const string legendTitle[nBins] = {"CSV>0.5","CSV>0.7","CSV>0.8"};
	const string legendTitleFull[nBins] = {"Fulljet CSV>0.5","Fulljet CSV>0.7","Fulljet CSV>0.8"};
	const int colorArr[nBins] = {1,2,4};
	
	TFile *ff[nBins];
	TFile *ff_light[nBins];
	TH1D *sfs[nBins];
	TH1D *dataEffs[nBins];
	TH1D *mcEffs[nBins];
	
	TH1D *dataPurs[nBins];
	TH1D *mcPurs[nBins];
	TH1D *sfsLight[nBins];
	
	TH1D *dataPurs_sub[nBins];
	TH1D *mcPurs_sub[nBins];
	
	TFile *ff_full[nBins];
	TH1D *sfs_full[nBins];
	TH1D *dataEff_full[nBins];
	TH1D *mcEff_full[nBins];
	for(int i=0; i<nBins; i++){
		if(i==1) ff[i] = new TFile("sf_subjets_zg__muFiltered_csv0.8_sv1p5_dataJPCalibOnMC.root");
		else ff[i] = new TFile(Form("sf_subjets_zg__muFiltered_csv%s_sv1p5.root",workpts[i].c_str()));
		//ff_light[i] = new TFile(Form("sf_lightJets_subjets_csv%s.root",workpts[i].c_str()));
		ff_full[i] = new TFile(Form("sf_fullJets_csv%s.root",workpts[i].c_str()));
		
		sfs[i] = (TH1D*)ff[i]->Get("scaleFactors_btagSF");
		dataEffs[i] = (TH1D*)ff[i]->Get("dataEff");
		dataEffs[i]->SetMarkerStyle(24); dataEffs[i]->SetLineColor(1); dataEffs[i]->SetMarkerColor(1);
		mcEffs[i] = (TH1D*)ff[i]->Get("mcEff");
		sfs_full[i] = (TH1D*)ff_full[i]->Get("scaleFactors_btagSF");
		dataEff_full[i] = (TH1D*)ff_full[i]->Get("dataEff");
		dataEff_full[i]->SetMarkerStyle(24); dataEff_full[i]->SetLineColor(1); dataEff_full[i]->SetMarkerColor(1);
		mcEff_full[i] = (TH1D*)ff_full[i]->Get("mcEff");
				
		mcPurs[i] = (TH1D*)ff_full[i]->Get("mcPur");
		sfsLight[i] = (TH1D*)ff_full[i]->Get("sfLight");
		dataPurs[i] = (TH1D*)ff_full[i]->Get("dataPur");
				
		mcPurs_sub[i] = (TH1D*)ff[i]->Get("mcPur");
		dataPurs_sub[i] = (TH1D*)ff[i]->Get("dataPur");
	}

	
	TLegend *ll = new TLegend(0.5,0.5,0.9,0.9);
	TLegend *ll2 = new TLegend(0.1,0.1,0.5,0.5);
	TLatex *llat[nBins];
	TCanvas *cc = new TCanvas("cc","",1200,800);
	cc->Divide(nBins,2);
	for(int i=0; i<nBins; i++){
		cc->cd(i+1); //1,2,3
		mcEffs[i]->SetMarkerColor(colorArr[i]);
		mcEffs[i]->SetLineColor(colorArr[i]);
		mcEff_full[i]->SetMarkerColor(colorArr[i]);
		mcEff_full[i]->SetLineColor(colorArr[i]);
		dataEffs[i]->SetMarkerStyle(24);
		dataEff_full[i]->SetMarkerStyle(25);
		mcEff_full[i]->SetMarkerStyle(21);
		mcEffs[i]->SetMarkerStyle(20);
		dataEffs[i]->SetYTitle("Tagging Efficiency");
		dataEffs[i]->SetXTitle("Full(sub)-jet p_{T}");
		dataEffs[i]->Draw("e");
		dataEffs[i]->GetXaxis()->SetNdivisions(505);
		dataEffs[i]->GetYaxis()->SetRangeUser(0.,1.);
		mcEffs[i]->Draw("same");
		//mcEff_full[i]->Draw("same");
		//dataEff_full[i]->Draw("e,same");
		llat[i] = new TLatex(0.1,0.9,legendTitle[i].c_str());
		llat[i]->Draw("same");
		if(i==1){
			ll->AddEntry(mcEffs[i], "MC Subjets");
			ll->AddEntry(dataEffs[i], "Data Subjets");
			//ll->AddEntry(mcEff_full[i], "MC Full jets");
			//ll->AddEntry(dataEff_full[i], "Data Full jets");
			ll->Draw();
		} 
		
		cc->cd(i+1+nBins); //4,5,6
		sfs[i]->SetMarkerColor(colorArr[i]);
		sfs[i]->SetLineColor(colorArr[i]);
		sfs[i]->SetMarkerSize(1);
		sfs[i]->SetMarkerStyle(20);
		sfs[i]->SetXTitle("Full(sub)-jet p_{T}");
		sfs[i]->SetYTitle("Tag Eff. Scale Factor (Data/MC)");
		//sfs[i]->Draw();
		sfs[i]->Draw();
		sfs[i]->GetXaxis()->SetNdivisions(505);
		sfs[i]->GetYaxis()->SetRangeUser(0.5,1.5);
		//ll->AddEntry(sfs[i],legendTitle[i].c_str(),"lp");
		
		sfs_full[i]->SetMarkerStyle(25);
		sfs_full[i]->SetMarkerSize(1);
		sfs_full[i]->SetMarkerColor(colorArr[i]);
		sfs_full[i]->SetLineColor(colorArr[i]);
		sfs_full[i]->GetYaxis()->SetRangeUser(0.5,1.5);
		//sfs_full[i]->Draw("same");
		if(i==nBins-1){
			ll2->AddEntry(sfs[i],"Subjet SF","lp");
			ll2->AddEntry(sfs_full[i],"Full Jet SF","lp");
			//ll2->Draw();
		}
		
	}
	//ll->Draw("Same");
	//ll2->Draw("same");
	
	TCanvas *cc4 = new TCanvas("cc4","",1200,600);
	cc4->Divide(nBins,1);
	for(int i=0; i<nBins; i++){
		cc4->cd(i+1);
		dataPurs_sub[i]->SetXTitle("zG");
		dataPurs_sub[i]->SetYTitle("Tagging Purity");
		dataPurs[i]->SetLineColor(1);
		dataPurs_sub[i]->SetLineColor(1);
		dataPurs_sub[i]->SetMarkerStyle(24);
		dataPurs[i]->SetMarkerStyle(25);
		mcPurs[i]->SetMarkerStyle(21);
		mcPurs_sub[i]->SetMarkerStyle(20);
		dataPurs_sub[i]->GetYaxis()->SetRangeUser(0,1);
		mcPurs[i]->SetLineColor(colorArr[i]);
		mcPurs_sub[i]->SetMarkerColor(colorArr[i]);
		mcPurs_sub[i]->SetLineColor(colorArr[i]);
		dataPurs_sub[i]->Draw();
		mcPurs[i]->SetMarkerColor(colorArr[i]);
		mcPurs[i]->Draw("same");
		mcPurs_sub[i]->Draw("same");
		dataPurs[i]->Draw("same");
		llat[i]->Draw("same");
		if(i==1) ll->Draw();
	}
	
	TCanvas *cc5 = new TCanvas("cc5","",1200,600);
	cc5->Divide(nBins,1);
	for(int i=0; i<nBins; i++){
		cc5->cd(i+1);
		sfs[i]->SetXTitle("Subjet pT");
		sfs[i]->SetYTitle("Scale Factors");
		sfs[i]->SetLineColor(1);
		sfs[i]->SetMarkerStyle(20);
		sfs[i]->Draw();
		//sfs[i]->GetYaxis()->SetRangeUser(0,1);
		sfsLight[i]->SetMarkerColor(6);
		sfsLight[i]->SetLineColor(6);
		sfsLight[i]->Draw("same");
		llat[i]->Draw("same");
	}
	
	
}