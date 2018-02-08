
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"

#include <iostream>

using std::cout;
using std::endl;

void normalizeZg(TH1D *input){
	input->Scale(1./input->Integral());
	input->Scale(1./input->GetBinWidth(1));
}

double getZg(double j1, double j2){
	return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}

void correctSubjetReco(TF1 *corrFactor, TH1D *input, bool doDivide=0){
	for(int ibin=1; ibin<=input->GetNbinsX(); ibin++){
		double binCenter = input->GetBinCenter(ibin);
		if(!doDivide){
			input->SetBinContent(ibin, input->GetBinContent(ibin)*corrFactor->Eval(binCenter));
			input->SetBinError(ibin, input->GetBinError(ibin)*corrFactor->Eval(binCenter));
		}
		else{
			input->SetBinContent(ibin, input->GetBinContent(ibin)/corrFactor->Eval(binCenter));
			input->SetBinError(ibin, input->GetBinError(ibin)/corrFactor->Eval(binCenter));
		}
	}
	
}

void subjetAngRes(bool ispp = true){
	
	const int ptBins = 2;
	const double xbins[ptBins+1] = {100, 140, 9999};
	
	TH1D *zgNominal[ptBins];
	TH1D *zgAngUp[ptBins];
	TH1D *zgAngDown[ptBins];
	
	TH1D *zgBNominal[ptBins];
	TH1D *zgBAngUp[ptBins];
	TH1D *zgBAngDown[ptBins];
	
	for(int i=0; i<ptBins; i++){
		zgNominal[i] = new TH1D(Form("zgNominal_pt%d",i),"",10,0,0.5); zgNominal[i]->Sumw2();
		zgAngUp[i] = new TH1D(Form("zgAngUp_pt%d",i),"",10,0,0.5); zgAngUp[i]->Sumw2();
		zgAngDown[i] = new TH1D(Form("zgAngDown_pt%d",i),"",10,0,0.5); zgAngDown[i]->Sumw2();
		
		zgBNominal[i] = new TH1D(Form("zgBNominal_pt%d",i),"",10,0,0.5); zgBNominal[i]->Sumw2();
		zgBAngUp[i] = new TH1D(Form("zgBAngUp_pt%d",i),"",10,0,0.5); zgBAngUp[i]->Sumw2();
		zgBAngDown[i] = new TH1D(Form("zgBAngDown_pt%d",i),"",10,0,0.5); zgBAngDown[i]->Sumw2();
	}
	
	TFile *fin = new TFile("gspTreeOut_ppData_explicitJTA_recalibJP_pt60.root");
	
	TTree *tin = (TTree*)fin->Get("gspTree");
	double jtpt, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, minJP, jtSubJetSvtxm1, jtSubJetSvtxm2;
	int hiBin;
	tin->SetBranchAddress("jtpt",&jtpt);
	tin->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tin->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tin->SetBranchAddress("minCSV",&minCSV);
	tin->SetBranchAddress("minJP",&minJP);
	tin->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	tin->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	tin->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	if(!ispp) tin->SetBranchAddress("hiBin",&hiBin);
	
	for(int ientry=0; ientry<tin->GetEntries(); ientry++){
		tin->GetEntry(ientry);
				
		if(jtpt<100) continue;
		int ipt=0;
		while(jtpt>xbins[ipt+1]) ipt++;
				
		if(jtSubJetdR>0.1) zgNominal[ipt]->Fill(getZg(jtSubJetPt1,jtSubJetPt2));
		if(jtSubJetdR>0.11) zgAngUp[ipt]->Fill(getZg(jtSubJetPt1,jtSubJetPt2));
		if(jtSubJetdR>0.09) zgAngDown[ipt]->Fill(getZg(jtSubJetPt1,jtSubJetPt2));
		
		if(minCSV>0.8 && jtSubJetSvtxm1>1.5 && jtSubJetSvtxm2>1.5){
			if(jtSubJetdR>0.1) zgBNominal[ipt]->Fill(getZg(jtSubJetPt1,jtSubJetPt2));
			if(jtSubJetdR>0.11) zgBAngUp[ipt]->Fill(getZg(jtSubJetPt1,jtSubJetPt2));
			if(jtSubJetdR>0.09) zgBAngDown[ipt]->Fill(getZg(jtSubJetPt1,jtSubJetPt2));
		}
	}
	
	TFile *fcorr = new TFile("zgCorrFactors_gsp.root");
	TF1 *recoCorr[ptBins];
	TF1 *tagCorr[ptBins];
	for(int i=0; i<ptBins; i++){
		
		recoCorr[i] = (TF1*)fcorr->Get(Form("fCorrFactor_pt%d",i));
		tagCorr[i] = (TF1*)fcorr->Get(Form("fTagCorrFactor_pt%d",i));
		correctSubjetReco(recoCorr[i], zgBNominal[i]);
		correctSubjetReco(recoCorr[i], zgBAngUp[i]);
		correctSubjetReco(recoCorr[i], zgBAngDown[i]);
		
		correctSubjetReco(tagCorr[i], zgBNominal[i], 1);
		correctSubjetReco(tagCorr[i], zgBAngUp[i], 1);
		correctSubjetReco(tagCorr[i], zgBAngDown[i], 1);
		
		normalizeZg(zgNominal[i]);
		normalizeZg(zgAngUp[i]);
		normalizeZg(zgAngDown[i]);
		normalizeZg(zgBNominal[i]);
		normalizeZg(zgBAngUp[i]);
		normalizeZg(zgBAngDown[i]);
	}
	
	TH1D *hRatiosUp[ptBins];
	TH1D *hRatiosDown[ptBins];
	TLatex *labels[ptBins];
	TLegend *tleg = new TLegend(0.5,0.5,0.9,0.9);
	TCanvas *cc = new TCanvas("cc","",1200,800);
	cc->Divide(ptBins,3);
	for(int i=0; i<ptBins; i++){
		cc->cd(i+1);
		zgNominal[i]->SetYTitle("1/N dN/dz_{g}");
		zgNominal[i]->SetXTitle("Inclusive z_{g}");
		zgNominal[i]->Draw();
		zgAngUp[i]->SetMarkerColor(4);
		zgAngUp[i]->SetMarkerStyle(33);
		zgAngDown[i]->SetMarkerColor(2);
		zgAngDown[i]->SetMarkerStyle(25);
		zgAngDown[i]->Draw("same");
		zgAngUp[i]->Draw("same");
		labels[i] = new TLatex(0.25,6.5,Form("%g < jet p_{T} < %g",xbins[i],xbins[i+1]));
		labels[i]->Draw("same");
		if(i==0){
			tleg->AddEntry(zgNominal[i],"#DeltaR > 0.1");
			tleg->AddEntry(zgAngDown[i],"#DeltaR > 0.09");
			tleg->AddEntry(zgAngUp[i],"#DeltaR > 0.11");
			tleg->Draw("same");
		}
	}
	for(int i=0; i<ptBins; i++){
		cc->cd(ptBins+i+1);
		zgBNominal[i]->SetYTitle("1/N dN/dz_{g}");
		zgBNominal[i]->SetXTitle("GSP z_{g}");
		zgBNominal[i]->Draw();
		zgBAngUp[i]->SetMarkerColor(4);
		zgBAngUp[i]->SetMarkerStyle(33);
		zgBAngDown[i]->SetMarkerColor(2);
		zgBAngDown[i]->SetMarkerStyle(25);
		zgBAngDown[i]->Draw("same");
		zgBAngUp[i]->Draw("same");
	}
	for(int i=0; i<ptBins; i++){
		cc->cd(ptBins*2+i+1);
		hRatiosUp[i] = (TH1D*)zgBNominal[i]->Clone(Form("hRatiosUp_pt%d",i));
		hRatiosUp[i]->Divide(hRatiosUp[i],zgBAngUp[i],1,1,"B");
		hRatiosUp[i]->SetMarkerColor(4);
		hRatiosUp[i]->SetMarkerStyle(33);
		hRatiosUp[i]->GetYaxis()->SetRangeUser(0.85,1.15);
		hRatiosUp[i]->SetYTitle("Ratios");
		hRatiosUp[i]->SetXTitle("GSP z_{g}");
		hRatiosUp[i]->Draw();
		hRatiosDown[i] = (TH1D*)zgBNominal[i]->Clone(Form("hRatiosDown_pt%d",i));
		hRatiosDown[i]->Divide(hRatiosDown[i],zgBAngDown[i],1,1,"B");
		hRatiosDown[i]->SetMarkerColor(2);
		hRatiosDown[i]->SetMarkerStyle(25);
		hRatiosDown[i]->Draw("same");
	}
	
	
	
	
}