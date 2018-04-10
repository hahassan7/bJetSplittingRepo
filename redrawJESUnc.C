

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include <iostream>

using std::cout;
using std::endl;

double getzg(double j1, double j2){
	return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}


void redrawJESUnc(){
	
	TFile *f1 = new TFile("ntuples/gspTreeOut_ppMC_QCDJet_noMuFilter_withRefInfo.root");
	TFile *f2 = new TFile("ntuples/gspTreeOut_ppMC_QCDJet_noMuFilter_drop3p9pctTracks.root");
	
	const int ptBins = 2;
	const int jetPtBin[ptBins+1] = {120, 160, 400};
	const int jetPtBinShift[ptBins+1] = {116, 155, 400};
	
	const bool doBJets = true;
	
	TH1D *zgSpecIncl_mc[ptBins], *zgSpecIncl_mc_JES[ptBins], *zgSpecIncl_refAxis_mc[ptBins], *hRatio[ptBins];
	TH1D *drSpecIncl_mc[ptBins], *drSpecIncl_mc_JES[ptBins], *drSpecIncl_refAxis_mc[ptBins], *hdRRatio[ptBins];
	
	for(int i=0; i<ptBins; i++){
		zgSpecIncl_mc[i] = new TH1D(Form("zgSpecIncl_mc%d",i),"",10,0,0.5); zgSpecIncl_mc[i]->Sumw2();
		zgSpecIncl_mc_JES[i] = new TH1D(Form("zgSpecIncl_mc_JES%d",i),"",10,0,0.5); zgSpecIncl_mc_JES[i]->Sumw2();
		zgSpecIncl_refAxis_mc[i] = new TH1D(Form("zgSpecIncl_refAxis_mc%d",i),"",10,0,0.5); zgSpecIncl_refAxis_mc[i]->Sumw2();
		hRatio[i] = new TH1D(Form("hRatio%d",i),"",10,0,0.5); hRatio[i]->Sumw2();
		
		drSpecIncl_mc[i] = new TH1D(Form("drSpecIncl_mc%d",i),"",20,0,0.4); drSpecIncl_mc[i]->Sumw2();
		drSpecIncl_mc_JES[i] = new TH1D(Form("drSpecIncl_mc_JES%d",i),"",20,0,0.4); drSpecIncl_mc_JES[i]->Sumw2();
		drSpecIncl_refAxis_mc[i] = new TH1D(Form("drSpecIncl_refAxis_mc%d",i),"",20,0,0.4); drSpecIncl_refAxis_mc[i]->Sumw2();
	}
	
	double jtpt, jteta, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, jtSubJetRefdR, minJP, jtSubJetSvtxm1, jtSubJetSvtxm2;
	int jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
	double jtSubJetHadronDR1, jtSubJetHadronDR2, jtSubJetHadronPt1, jtSubJetHadronPt2;
	double weight;
	
	TF1 *randSmearGaus = new TF1("randSmearGaus","gaus",-1,1);
	randSmearGaus->SetParameter(0,1);
	randSmearGaus->SetParameter(1,0);
	
	TTree *t1 = (TTree*)f1->Get("gspTree");
	t1->SetBranchAddress("matchedJtpt",&jtpt);
	t1->SetBranchAddress("jteta",&jteta);
	t1->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	t1->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	t1->SetBranchAddress("minCSV",&minCSV);
	t1->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	t1->SetBranchAddress("jtSubJetRefdR",&jtSubJetRefdR);
	t1->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	t1->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
	t1->SetBranchAddress("jtSubJetHadronDR1",&jtSubJetHadronDR1);
	t1->SetBranchAddress("jtSubJetHadronDR2",&jtSubJetHadronDR2);
	t1->SetBranchAddress("jtSubJetHadronPt1",&jtSubJetHadronPt1);
	t1->SetBranchAddress("jtSubJetHadronPt2",&jtSubJetHadronPt2);
	t1->SetBranchAddress("minJP",&minJP);
	t1->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	t1->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	t1->SetBranchAddress("weight",&weight);
	for(int ientry=0; ientry<t1->GetEntries(); ientry++){
		t1->GetEntry(ientry);
		
		double resolution = 0.06 + 0.9525/jtpt - 1.709e-05/jtpt/jtpt; //from AN2015-257-v18 p26
		randSmearGaus->SetParameter(2,resolution);
		double randSmear = randSmearGaus->GetRandom(-0.5,0.5);
		double jtptSmeared=(1.+randSmear)*jtpt;
		
		double zg = getzg(jtSubJetPt1, jtSubJetPt2);
		if(zg<0.1) continue;
		if(doBJets && (minCSV<0.8 || jtSubJetSvtxm1<1.2 || jtSubJetSvtxm2<1.2)) continue;
		
		//cout << "jtpt: "<< jtpt << " jtptSmeared: "<< jtptSmeared << endl;
		if(jtpt>120 && jtpt<400){
			int ibin=0;
			while(jtpt>jetPtBin[ibin+1]) ibin++;
			
			if(jtSubJetdR>0.1) zgSpecIncl_mc[ibin]->Fill(zg, weight);
			if(jtSubJetdR>0.1) drSpecIncl_mc[ibin]->Fill(jtSubJetdR, weight);
			if(jtSubJetRefdR>0.1) drSpecIncl_refAxis_mc[ibin]->Fill(jtSubJetRefdR, weight);
		}
		
		if(jtptSmeared>120 && jtptSmeared<400){
			int ibin=0;
			while(jtptSmeared>jetPtBin[ibin+1]) ibin++;
			if(jtSubJetdR>0.1) zgSpecIncl_mc_JES[ibin]->Fill(zg, weight);
			if(jtSubJetdR>0.1) drSpecIncl_mc_JES[ibin]->Fill(jtSubJetdR, weight);
		}
	}
		
	for(int i=0; i<ptBins; i++){
		zgSpecIncl_mc[i]->Scale(1./zgSpecIncl_mc[i]->Integral());
		zgSpecIncl_mc[i]->Scale(1./zgSpecIncl_mc[i]->GetBinWidth(4));
		
		drSpecIncl_mc[i]->Scale(1./drSpecIncl_mc[i]->Integral());
		drSpecIncl_mc[i]->Scale(1./drSpecIncl_mc[i]->GetBinWidth(4));
		
		drSpecIncl_refAxis_mc[i]->Scale(1./drSpecIncl_refAxis_mc[i]->Integral());
		drSpecIncl_refAxis_mc[i]->Scale(1./drSpecIncl_refAxis_mc[i]->GetBinWidth(4));
	}
	
	//////////////////////////////////////////
	
	TTree *t2 = (TTree*)f2->Get("gspTree");
	t2->SetBranchAddress("matchedJtpt",&jtpt);
	t2->SetBranchAddress("jteta",&jteta);
	t2->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	t2->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	t2->SetBranchAddress("minCSV",&minCSV);
	t2->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	t2->SetBranchAddress("jtSubJetRefdR",&jtSubJetRefdR);
	t2->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	t2->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
	t2->SetBranchAddress("jtSubJetHadronDR1",&jtSubJetHadronDR1);
	t2->SetBranchAddress("jtSubJetHadronDR2",&jtSubJetHadronDR2);
	t2->SetBranchAddress("jtSubJetHadronPt1",&jtSubJetHadronPt1);
	t2->SetBranchAddress("jtSubJetHadronPt2",&jtSubJetHadronPt2);
	t2->SetBranchAddress("minJP",&minJP);
	t2->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	t2->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	t2->SetBranchAddress("weight",&weight);
	for(int ientry=0; ientry<t2->GetEntries(); ientry++){
		t2->GetEntry(ientry);
		
		if(jtpt<120) continue;
		if(jtpt>400) continue;
		if(doBJets && (minCSV<0.8 || jtSubJetSvtxm1<1.2 || jtSubJetSvtxm2<1.2)) continue;
		
		int ibin=0;
		while(jtpt>jetPtBin[ibin+1]) ibin++;
		
		double zg = getzg(jtSubJetPt1, jtSubJetPt2);
		if(zg<0.1) continue;
		//if(jtSubJetdR>0.1 && minCSV>0.8 && jtSubJetSvtxm1>1.2 && jtSubJetSvtxm2>1.2) zgSpecIncl_mc_JES[ibin]->Fill(zg, weight);
		//if(jtSubJetdR>0.1 && minCSV>0.8 && jtSubJetSvtxm1>1.2 && jtSubJetSvtxm2>1.2) drSpecIncl_mc_JES[ibin]->Fill(jtSubJetdR, weight);
		
	}
	
	for(int i=0; i<ptBins; i++){
		zgSpecIncl_mc_JES[i]->Scale(1./zgSpecIncl_mc_JES[i]->Integral());
		zgSpecIncl_mc_JES[i]->Scale(1./zgSpecIncl_mc_JES[i]->GetBinWidth(4));
		
		drSpecIncl_mc_JES[i]->Scale(1./drSpecIncl_mc_JES[i]->Integral());
		drSpecIncl_mc_JES[i]->Scale(1./drSpecIncl_mc_JES[i]->GetBinWidth(4));
	}
	
	TCanvas *cJES = new TCanvas("cJES","",1200,600);
	TPad *p1[2], *p2[2];
	TLatex *labelsZg[ptBins];
	TLegend *tl1 = new TLegend(0.455,0.555,0.855,0.735);
	double meanBin[2] = {0,0};
	cJES->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cJES->cd(i+1);
		p1[i] = new TPad(Form("p1_%d",i),"",0,0.3,1.0,1.0);
		p1[i]->SetBottomMargin(0.01);
		p2[i] = new TPad(Form("p2_%d",i),"",0,0,1.0,0.3);
		p2[i]->SetTopMargin(0.01);
		p2[i]->SetBottomMargin(0.4);
		p1[i]->Draw();
		p2[i]->Draw();
		p1[i]->cd();
		zgSpecIncl_mc[i]->SetYTitle("1/N dN/dz_{g}");
		zgSpecIncl_mc[i]->Draw();
		zgSpecIncl_mc_JES[i]->SetLineColor(2);
		zgSpecIncl_mc_JES[i]->SetMarkerColor(2);
		zgSpecIncl_mc_JES[i]->Draw("same");
		if(!i){
			tl1->AddEntry(zgSpecIncl_mc[i],"Nominal Distribution");
			tl1->AddEntry(zgSpecIncl_mc_JES[i],"Dropped 3.9\% of Tracks");
			tl1->Draw("same");
		}
		labelsZg[i] = new TLatex(0.2,7.0,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		labelsZg[i]->Draw("same");
		p2[i]->cd();
		hRatio[i] = (TH1D*)zgSpecIncl_mc[i]->Clone(Form("hRatio_%d",i));
		hRatio[i]->Divide(zgSpecIncl_mc[i], zgSpecIncl_mc_JES[i], 1, 1, "B");
		hRatio[i]->SetMaximum(1.15);
		hRatio[i]->SetMinimum(0.85);
		hRatio[i]->SetLabelSize(0.12,"X");
		hRatio[i]->SetLabelSize(0.12,"Y");
		hRatio[i]->SetNdivisions(505,"Y");
		hRatio[i]->SetTitleSize(0.13,"X");
		hRatio[i]->SetTitleSize(0.13,"Y");
		hRatio[i]->SetTitleOffset(0.35,"Y");
		hRatio[i]->SetTitleOffset(0.85,"X");
		hRatio[i]->SetTickLength(0.1,"X");
		hRatio[i]->SetYTitle("Ratio");
		hRatio[i]->SetXTitle("z_{g}");
		for(int ibin=3; ibin<=10; ibin++){
			meanBin[i]+=hRatio[i]->GetBinContent(ibin);
		}
		cout << " mean, bin " << i << " " << meanBin[i]/8. << endl;
		hRatio[i]->Draw();
	}
	
	TCanvas *cRefAxis = new TCanvas("cRefAxis","",1200,600);
	TPad *p3[2], *p4[2];
	TLatex *labelsdr[ptBins];
	TLegend *tl2 = new TLegend(0.21,0.056,0.615,0.279);
	cRefAxis->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cRefAxis->cd(i+1);
		p3[i] = new TPad(Form("p3_%d",i),"",0,0.3,1.0,1.0);
		p3[i]->SetBottomMargin(0.01);
		p4[i] = new TPad(Form("p4_%d",i),"",0,0,1.0,0.3);
		p4[i]->SetTopMargin(0.01);
		p4[i]->SetBottomMargin(0.4);
		p3[i]->Draw();
		p4[i]->Draw();
		p3[i]->cd();
		drSpecIncl_mc[i]->SetYTitle("1/N dN/d#Delta R");
		drSpecIncl_mc[i]->Draw();
		drSpecIncl_mc_JES[i]->SetLineColor(2);
		drSpecIncl_mc_JES[i]->SetMarkerColor(2);
		drSpecIncl_mc_JES[i]->Draw("same");
		if(!i){
			tl2->AddEntry(drSpecIncl_mc[i],"Nominal Distribution");
			tl2->AddEntry(drSpecIncl_mc_JES[i],"Dropped 3.9\% of Tracks");
			tl2->Draw("same");
		}
		labelsdr[i] = new TLatex(0.2,6.0,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		labelsdr[i]->Draw("same");
		p4[i]->cd();
		hdRRatio[i] = (TH1D*)drSpecIncl_mc[i]->Clone(Form("hRatio_%d",i));
		hdRRatio[i]->Divide(drSpecIncl_mc[i], drSpecIncl_mc_JES[i], 1, 1, "B");
		hdRRatio[i]->SetMaximum(1.15);
		hdRRatio[i]->SetMinimum(0.85);
		hdRRatio[i]->SetLabelSize(0.12,"X");
		hdRRatio[i]->SetLabelSize(0.12,"Y");
		hdRRatio[i]->SetNdivisions(505,"Y");
		hdRRatio[i]->SetTitleSize(0.13,"X");
		hdRRatio[i]->SetTitleSize(0.13,"Y");
		hdRRatio[i]->SetTitleOffset(0.35,"Y");
		hdRRatio[i]->SetTitleOffset(0.85,"X");
		hdRRatio[i]->SetTickLength(0.1,"X");
		hdRRatio[i]->SetYTitle("Ratio");
		hdRRatio[i]->SetXTitle("#Delta R");
		hdRRatio[i]->Draw();
	}
	
}