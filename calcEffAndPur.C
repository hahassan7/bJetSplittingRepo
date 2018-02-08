
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include <iostream>

using std::endl;
using std::cout;

void calcEffAndPur(){
	
	TFile *fp6 = new TFile("gspTreeOut_PbPbMC_pthat120_akPu4CaloJets.root");
	TFile *f1 = new TFile("gspTreeOut_PbPbMC_explicitJTA.root");
	TTree *t = (TTree*)f1->Get("gspTree");
	TTree *tp6 = (TTree*)fp6->Get("gspTree");
	
	/*TH1D *taggedB = new TH1D("taggedB","taggedB",50,100,400); taggedB->Sumw2();
	TH1D *mcB = new TH1D("mcB","mcB",50,100,400); mcB->Sumw2();
	TH1D *tagged = new TH1D("tagged","tagged",50,100,400); tagged->Sumw2();
	
	t->Project("taggedB","jtpt","weight*(discr_csvV2>0.95 && abs(jteta)<2 && jtHadronFlavor==5)");
	t->Project("mcB","jtpt","weight*(abs(jteta)<2 && jtHadronFlavor==5)");
	t->Project("tagged","jtpt","weight*(discr_csvV2>0.95 && abs(jteta)<2)");*/
	
	TH1D *allB = new TH1D("allB","",440,-10,1);
	TH1D *all = new TH1D("all","",440,-10,1);
	
	TH1D *allB_p6 = new TH1D("allB_p6","",440,-10,1);
	TH1D *all_p6 = new TH1D("all_p6","",440,-10,1);
	
	t->Project("allB","discr_csvV2","(abs(jteta)<1.6 && jtpt>120 && jtHadronFlavor==5 && hiBin<60)");
	t->Project("all","discr_csvV2","(abs(jteta)<1.6 && jtpt>120 && hiBin<60)");
	
	tp6->Project("allB_p6","discr_csvV2","weight*(abs(jteta)<1.6 && jtpt>120 && abs(refparton_flavorForB)==5 && hiBin<60)");
	tp6->Project("all_p6","discr_csvV2","weight*(abs(jteta)<1.6 && jtpt>120 && hiBin<60)");
	
	cout << "pf bjets: "<< allB->GetEntries() << " calo bjets: "<<  allB_p6->GetEntries() << endl;
	
	double eff[40], pur[40];
	double effP6[40], purP6[40];
	double allBjets = allB->Integral();
	double allBjets_p6 = allB_p6->Integral();
	for(int i=401; i<=allB->GetNbinsX(); i++){
		double allJ=0, bCut=0;
		double allJ_p6=0, bCut_p6=0;
		for(int j=i; j<=all->GetNbinsX(); j++){
			allJ+=all->GetBinContent(j);
			bCut+=allB->GetBinContent(j);
			
			allJ_p6+=all_p6->GetBinContent(j);
			bCut_p6+=allB_p6->GetBinContent(j);
		}
		eff[i-401] = bCut/allBjets;
		pur[i-401] = bCut/allJ;
		
		effP6[i-401] = bCut_p6/allBjets_p6;
		purP6[i-401] = bCut_p6/allJ_p6;
		
	}
	
	for(int i=0; i<40; i++){
		cout << purP6[i] << ", " << effP6[i] << ", " << (double)i/40. << endl;
	}
	
	TGraph *roc = new TGraph(40, pur, eff);
	TGraph *rocP6 = new TGraph(40, purP6, effP6);
	roc->SetLineColor(kGreen+2);
	roc->SetMarkerColor(kGreen+2);
	roc->Draw("alp");
	roc->GetYaxis()->SetTitle("Tagging Efficiency");
	roc->GetXaxis()->SetTitle("Tagging Purity");
	
	rocP6->SetLineColor(kRed+1);
	rocP6->SetMarkerColor(kRed+1);
	rocP6->Draw("pl,same");
	
	double purWorkPt[1] = {pur[(int)(0.9/(1./40.))]};
	double effWorkPt[1] = {eff[(int)(0.9/(1./40.))]};
	
	double purWorkPtP6[1] = {purP6[(int)(0.9/(1./40.))]};
	double effWorkPtP6[1] = {effP6[(int)(0.9/(1./40.))]};
	
	//cout << "working point eff: "<< effWorkPt[0] << " pur: "<< purWorkPt[0] << endl;
	
	TGraph *workPt = new TGraph(1, purWorkPt, effWorkPt);
	workPt->SetMarkerColor(4);
	workPt->SetMarkerStyle(25);
	workPt->SetMarkerSize(1.4);
	workPt->Draw("p,same");
	
	TGraph *workPtP6 = new TGraph(1, purWorkPtP6, effWorkPtP6);
	workPtP6->SetMarkerColor(4);
	workPtP6->SetMarkerStyle(25);
	workPtP6->SetMarkerSize(1.4);
	workPtP6->Draw("p,same");
	
	
	//taggedB->Draw();
	//mcB->Draw("same");
	//tagged->Draw("same");
	
	/*TH1D *effHisto = (TH1D*)taggedB->Clone("eff");
	effHisto->Divide(mcB);
	
	TH1D *purHisto = (TH1D*)taggedB->Clone("pur");
	purHisto->Divide(tagged);
	
	//effHisto->Draw();
	purHisto->SetMarkerColor(2);
	purHisto->SetLineColor(2);*/
	//purHisto->Draw("same");
	
}