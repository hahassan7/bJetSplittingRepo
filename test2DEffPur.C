

#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

double getzg(double j1, double j2){
	return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}

void test2DEffPur(){
	
	TH2D *eff = new TH2D("eff","",10,0,0.5,20,0,3);
	TH2D *effDen = new TH2D("effDen","",10,0,0.5,20,0,3);
	TH2D *pur = new TH2D("pur","",10,0,0.5,20,0,3);
	TH2D *purDen = new TH2D("purDen","",10,0,0.5,20,0,3);
	
	TFile *fin = new TFile("gspTreeOut_ppMC_Pythia6_officialQCDMC_explicitJTA_recalibJP_pthat80plus_trgFiltered.root");
	TTree *tin = (TTree*)fin->Get("gspTree");
	
	double jtpt, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, minJP, jtSubJetSvtxm1, jtSubJetSvtxm2, jtSubJetJP_1, jtSubJetJP_2;
	int jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
	double weight;
	tin->SetBranchAddress("jtpt",&jtpt);
	tin->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tin->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tin->SetBranchAddress("minCSV",&minCSV);
	tin->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	tin->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	tin->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
	tin->SetBranchAddress("minJP",&minJP);
	tin->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	tin->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	tin->SetBranchAddress("jtSubJetJP_1",&jtSubJetJP_1);
	tin->SetBranchAddress("jtSubJetJP_2",&jtSubJetJP_2);
	tin->SetBranchAddress("weight",&weight);
	
	for(int ientry=0; ientry<tin->GetEntries(); ientry++){
		tin->GetEntry(ientry);
		
		if(jtpt<100) continue;
			
		//int ibin=0;
		//while(jtpt>jetPtBin[ibin+1]) ibin++;
		
		double zg = getzg(jtSubJetPt1, jtSubJetPt2);
		
		bool tag = (minCSV>0.8) ? true : false; 
		
		if(jtSubJetHadronFlavor1==5 && jtSubJetHadronFlavor2==5){
			effDen->Fill(zg,jtSubJetJP_1,weight);
			effDen->Fill(zg,jtSubJetJP_2,weight);
			if(tag){
				eff->Fill(zg,jtSubJetJP_1,weight);
				eff->Fill(zg,jtSubJetJP_2,weight);
				pur->Fill(zg,jtSubJetJP_1,weight);
				pur->Fill(zg,jtSubJetJP_2,weight);
			}
		}
		if(tag){
			purDen->Fill(zg,jtSubJetJP_1,weight);
			purDen->Fill(zg,jtSubJetJP_2,weight);
		}
	}
	eff->Divide(effDen);
	pur->Divide(purDen);
	
	TCanvas *c2 = new TCanvas("c2","",1200,600);
	c2->Divide(2,1);
	c2->cd(1);
	eff->SetXTitle("zG");
	eff->SetYTitle("Single JP Discriminator");
	eff->SetTitle("Efficiency");
	eff->Draw("colz");
	c2->cd(2);
	pur->SetXTitle("zG");
	pur->SetYTitle("Single JP Discriminator");
	pur->SetTitle("Purity");
	pur->Draw("colz");
}