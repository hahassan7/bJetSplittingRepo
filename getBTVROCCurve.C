
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include <iostream>

using std::cout;
using std::endl;

void getBTVROCCurve(){
	
	bool doSubjets = true;
	
	TFile *ff = new TFile("gspTreeOut_ppMC_Pythia6_OfficialQCDJetMC_explicitJTA_recalibJP_withNegCSV.root");
	
	TTree *mcT = (TTree*)ff->Get("gspTree");
	double jtpt, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR;
	int jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
	double weight;
	int jtHadronFlavor;
	double discr_csvV2;
	mcT->SetBranchAddress("jtpt",&jtpt);
	mcT->SetBranchAddress("jtHadronFlavor",&jtHadronFlavor);
	mcT->SetBranchAddress("discr_csvV2",&discr_csvV2);
	mcT->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	mcT->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	mcT->SetBranchAddress("minCSV",&minCSV);
	mcT->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	mcT->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	mcT->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
	mcT->SetBranchAddress("weight",&weight);
	
	double nBjetsTag[30];
	double nJetsTag[30];
	double nLightJetsTag[30];
	double nLightTotal=0;
	double nBTotal=0;
	
	for(int i=0; i<30; i++){
		nBjetsTag[i]=0;
		nJetsTag[i]=0;
		nLightJetsTag[i]=0;
	}
	
	for(int ientry=0; ientry<mcT->GetEntries(); ientry++){
		mcT->GetEntry(ientry);
		
		if(jtpt<80) continue;
		
		int icsv=0;
		if(doSubjets) while(minCSV>(icsv)/30.) icsv++;
		else while(discr_csvV2>(icsv)/30.) icsv++;
		
		if(doSubjets){
			if(abs(jtSubJetHadronFlavor1)!=5 && abs(jtSubJetHadronFlavor2)!=5 && abs(jtSubJetHadronFlavor1)!=4 && abs(jtSubJetHadronFlavor2)!=4) nLightTotal+=weight;
			if(abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5) nBTotal +=weight;
		}
		else{
			if(abs(jtHadronFlavor)!=4 || abs(jtHadronFlavor)!=5) nLightTotal+=weight;
			if(abs(jtHadronFlavor)==5) nBTotal +=weight;
		}
		
		for(int cc=0; cc<icsv; cc++){
			nJetsTag[cc]+=weight;
			if(doSubjets){
				if(abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5) nBjetsTag[cc]+=weight;				
				if(abs(jtSubJetHadronFlavor1)!=5 && abs(jtSubJetHadronFlavor2)!=5 && abs(jtSubJetHadronFlavor1)!=4 && abs(jtSubJetHadronFlavor2)!=4) nLightJetsTag[cc]+=weight;
			}
			else{
				if(abs(jtHadronFlavor)==5) nBjetsTag[cc]+=weight;				
				else if(abs(jtHadronFlavor)!=4) nLightJetsTag[cc]+=weight;
			}
		}	
	}
	
	double bJetPurity[30];
	double bJetPurityErr[30];
	double lightRej[30];
	double lightRejErr[30];
	for(int i=0; i<30; i++){
		bJetPurity[i] = nBjetsTag[i]/nJetsTag[i];
		bJetPurityErr[i] = bJetPurity[i]*sqrt(nBjetsTag[i]+nJetsTag[i]);
		lightRej[i] = nLightJetsTag[i]/nLightTotal;
		lightRejErr[i] = lightRej[i]*sqrt(nLightJetsTag[i]+nLightTotal);
		
		cout << "bin " << i << " nJetsTag: "<< nJetsTag[i] << " nBjetsTag: "<< nBjetsTag[i] << " nLightTag: "<< nLightJetsTag[i] << endl;
	}
	
	TGraphErrors *hROC = new TGraphErrors(30,bJetPurity,lightRej,bJetPurityErr,lightRejErr);
	hROC->Draw("alp");
	hROC->GetXaxis()->SetTitle("di-b-jet Purity");
	hROC->GetYaxis()->SetTitle("di-usdg Purity");
	
}