#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

void getMCEff(){
	
	const int ptBins = 2;
	TH1D *tagged[ptBins];
	TH1D *total[ptBins];
	const double xbins[9] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4};
	for(int i=0; i<ptBins; i++){
		tagged[i] = new TH1D(Form("tagged_%d",i),"",8,xbins);
		total[i] = new TH1D(Form("total_%d",i),"",8,xbins);
	}
		
	TFile *fin_mc = new TFile("gspTreeOut_ppMC_muFilteredMC_Pythia8_fullSample.root");
	TTree *tmc = (TTree*)fin_mc->Get("gspTree");
	
	double jtpt, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, minJP, discr_csvV2, ndiscr_csvV2, discr_prob, svtxm;
	double jtSubJetJP1, jtSubJetJP2, jtSubJetCSV1, jtSubJetCSV2, jtSubJetNegCSV1, jtSubJetNegCSV2, jtSubJetSvtxm1, jtSubJetSvtxm2;
	int hiBin;
	double weight;
	int jtHadronFlavor, jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
	tmc->SetBranchAddress("jtpt",&jtpt);
	tmc->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tmc->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tmc->SetBranchAddress("minCSV",&minCSV);
	tmc->SetBranchAddress("minJP",&minJP);
	tmc->SetBranchAddress("svtxm",&svtxm);
	tmc->SetBranchAddress("discr_prob", &discr_prob);
	tmc->SetBranchAddress("ndiscr_csvV2",&ndiscr_csvV2);
	tmc->SetBranchAddress("discr_csvV2",&discr_csvV2);
	tmc->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	tmc->SetBranchAddress("jtHadronFlavor",&jtHadronFlavor);
	tmc->SetBranchAddress("weight",&weight);
	tmc->SetBranchAddress("jtSubJetJP_1",&jtSubJetJP1);
	tmc->SetBranchAddress("jtSubJetJP_2",&jtSubJetJP2);
	tmc->SetBranchAddress("jtSubJetcsvV2_1",&jtSubJetCSV1);
	tmc->SetBranchAddress("jtSubJetcsvV2_2",&jtSubJetCSV2);
	tmc->SetBranchAddress("jtSubJetNegCsvV2_1",&jtSubJetNegCSV1);
	tmc->SetBranchAddress("jtSubJetNegCsvV2_2",&jtSubJetNegCSV2);
	tmc->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	tmc->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
	tmc->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	tmc->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	tmc->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	for(int ievt=0; ievt<tmc->GetEntries(); ievt++){
		tmc->GetEntry(ievt);
		
		if(weight>1e9) continue;
		if(jtSubJetdR<0.1) continue;
		
		int ptBin=0;
		if(jtpt<100) continue;
		if(jtpt>140) ptBin=1;
		
		if(jtSubJetCSV1 > 0.8 && jtSubJetCSV2 > 0.8 && jtSubJetSvtxm1>1.5 && jtSubJetSvtxm2 > 1.5 && abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5) tagged[ptBin]->Fill(jtSubJetdR, weight);
		if(abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5) total[ptBin]->Fill(jtSubJetdR, weight);
	}
	
	for(int i=0; i<ptBins; i++){
		tagged[i]->Divide(total[i]);
	}
	
	TFile *fout = new TFile("drMCTagEff.root","recreate");
	fout->cd();
	for(int i=0; i<ptBins; i++){
		tagged[i]->Write();
	}
	
}