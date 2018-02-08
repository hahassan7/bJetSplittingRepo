#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"

void drawInvMass(){
	
	TFile *fppdata = new TFile("gspTreeOut_ppData_explicitJTA_recalibJP_pt60.root");
	TFile *fppmc = new TFile("gspTreeOut_ppMC_muFilteredMC_Pythia8_fullSample.root");
	
	const int ptBins = 2;
	const int jetPtBin[ptBins+1] = {100, 140, 9999};
	TH1D *invMassDataIncl[ptBins], *invMassDataGSP[ptBins];
	TH1D *invMassMCIncl[ptBins], *invMassMCGSP[ptBins];
	
	for(int i=0; i<ptBins; i++){
		invMassDataIncl[i] = new TH1D(Form("invMassDataIncl_%d",i),"",20,0,40); invMassDataIncl[i]->Sumw2();
		invMassDataGSP[i] = new TH1D(Form("invMassDataGSP_%d",i),"",20,0,40); invMassDataGSP[i]->Sumw2();
		invMassMCIncl[i] = new TH1D(Form("invMassMCIncl_%d",i),"",20,0,40); invMassMCIncl[i]->Sumw2();
		invMassMCGSP[i] = new TH1D(Form("invMassMCGSP_%d",i),"",20,0,40); invMassMCGSP[i]->Sumw2();
	}
	
	
	//Load data tree
	TTree *dataT = (TTree*)fppdata->Get("gspTree");
	double jtpt, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, minJP, jtSubJetSvtxm1, jtSubJetSvtxm2;
	double jtSubJetEta1, jtSubJetEta2, jtSubJetPhi1, jtSubJetPhi2;
	int hiBin;
	dataT->SetBranchAddress("jtpt",&jtpt);
	dataT->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	dataT->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	dataT->SetBranchAddress("jtSubJetEta1",&jtSubJetEta1);
	dataT->SetBranchAddress("jtSubJetEta2",&jtSubJetEta2);
	dataT->SetBranchAddress("jtSubJetPhi1",&jtSubJetPhi1);
	dataT->SetBranchAddress("jtSubJetPhi2",&jtSubJetPhi2);
	dataT->SetBranchAddress("minCSV",&minCSV);
	dataT->SetBranchAddress("minJP",&minJP);
	dataT->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	dataT->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	dataT->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	for(int ientry=0; ientry<dataT->GetEntries(); ientry++){
		dataT->GetEntry(ientry);
		
		if(jtpt<100) continue;
		if(jtSubJetdR<0.1) continue;
		
		int ibin=0;
		while(jtpt>jetPtBin[ibin+1]) ibin++;
		
		//M^2 = 2*pt1*pt2( cosh(deta) - cos(dphi) )
		double dphi = TMath::ACos(TMath::Cos(jtSubJetPhi1-jtSubJetPhi2));
		double invMass2 = 2*jtSubJetPt1*jtSubJetPt2*( TMath::CosH(jtSubJetEta1-jtSubJetEta2) - TMath::Cos(dphi));
		 
		if(minCSV>0.8 && jtSubJetSvtxm1>1.5 && jtSubJetSvtxm2>1.5 && jtSubJetdR>0.1) invMassDataGSP[ibin]->Fill(sqrt(invMass2));
		invMassDataIncl[ibin]->Fill(sqrt(invMass2));
	}
	
	
	//finish with MC
	TTree *mcT = (TTree*)fppmc->Get("gspTree");
	int jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
	double jtSubJetHadronDR1, jtSubJetHadronDR2, jtSubJetHadronPt1, jtSubJetHadronPt2;
	double weight;
	mcT->SetBranchAddress("jtpt",&jtpt);
	mcT->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	mcT->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	mcT->SetBranchAddress("jtSubJetEta1",&jtSubJetEta1);
	mcT->SetBranchAddress("jtSubJetEta2",&jtSubJetEta2);
	mcT->SetBranchAddress("jtSubJetPhi1",&jtSubJetPhi1);
	mcT->SetBranchAddress("jtSubJetPhi2",&jtSubJetPhi2);
	mcT->SetBranchAddress("minCSV",&minCSV);
	mcT->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	mcT->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	mcT->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
	mcT->SetBranchAddress("jtSubJetHadronDR1",&jtSubJetHadronDR1);
	mcT->SetBranchAddress("jtSubJetHadronDR2",&jtSubJetHadronDR2);
	mcT->SetBranchAddress("jtSubJetHadronPt1",&jtSubJetHadronPt1);
	mcT->SetBranchAddress("jtSubJetHadronPt2",&jtSubJetHadronPt2);
	mcT->SetBranchAddress("minJP",&minJP);
	mcT->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	mcT->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	mcT->SetBranchAddress("weight",&weight);
	
	for(int ientry=0; ientry<mcT->GetEntries(); ientry++){
		mcT->GetEntry(ientry);
		
		if(jtpt<100) continue;
	
		if(jtSubJetdR<0.1) continue;
		
		int ibin=0;
		while(jtpt>jetPtBin[ibin+1]) ibin++;
		
		//M^2 = 2*pt1*pt2( cosh(deta) - cos(dphi) )
		double dphi = TMath::ACos(TMath::Cos(jtSubJetPhi1-jtSubJetPhi2));
		double invMass2 = 2*jtSubJetPt1*jtSubJetPt2*( TMath::CosH(jtSubJetEta1-jtSubJetEta2) - TMath::Cos(dphi));
		 
		if(jtSubJetHadronFlavor1==5 && jtSubJetHadronFlavor2==5) invMassMCGSP[ibin]->Fill(sqrt(invMass2), weight);
		invMassMCIncl[ibin]->Fill(sqrt(invMass2), weight);
	}
	
	for(int i=0; i<ptBins; i++){
		invMassDataGSP[i]->Scale(1./invMassDataGSP[i]->Integral());
		invMassDataIncl[i]->Scale(1./invMassDataIncl[i]->Integral());
		invMassMCGSP[i]->Scale(1./invMassMCGSP[i]->Integral());
		invMassMCIncl[i]->Scale(1./invMassMCIncl[i]->Integral());
		
		invMassDataGSP[i]->Scale(1./invMassDataGSP[i]->GetBinWidth(1));
		invMassDataIncl[i]->Scale(1./invMassDataIncl[i]->GetBinWidth(1));
		invMassMCGSP[i]->Scale(1./invMassMCGSP[i]->GetBinWidth(1));
		invMassMCIncl[i]->Scale(1./invMassMCIncl[i]->GetBinWidth(1));
	}
	
	//now correct for tagging efficiencies
	TFile *fEff = new TFile("invMassMCTagEff.root");
	TH1D *mcEff[ptBins];
	TF1 *mcEffFits[ptBins];
	for(int i=0; i<ptBins; i++){
		mcEff[i] = (TH1D*)fEff->Get(Form("tagged_%d",i));
		mcEffFits[i] = new TF1(Form("mcEffFits_%d",i),"pol1",5,40);
		if(i==0) mcEff[i]->Fit(Form("mcEffFits_%d",i),"","",5,25);
		else mcEff[i]->Fit(Form("mcEffFits_%d",i),"","",8,40);
		
		for(int ibin=1; ibin<=invMassDataGSP[i]->GetNbinsX(); ibin++){
			invMassDataGSP[i]->SetBinContent(ibin, invMassDataGSP[i]->GetBinContent(ibin)/mcEffFits[i]->Eval(invMassDataGSP[i]->GetBinCenter(ibin)));
			invMassDataGSP[i]->SetBinError(ibin, invMassDataGSP[i]->GetBinError(ibin)/mcEffFits[i]->Eval(invMassDataGSP[i]->GetBinCenter(ibin)));		
		}
		
		invMassDataGSP[i]->Scale(1./invMassDataGSP[i]->Integral());		
		invMassDataGSP[i]->Scale(1./invMassDataGSP[i]->GetBinWidth(1));
	}
	
	
	TLegend *leg2 = new TLegend(0.217,0.6,0.62,0.775);
	TCanvas *cc = new TCanvas("cc","",1200,600);
	TLatex *tl1[ptBins];
	cc->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cc->cd(i+1);
		invMassDataIncl[i]->SetXTitle("Mass (GeV)");
		invMassDataIncl[i]->SetYTitle("1/N dN/dm (1/GeV)");
		invMassDataIncl[i]->SetMaximum(0.16);
		invMassDataIncl[i]->Draw();
		invMassDataGSP[i]->SetMarkerColor(kMagenta+2);
		invMassDataGSP[i]->SetLineColor(kMagenta+2);
		invMassDataGSP[i]->Draw("same");
		
		invMassMCGSP[i]->SetLineColor(4);
		invMassMCGSP[i]->SetMarkerColor(4);
		invMassMCGSP[i]->Draw("same");
		
		if(i<ptBins-1) tl1[i] = new TLatex(3,0.14,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		else tl1[i] = new TLatex(5,0.14,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		tl1[i]->Draw("Same");
		
		if(i==1){
			leg2->AddEntry(invMassDataIncl[i], "Inclusive pp", "lpf");
			leg2->AddEntry(invMassDataGSP[i], "pp g->bb", "lp");
			leg2->AddEntry(invMassMCGSP[i], "Pythia 8 g->bb", "lp");
			leg2->Draw("same");
		}
	}
}
	
	



