

#include "TH2D.h"
#include "TH1D.h"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"

#include <iostream>

//Need to load the RooUnfold Library first before compiling!!
//gSystem->Load("/Users/kjung/Downloads/RooUnfold/libRooUnfold.so");
#include "/Users/kjung/Downloads/RooUnfold/src/RooUnfold.h"
#include "/Users/kjung/Downloads/RooUnfold/src/RooUnfoldBayes.h"

using std::cout;
using std::endl;

const double PI=3.14159;

double getDphi(double phi1, double phi2){
	
	//project phi's to 0->pi
	phi1 = TMath::ACos(TMath::Cos(phi1));
	phi2 = TMath::ACos(TMath::Cos(phi2));
	
	return(std::abs(phi1-phi2));
}

void rescaleAj(TH1D *out, TH2D* preUnfold, TH2D* postUnfold, double ptLow, double ptHigh){
	
	TH2D *temp = (TH2D*)postUnfold->Clone("temp");
	temp->Divide(preUnfold);
	
	TH1D *temp1D = (TH1D*)out->Clone("temp1D");
	temp1D->Reset();
	
	for(int ixbin=1; ixbin<=temp->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=temp->GetYaxis()->FindBin(temp->GetXaxis()->GetBinLowEdge(ixbin)); iybin++){
			if(temp->GetXaxis()->GetBinCenter(ixbin)>ptLow && temp->GetXaxis()->GetBinCenter(ixbin)<ptHigh){
				
				//for each entry in the 2D histogram, mock filling the 1D pt2/pt1 histogram N number of times
				temp1D->Fill(temp->GetYaxis()->GetBinCenter(iybin)/(temp->GetXaxis()->GetBinCenter(ixbin)+temp->GetYaxis()->GetBinCenter(ixbin)), temp->GetBinContent(ixbin,iybin));
			}
		}
	}
	out->Add(temp1D);
	temp1D->Delete();
}

void fillAj(TH1D *out, TH2D* unfold, double ptLow, double ptHigh){
	for(int ixbin=1; ixbin<=unfold->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=unfold->GetYaxis()->FindBin(unfold->GetXaxis()->GetBinLowEdge(ixbin)); iybin++){
			if(unfold->GetXaxis()->GetBinCenter(ixbin)>ptLow && unfold->GetXaxis()->GetBinCenter(ixbin)<ptHigh){
				
				//for each entry in the 2D histogram, mock filling the 1D pt2/pt1 histogram N number of times
				out->Fill(unfold->GetYaxis()->GetBinCenter(iybin)/unfold->GetXaxis()->GetBinCenter(ixbin), unfold->GetBinContent(ixbin,iybin));
			}
		}
	}
}

void normalizeHisto(TH1D *histo){
	for(int i=1; i<=histo->GetNbinsX(); i++){
		histo->SetBinContent(i, histo->GetBinContent(i)/histo->GetBinWidth(i));
		histo->SetBinError(i, histo->GetBinError(i)/histo->GetBinWidth(i));
	}
	histo->Scale(1./histo->Integral());
}

void test2Dunfold(){
		
	//switch to fill the reco distribution from MC for testing...
	bool fillMC = true; 
	bool doSymmeterize = true;
	bool doBTagging = false;
	
	const int nPtBins = 2;
	const double ptBins[nPtBins+1] = {100, 140, 9999};
	const double recoBins[31] = {0, 5, 8, 10, 12, 15, 17.5, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 150, 170, 190, 205, 220, 250, 350, 500};
	/*const int nGenBins = 40;
	double genBins[nGenBins+1];
	genBins[0] = 30;
	for(int i=1; i<nGenBins+1; i++){ genBins[i] = 30. + (double)i*(470./(double)nGenBins); }*/
	
	TH1D *leadingSubjetPt[nPtBins], *leadingSubjetRefPt[nPtBins];
	TH1D *subleadingSubjetPt[nPtBins], *subleadingSubjetRefPt[nPtBins];
	RooUnfoldResponse *leadingResponse[nPtBins];
	RooUnfoldResponse *subleadingResponse[nPtBins];
	for(int i=0; i<nPtBins; i++){
		leadingSubjetPt[i] = new TH1D(Form("leadingSubjetPt_%d",i),"",30,recoBins); leadingSubjetPt[i]->Sumw2();
		subleadingSubjetPt[i] = new TH1D(Form("subleadingSubjetPt_%d",i),"",30,recoBins); subleadingSubjetPt[i]->Sumw2();
		
		leadingSubjetRefPt[i] = new TH1D(Form("leadingSubjetRefPt_%d",i),"",30,recoBins); leadingSubjetRefPt[i]->Sumw2();
		subleadingSubjetRefPt[i] = new TH1D(Form("subleadingSubjetRefPt_%d",i),"",30,recoBins); subleadingSubjetRefPt[i]->Sumw2();
		
		leadingResponse[i] = new RooUnfoldResponse(leadingSubjetPt[i], leadingSubjetRefPt[i]);
		subleadingResponse[i] = new RooUnfoldResponse(subleadingSubjetPt[i], subleadingSubjetRefPt[i]);
	}
		
	
	TH2D *jetReco2D[nPtBins], *jetPrior2D[nPtBins];
	RooUnfoldResponse *resp[nPtBins];
	for(int i=0; i<nPtBins; i++){
		jetReco2D[i] = new TH2D(Form("jetReco2D_%d",i),"",30,recoBins,30,recoBins); jetReco2D[i]->Sumw2();
		jetPrior2D[i] = new TH2D(Form("jetPrior2D_%d",i),"",30,recoBins,30,recoBins); jetPrior2D[i]->Sumw2();
		resp[i] = new RooUnfoldResponse(jetReco2D[i], jetPrior2D[i]);
	}
	
	TH1D *recoZg[nPtBins], *genZg[nPtBins], *unfoldedZg[nPtBins];
	for(int i=0; i<nPtBins; i++){
		recoZg[i] = new TH1D(Form("recoXj_pt%d",i),"",10,0,0.5); recoZg[i]->Sumw2();
		genZg[i] = new TH1D(Form("genXj_pt%d",i),"",10,0,0.5); genZg[i]->Sumw2();
		unfoldedZg[i] = new TH1D(Form("unfoldedZg_pt%d",i),"",10,0,0.5); unfoldedZg[i]->Sumw2();

	}
	
	TFile *fin = new TFile("gspTreeOut_ppMC_muFilteredMC_Pythia8_fullSample.root");
	TTree *tin = (TTree*)fin->Get("gspTree");
	
	TFile *finData = new TFile("gspTreeOut_ppData_pt60cut_explicitJTA_recalibJP_muFiltered_perEvent.root");
	TTree *tinData = (TTree*)finData->Get("gspTree");
	
	double jtpt, subJetPt1, subJetPt2, subJetRefPt1, subJetRefPt2;
	double jtSubJetSvtxm1, jtSubJetSvtxm2, jtSubJetcsvV2_1, jtSubJetcsvV2_2;
	double weight;
	int hiBin;
	tin->SetBranchAddress("jtpt",&jtpt);
	tin->SetBranchAddress("jtSubJetPt1",&subJetPt1);
	tin->SetBranchAddress("jtSubJetPt2",&subJetPt2);
	tin->SetBranchAddress("jtSubJetRefPt1",&subJetRefPt1);
	tin->SetBranchAddress("jtSubJetRefPt2",&subJetRefPt2);
	tin->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
	tin->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
	tin->SetBranchAddress("jtSubJetcsvV2_1",&jtSubJetcsvV2_1);
	tin->SetBranchAddress("jtSubJetcsvV2_2",&jtSubJetcsvV2_2);
	tin->SetBranchAddress("weight",&weight);
	
	cout << "starting MC" << endl;
	for(int ientry=0; ientry<tin->GetEntries(); ientry++){
		
		tin->GetEntry(ientry);
		
		if(jtpt<100) continue;
		//if(subJetPt1<0 || subJetRefPt1<0 || subJetPt2<0 || subJetRefPt2<0) continue;
		if(doBTagging){
			if(jtSubJetSvtxm1<1.5 || jtSubJetSvtxm2<1.5 || jtSubJetcsvV2_1<0.8 || jtSubJetcsvV2_2<0.8) continue;
		}
		
		int ptBin=0;
		while(jtpt>ptBins[ptBin+1]) ptBin++;
		
		leadingSubjetRefPt[ptBin]->Fill(subJetRefPt1, weight);
		subleadingSubjetRefPt[ptBin]->Fill(subJetRefPt2, weight);
		if(fillMC){
			leadingSubjetPt[ptBin]->Fill(subJetPt1, weight);
			subleadingSubjetPt[ptBin]->Fill(subJetPt2, weight);
		}
		leadingResponse[ptBin]->Fill(subJetPt1, subJetRefPt1, weight);
		subleadingResponse[ptBin]->Fill(subJetPt2, subJetRefPt2, weight);

		if(fillMC){
			jetReco2D[ptBin]->Fill(subJetPt1,subJetPt2, weight);
			if(doSymmeterize) jetReco2D[ptBin]->Fill(subJetPt2,subJetPt1,weight);
		}
		
		jetPrior2D[ptBin]->Fill(subJetRefPt1,subJetRefPt2,weight);
		if(doSymmeterize) jetPrior2D[ptBin]->Fill(subJetRefPt2,subJetRefPt1,weight);
		
		genZg[ptBin]->Fill(subJetRefPt2/(subJetRefPt1+subJetRefPt2),weight);
		if(fillMC) recoZg[ptBin]->Fill(subJetPt2/(subJetPt1+subJetPt2),weight);
		
		resp[ptBin]->Fill(subJetPt1, subJetPt2, subJetRefPt1, subJetRefPt2,weight);
		if(doSymmeterize) resp[ptBin]->Fill(subJetPt2, subJetPt1, subJetRefPt2, subJetRefPt1,weight);
	
	}
	
	cout << "starting data" << endl;
	tinData->SetBranchAddress("jtpt",&jtpt);
	tinData->SetBranchAddress("jtSubJetPt1",&subJetPt1);
	tinData->SetBranchAddress("jtSubJetPt2",&subJetPt2);
	if(!fillMC){
		for(int ientry=0; ientry<tinData->GetEntries(); ientry++){
			
			tinData->GetEntry(ientry);
			
			if(jtpt<100) continue;
			if(subJetPt1<0 || subJetPt2<0) continue;
			
			int ptBin=0;
			while(jtpt>ptBins[ptBin+1]) ptBin++;
	
			jetReco2D[ptBin]->Fill(subJetPt1,subJetPt2);
			if(doSymmeterize) jetReco2D[ptBin]->Fill(subJetPt2,subJetPt1);
			recoZg[ptBin]->Fill(subJetPt2/(subJetPt1+subJetPt2));
			
			leadingSubjetPt[ptBin]->Fill(subJetPt1, weight);
			subleadingSubjetPt[ptBin]->Fill(subJetPt2, weight);
			
		}
	}
	
	if(doSymmeterize){
		for(int ibin=0; ibin<nPtBins; ibin++){
			jetReco2D[ibin]->Scale(0.5);
			jetPrior2D[ibin]->Scale(0.5);
		}
	}
	
	//do the unfolding
	
	RooUnfoldBayes *ruLeading[nPtBins];
	RooUnfoldBayes *ruSubleading[nPtBins];
	
	RooUnfoldBayes *ru20[nPtBins];
	TH2D *unfolded[nPtBins];
	TH1D *leadingUnfolded[nPtBins], *subleadingUnfolded[nPtBins];
	for(int i=0; i<nPtBins; i++){
		cout << "starting unfold, ptBin " << i << endl;
		
		ruLeading[i] = new RooUnfoldBayes(leadingResponse[i], leadingSubjetPt[i], 10);
		ruSubleading[i] = new RooUnfoldBayes(subleadingResponse[i], subleadingSubjetPt[i], 10);
		
		leadingUnfolded[i] = (TH1D*)ruLeading[i]->Hreco(RooUnfold::kCovToy);
		subleadingUnfolded[i] = (TH1D*)ruSubleading[i]->Hreco(RooUnfold::kCovToy);
		
		ru20[i] = new RooUnfoldBayes(resp[i], jetReco2D[i], 5);
		ru20[i]->SetVerbose(1);
		unfolded[i] = (TH2D*)ru20[i]->Hreco(RooUnfold::kCovToy);
		
		fillAj(unfoldedZg[i], unfolded[i], 0, 1000);
	}
	
	//unfolding done
	
	TCanvas *c0 = new TCanvas("c0","",1200,600);
	c0->Divide(nPtBins,1);
	for(int i=0; i<nPtBins; i++){
		c0->cd(i+1);
		leadingSubjetPt[i]->SetXTitle("Leading Subjet pT");
		normalizeHisto(leadingSubjetPt[i]);
		leadingSubjetPt[i]->Draw();
		leadingSubjetRefPt[i]->SetMarkerColor(2);
		leadingSubjetRefPt[i]->SetLineColor(2);
		normalizeHisto(leadingSubjetRefPt[i]);
		leadingSubjetRefPt[i]->Draw("Same");
		leadingUnfolded[i]->SetMarkerColor(4);
		leadingUnfolded[i]->SetLineColor(4);
		normalizeHisto(leadingUnfolded[i]);
		leadingUnfolded[i]->Draw("same");
	}
	
	TCanvas *c02 = new TCanvas("c02","",1200,600);
	c02->Divide(nPtBins,1);
	for(int i=0; i<nPtBins; i++){
		c02->cd(i+1);
		subleadingSubjetPt[i]->SetXTitle("Subleading Subjet pT");
		subleadingSubjetPt[i]->Draw();
		normalizeHisto(subleadingSubjetPt[i]);
		subleadingSubjetRefPt[i]->SetMarkerColor(2);
		subleadingSubjetRefPt[i]->SetLineColor(2);
		normalizeHisto(subleadingSubjetRefPt[i]);
		subleadingSubjetRefPt[i]->Draw("Same");
		subleadingUnfolded[i]->SetMarkerColor(4);
		subleadingUnfolded[i]->SetLineColor(4);
		normalizeHisto(subleadingUnfolded[i]);
		subleadingUnfolded[i]->Draw("same");
	}
	
		
	
	TCanvas *c1 = new TCanvas("c1","",800,800);
	c1->cd();
			
	//jetPrior2D[0]->Divide(unfolded[0]);
	if(doSymmeterize){
		jetPrior2D[0]->SetXTitle("Gen, p_{T,1+2}");
		jetPrior2D[0]->SetYTitle("Gen, p_{T,2+1}");
	}
	else{
		jetPrior2D[0]->SetXTitle("Gen, p_{T,1}");
		jetPrior2D[0]->SetYTitle("Gen, p_{T,2}");
	}
	jetPrior2D[0]->Draw("colz");
	
	TCanvas *c2 = new TCanvas("c2","",1200,600);
	c2->Divide(nPtBins,1);
	for(int i=0; i<nPtBins; i++){
		c2->cd(i+1);
		recoZg[i]->Scale(1./recoZg[i]->GetBinWidth(1));
		recoZg[i]->Scale(1./recoZg[i]->Integral());
		recoZg[i]->Draw();

		genZg[i]->SetLineColor(2);
		genZg[i]->SetMarkerColor(2);
		genZg[i]->Scale(1./genZg[i]->GetBinWidth(1));
		genZg[i]->Scale(1./genZg[i]->Integral());
		genZg[i]->Draw("same");

		unfoldedZg[i]->SetLineColor(4);
		unfoldedZg[i]->SetMarkerColor(4);
		unfoldedZg[i]->Scale(1./unfoldedZg[i]->GetBinWidth(1));
		unfoldedZg[i]->Scale(1./unfoldedZg[i]->Integral());
		unfoldedZg[i]->Draw("same");
	}
	
	/*TCanvas *c3 = new TCanvas("c3","",800,800);
	c3->Divide(2,2);
	for(int i=0; i<3; i++){
		c3->cd(i+1);
		
		unfoldedXj_20[i]->SetLineColor(4);
		unfoldedXj_20[i]->SetMarkerColor(4);
		//unfoldedXj_20[i]->Scale(1./unfoldedXj_20[i]->Integral(1,hiBin)/unfoldedXj_20[i]->GetBinWidth(3));
		//unfoldedXj_20[i]->Draw("same");
		
		unfoldedXj_10[i]->SetLineColor(kBlue-1);
		unfoldedXj_10[i]->SetMarkerColor(kBlue-1);
		unfoldedXj_10[i]->Scale(1./unfoldedXj_10[i]->Integral(1,hiBin,"width"));
		//unfoldedXj_10[i]->Draw("same");
		
		unfoldedXj_5[i]->SetLineColor(kBlue-3);
		unfoldedXj_5[i]->SetMarkerColor(kBlue-3);
		unfoldedXj_5[i]->Scale(1./unfoldedXj_5[i]->Integral(1,hiBin,"width"));
		//unfoldedXj_5[i]->Draw("same");
	}*/
}