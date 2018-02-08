

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TSystem.h"
#include "headers.h"
#include <iostream>

using std::cout;
using std::endl;

void setupCFIT(double taggerVal=0.9){

	//const int ptBins = 5;
	//const double xbins[ptBins+1] = {80, 100, 120, 170, 250, 9999};
	const bool ispp = true;

	//TFile *fin_data = new TFile("gspTreeOut_ppData_explicitJTA_pt60.root");
	TFile *fin_data = new TFile("/mnt/hadoop/store/user/kjung/bJetSplittingTrees/gspTreeOut_ppData_explicitJTA_JPRecalibration_pt60.root");
        //TFile *fin_mc = new TFile("gspTreeOut_ppMC_Pythia8_explicitJTA_allPthat.root");
	TFile *fin_mc = new TFile("gspTreeOut_ppMC_Pythia8_explicitJTA_recalibJP_withNegCSV.root");
	TFile *fin_b_mc = new TFile("gspTreeOut_ppMC_Pythia6_bjetMC_explicitJTA_recalibJP_withNegCSV.root");

	TTree *tdata = (TTree*)fin_data->Get("gspTree");
	TTree *tmc = (TTree*)fin_mc->Get("gspTree");
	TTree *tbmc = (TTree*)fin_b_mc->Get("gspTree");
	
	TH1D *h_data[ptBins], *h_data_tag[ptBins], *h_data_untag[ptBins];
	TH1D *h_usdg_mc[ptBins], *h_c_mc[ptBins], *h_b_mc[ptBins];
	TH1D *h_usdg_mc_tag[ptBins], *h_c_mc_tag[ptBins], *h_b_mc_tag[ptBins];
	TH1D *h_usdg_mc_untag[ptBins], *h_c_mc_untag[ptBins], *h_b_mc_untag[ptBins];

	TH1D *h_usdg_mc_negTag[ptBins], *h_c_mc_negTag[ptBins], *h_b_mc_negTag[ptBins];
	TH1D *h_usdg_mc_negTagCut[ptBins], *h_c_mc_negTagCut[ptBins], *h_b_mc_negTagCut[ptBins];
	TH1D *h_data_negTag[ptBins];
	TH1D *h_data_negTagCut[ptBins];

	for(int i=0; i<ptBins; i++){
		h_data[i] = new TH1D(Form("h_data_%d",i),"",40,0,3); h_data[i]->Sumw2();
		h_data_tag[i] = new TH1D(Form("h_data_tag_%d",i),"",40,0,3); h_data_tag[i]->Sumw2();
		h_data_untag[i] = new TH1D(Form("h_data_untag_%d",i),"",40,0,3); h_data_untag[i]->Sumw2();
		
		h_usdg_mc[i] = new TH1D(Form("h_usdg_mc_%d",i),"",40,0,3); h_usdg_mc[i]->Sumw2();
		h_c_mc[i] = new TH1D(Form("h_c_mc_%d",i),"",40,0,3); h_c_mc[i]->Sumw2();
		h_b_mc[i] = new TH1D(Form("h_b_mc_%d",i),"",40,0,3); h_b_mc[i]->Sumw2();
		
		h_usdg_mc_tag[i] = new TH1D(Form("h_usdg_mc_tag_%d",i),"",40,0,3); h_usdg_mc_tag[i]->Sumw2();
		h_c_mc_tag[i] = new TH1D(Form("h_c_mc_tag_%d",i),"",40,0,3); h_c_mc_tag[i]->Sumw2();
		h_b_mc_tag[i] = new TH1D(Form("h_b_mc_tag_%d",i),"",40,0,3); h_b_mc_tag[i]->Sumw2();
		
		h_usdg_mc_untag[i] = new TH1D(Form("h_usdg_mc_untag_%d",i),"",40,0,3); h_usdg_mc_untag[i]->Sumw2();
		h_c_mc_untag[i] = new TH1D(Form("h_c_mc_untag_%d",i),"",40,0,3); h_c_mc_untag[i]->Sumw2();
		h_b_mc_untag[i] = new TH1D(Form("h_b_mc_untag_%d",i),"",40,0,3); h_b_mc_untag[i]->Sumw2();
		
		h_usdg_mc_negTag[i] = new TH1D(Form("h_usdg_mc_negTag_%d",i),"",30,-10,3); h_usdg_mc_negTag[i]->Sumw2();
		h_c_mc_negTag[i] = new TH1D(Form("h_c_mc_negTag_%d",i),"",30,-10,3); h_c_mc_negTag[i]->Sumw2();
		h_b_mc_negTag[i] = new TH1D(Form("h_b_mc_negTag_%d",i),"",30,-10,3); h_b_mc_negTag[i]->Sumw2();

		h_usdg_mc_negTagCut[i] = new TH1D(Form("h_usdg_mc_negTagCut_%d",i),"",30,-10,3); h_usdg_mc_negTagCut[i]->Sumw2();
		h_c_mc_negTagCut[i] = new TH1D(Form("h_c_mc_negTagCut_%d",i),"",30,-10,3); h_c_mc_negTagCut[i]->Sumw2();
		h_b_mc_negTagCut[i] = new TH1D(Form("h_b_mc_negTagCut_%d",i),"",30,-10,3); h_b_mc_negTagCut[i]->Sumw2();

		h_data_negTag[i] = new TH1D(Form("h_data_negTag_%d",i),"",30,-10,3); h_data_negTag[i]->Sumw2();
		h_data_negTagCut[i] = new TH1D(Form("h_data_negTagCut_%d",i),"",30,-10,3); h_data_negTagCut[i]->Sumw2();
	}
	
	double jtpt, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, minJP, discr_csvV2, ndiscr_csvV2, discr_prob;
	double jtSubJetJP1, jtSubJetJP2, jtSubJetCSV1, jtSubJetCSV2, jtSubJetNegCSV1, jtSubJetNegCSV2;
	int hiBin;
	double weight;
	
	//start with data
	tdata->SetBranchAddress("jtpt",&jtpt);
	tdata->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tdata->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tdata->SetBranchAddress("minCSV",&minCSV);
	tdata->SetBranchAddress("discr_prob", &discr_prob);
	tdata->SetBranchAddress("minJP",&minJP);
	tdata->SetBranchAddress("ndiscr_csvV2",&ndiscr_csvV2);
	tdata->SetBranchAddress("discr_csvV2",&discr_csvV2);
	tdata->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	tdata->SetBranchAddress("jtSubJetJP_1",&jtSubJetJP1);
	tdata->SetBranchAddress("jtSubJetJP_2",&jtSubJetJP2);
	tdata->SetBranchAddress("jtSubJetcsvV2_1",&jtSubJetCSV1);
	tdata->SetBranchAddress("jtSubJetcsvV2_2",&jtSubJetCSV2);
	tdata->SetBranchAddress("jtSubJetNegCsvV2_1",&jtSubJetNegCSV1);
	tdata->SetBranchAddress("jtSubJetNegCsvV2_2",&jtSubJetNegCSV2);

	tdata->SetBranchAddress("weight",&weight);
	
	if(!ispp) tdata->SetBranchAddress("hiBin",&hiBin);
	for(int ievt=0; ievt<tdata->GetEntries(); ievt++){
		tdata->GetEntry(ievt);
		if((!doSubjets && jtpt>80) || (doSubjets && jtSubJetPt1>30 && jtSubJetPt2>30)){

                        //weight*=1e11;
			int ptBin=0;
			while(jtpt>xbins[ptBin+1] && ptBin<ptBins-1) ptBin++;
			
			int ptBin1=0, ptBin2=0;
			while(jtSubJetPt1>xbins[ptBin1+1] && ptBin1<ptBins-1) ptBin1++;
			while(jtSubJetPt2>xbins[ptBin2+1] && ptBin2<ptBins-1) ptBin2++;

			if(!doSubjets){
				h_data[ptBin]->Fill(discr_prob, weight);
				if(discr_csvV2>taggerVal) h_data_tag[ptBin]->Fill(discr_prob, weight);
				else h_data_untag[ptBin]->Fill(discr_prob, weight);
				
				if(ndiscr_csvV2>0) h_data_negTag[ptBin]->Fill(discr_prob, weight);
				if(ndiscr_csvV2>taggerVal) h_data_negTagCut[ptBin]->Fill(discr_prob, weight); 
			}
			else{
				h_data[ptBin1]->Fill(jtSubJetJP1, weight);
				h_data[ptBin2]->Fill(jtSubJetJP2, weight);
				if(jtSubJetCSV1>taggerVal) h_data_tag[ptBin1]->Fill(jtSubJetJP1, weight);
				else h_data_untag[ptBin1]->Fill(jtSubJetJP1, weight);
				if(jtSubJetCSV2>taggerVal) h_data_tag[ptBin2]->Fill(jtSubJetJP2, weight);
				else h_data_untag[ptBin2]->Fill(jtSubJetJP2, weight);
				
				if(jtSubJetNegCSV1>0) h_data_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
				if(jtSubJetNegCSV1>taggerVal) h_data_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);
				if(jtSubJetNegCSV2>0) h_data_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
				if(jtSubJetNegCSV2>taggerVal) h_data_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);

			}
		}	
	}
	
	//then do MC
	int jtHadronFlavor, jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
	tmc->SetBranchAddress("jtpt",&jtpt);
	tmc->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tmc->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tmc->SetBranchAddress("minCSV",&minCSV);
	tmc->SetBranchAddress("minJP",&minJP);
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
	if(!ispp) tmc->SetBranchAddress("hiBin",&hiBin);
	for(int ievt=0; ievt<tmc->GetEntries(); ievt++){
		tmc->GetEntry(ievt);
		if((!doSubjets && jtpt>80) || (doSubjets && jtSubJetPt1>30 && jtSubJetPt2>30)){

                        //weight*=1e11;
			
			int ptBin=0;
			while(jtpt>xbins[ptBin+1] && ptBin<ptBins-1) ptBin++;

			int ptBin1=0, ptBin2=0;
			while(jtSubJetPt1>xbins[ptBin1+1] && ptBin1<ptBins-1) ptBin1++;
			while(jtSubJetPt2>xbins[ptBin2+1] && ptBin2<ptBins-1) ptBin2++;

			if(!doSubjets){
				/*if(abs(jtHadronFlavor)==4){
                              h_c_mc[ptBin]->Fill(discr_prob, weight);
                              if(minCSV>taggerVal) h_c_mc_tag[ptBin]->Fill(discr_prob, weight);
                              else h_c_mc_untag[ptBin]->Fill(discr_prob, weight);
                              }*/
				//else{
					h_usdg_mc[ptBin]->Fill(discr_prob, weight);
					if(discr_csvV2>taggerVal) h_usdg_mc_tag[ptBin]->Fill(discr_prob, weight);
					else h_usdg_mc_untag[ptBin]->Fill(discr_prob, weight);
					
					if(ndiscr_csvV2>0) h_usdg_mc_negTag[ptBin]->Fill(discr_prob, weight);
					if(ndiscr_csvV2>taggerVal) h_usdg_mc_negTagCut[ptBin]->Fill(discr_prob, weight);
				//}
			}
			else{
					h_usdg_mc[ptBin1]->Fill(jtSubJetJP1, weight);
					h_usdg_mc[ptBin2]->Fill(jtSubJetJP2, weight);
					if(jtSubJetCSV1>taggerVal) h_usdg_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
					else h_usdg_mc_untag[ptBin1]->Fill(jtSubJetJP1, weight);
					if(jtSubJetCSV2>taggerVal) h_usdg_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
					else h_usdg_mc_untag[ptBin2]->Fill(jtSubJetJP2, weight);   
					
					if(jtSubJetNegCSV1>0) h_usdg_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
					if(jtSubJetNegCSV1>taggerVal) h_usdg_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);
					if(jtSubJetNegCSV2>0) h_usdg_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
					if(jtSubJetNegCSV2>taggerVal) h_usdg_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);
			}	
		}
	}
	
	
	//then do bMC
	tbmc->SetBranchAddress("jtpt",&jtpt);
	tbmc->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tbmc->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tbmc->SetBranchAddress("minCSV",&minCSV);
	tbmc->SetBranchAddress("minJP",&minJP);
	tbmc->SetBranchAddress("discr_prob", &discr_prob);
	tbmc->SetBranchAddress("ndiscr_csvV2",&ndiscr_csvV2);
	tbmc->SetBranchAddress("discr_csvV2",&discr_csvV2);
	tbmc->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	tbmc->SetBranchAddress("jtHadronFlavor",&jtHadronFlavor);
	tbmc->SetBranchAddress("weight",&weight);
	tbmc->SetBranchAddress("jtSubJetJP_1",&jtSubJetJP1);
	tbmc->SetBranchAddress("jtSubJetJP_2",&jtSubJetJP2);
	tbmc->SetBranchAddress("jtSubJetcsvV2_1",&jtSubJetCSV1);
	tbmc->SetBranchAddress("jtSubJetcsvV2_2",&jtSubJetCSV2);
	tbmc->SetBranchAddress("jtSubJetNegCsvV2_1",&jtSubJetNegCSV1);
	tbmc->SetBranchAddress("jtSubJetNegCsvV2_2",&jtSubJetNegCSV2);
	tbmc->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	tbmc->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
	if(!ispp) tbmc->SetBranchAddress("hiBin",&hiBin);
	for(int ievt=0; ievt<tbmc->GetEntries(); ievt++){
		tbmc->GetEntry(ievt);
		if((!doSubjets && jtpt>80) || (doSubjets && jtSubJetPt1>30 && jtSubJetPt2>30)){

                        //weight*=1e11;
			
			int ptBin=0;
			while(jtpt>xbins[ptBin+1] && ptBin<ptBins-1) ptBin++;

			int ptBin1=0, ptBin2=0;
			while(jtSubJetPt1>xbins[ptBin1+1] && ptBin1<ptBins-1) ptBin1++;
			while(jtSubJetPt2>xbins[ptBin2+1] && ptBin2<ptBins-1) ptBin2++;

			if(!doSubjets){
				if(abs(jtHadronFlavor)==5){
					h_b_mc[ptBin]->Fill(discr_prob, weight);
					if(discr_csvV2>taggerVal) h_b_mc_tag[ptBin]->Fill(discr_prob, weight);
					else h_b_mc_untag[ptBin]->Fill(discr_prob, weight);
					
					if(ndiscr_csvV2>0) h_b_mc_negTag[ptBin]->Fill(discr_prob, weight);
					if(ndiscr_csvV2>taggerVal) h_b_mc_negTagCut[ptBin]->Fill(discr_prob, weight);
				}
			}
			else{
				if(abs(jtHadronFlavor)==5){
					h_b_mc[ptBin1]->Fill(jtSubJetJP1, weight);
					h_b_mc[ptBin2]->Fill(jtSubJetJP2, weight);
					if(jtSubJetCSV1>taggerVal) h_b_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
					else h_b_mc_untag[ptBin1]->Fill(jtSubJetCSV1, weight);
					if(jtSubJetCSV2>taggerVal) h_b_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
					else h_b_mc_untag[ptBin2]->Fill(jtSubJetCSV2, weight);

					if(jtSubJetNegCSV1>0) h_b_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
					if(jtSubJetNegCSV1>taggerVal) h_b_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);
					if(jtSubJetNegCSV2>0) h_b_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
					if(jtSubJetNegCSV2>taggerVal) h_b_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);
				}
			}	
		}
	}
	
	string outputName = "cfitInput_fullJets.root";
	if(doSubjets) outputName = "cfitInput_subjets.root";
	
	TFile *out = new TFile(outputName.c_str(),"recreate");
	out->cd();
	for(int i=0; i<ptBins; i++){
		h_data[i]->Write();
		h_data_tag[i]->Write();
		h_data_untag[i]->Write();
		
		h_usdg_mc[i]->Write();
		h_c_mc[i]->Write();
		h_b_mc[i]->Write();
		
		h_usdg_mc_tag[i]->Write();
		h_c_mc_tag[i]->Write();
		h_b_mc_tag[i]->Write();
		
		h_usdg_mc_untag[i]->Write();
		h_c_mc_untag[i]->Write();
		h_b_mc_untag[i]->Write();
		
		h_usdg_mc_negTag[i]->Write();
		h_c_mc_negTag[i]->Write();
		h_b_mc_negTag[i]->Write();

		h_usdg_mc_negTagCut[i]->Write();
		h_c_mc_negTagCut[i]->Write();
		h_b_mc_negTagCut[i]->Write();
		
		h_data_negTag[i]->Write();
		h_data_negTagCut[i]->Write();
	}
	out->Close();
}
