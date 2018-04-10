

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TSystem.h"
#include "headers.h"
#include <iostream>

using std::cout;
using std::endl;

double getZg(double j1, double j2){
        return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}

void setupCFIT(double taggerVal=0.8, double drSelection=-1e99, double svtxCut=1.2){

	//const int ptBins = 5;
	//const double xbins[ptBins+1] = {80, 100, 120, 170, 250, 9999};
	const bool ispp = true;
        const double jtptCutMin = 120;
        const double jtptCutMax = 400;

        const double switchToDR = false;

	TFile *fin_data = new TFile("gspTreeOut_ppData_muFilterPerJet_muPt4.root");
	//TFile *fin_data = new TFile("gspTreeOut_ppData_pt60cut_explicitJTA_recalibJP_muFiltered_perEvent.root");
	//TFile *fin_data = new TFile("gspTreeOut_ppData_noMuFilter.root");
        TFile *fin_mc = new TFile("gspTreeOut_ppMC_pythia8MuFilteredInJet_muPt4.root");
	//TFile *fin_mc = new TFile("gspTreeOut_pthat100_bbFiltered_dr0p1_4.root");
        //TFile *fin_mc = new TFile("gspTreeOut_pt50BHadron_muPt4Filtered.root");
        //TFile *fin_mc = new TFile("gspTreeOut_ppMC_QCDJet_noMuFilter_withRefInfo.root");
        //TFile *fin_light_mc = new TFile("gspTreeOut_ppMC_QCDJet_noMuFilter_withRefInfo.root");
        TFile *fin_light_mc = new TFile("gspTreeOut_ppMC_qcdPythia8_fullSample_withFullJets.root");
        //TFile *fin_b_mc = new TFile("gspTreeOut_ppMC_Pythia6_bjetMC_explicitJTA_recalibJP_withNegCSV.root");

	TTree *tdata = (TTree*)fin_data->Get("gspTree");
	TTree *tmc = (TTree*)fin_mc->Get("gspTree");
	TTree *tlightmc = (TTree*)fin_light_mc->Get("gspTree");
	
	TH1D *h_data[ptBins], *h_data_tag[ptBins], *h_data_untag[ptBins];
	TH1D *h_usdg_mc[ptBins], *h_c_mc[ptBins], *h_b_mc[ptBins];
	TH1D *h_usdg_mc_tag[ptBins], *h_c_mc_tag[ptBins], *h_b_mc_tag[ptBins];
	TH1D *h_usdg_mc_untag[ptBins], *h_c_mc_untag[ptBins], *h_b_mc_untag[ptBins];

	TH1D *h_usdg_mc_negTag[ptBins], *h_c_mc_negTag[ptBins], *h_b_mc_negTag[ptBins];
	TH1D *h_usdg_mc_negTagCut[ptBins], *h_c_mc_negTagCut[ptBins], *h_b_mc_negTagCut[ptBins];
	TH1D *h_data_negTag[ptBins];
	TH1D *h_data_negTagCut[ptBins];

        TH1D *h_data_zg[zgBins], *h_data_tag_zg[zgBins], *h_data_untag_zg[zgBins];
        TH1D *h_usdg_mc_zg[zgBins], *h_c_mc_zg[zgBins], *h_b_mc_zg[zgBins];
        TH1D *h_usdg_mc_tag_zg[zgBins], *h_c_mc_tag_zg[zgBins], *h_b_mc_tag_zg[zgBins];
        TH1D *h_usdg_mc_untag_zg[zgBins], *h_c_mc_untag_zg[zgBins], *h_b_mc_untag_zg[zgBins];

        TH1D *h_usdg_mc_negTag_zg[zgBins], *h_c_mc_negTag_zg[zgBins], *h_b_mc_negTag_zg[zgBins];
        TH1D *h_usdg_mc_negTagCut_zg[zgBins], *h_c_mc_negTagCut_zg[zgBins], *h_b_mc_negTagCut_zg[zgBins];
        TH1D *h_data_negTag_zg[zgBins];
        TH1D *h_data_negTagCut_zg[zgBins];

        for(int i=0; i<zgBins; i++){
            h_data_zg[i] = new TH1D(Form("h_data_zg_%d",i),"",40,0,3); h_data_zg[i]->Sumw2();
            h_data_tag_zg[i] = new TH1D(Form("h_data_tag_zg_%d",i),"",40,0,3); h_data_tag_zg[i]->Sumw2();
            h_data_untag_zg[i] = new TH1D(Form("h_data_untag_zg_%d",i),"",40,0,3); h_data_untag_zg[i]->Sumw2();

            h_usdg_mc_zg[i] = new TH1D(Form("h_usdg_mc_zg_%d",i),"",40,0,3); h_usdg_mc_zg[i]->Sumw2();
            h_c_mc_zg[i] = new TH1D(Form("h_c_mc_zg_%d",i),"",40,0,3); h_c_mc_zg[i]->Sumw2();
            h_b_mc_zg[i] = new TH1D(Form("h_b_mc_zg_%d",i),"",40,0,3); h_b_mc_zg[i]->Sumw2();

            h_usdg_mc_tag_zg[i] = new TH1D(Form("h_usdg_mc_tag_zg_%d",i),"",40,0,3); h_usdg_mc_tag_zg[i]->Sumw2();
            h_c_mc_tag_zg[i] = new TH1D(Form("h_c_mc_tag_zg_%d",i),"",40,0,3); h_c_mc_tag_zg[i]->Sumw2();
            h_b_mc_tag_zg[i] = new TH1D(Form("h_b_mc_tag_zg_%d",i),"",40,0,3); h_b_mc_tag_zg[i]->Sumw2();
    
            h_usdg_mc_untag_zg[i] = new TH1D(Form("h_usdg_mc_untag_zg_%d",i),"",40,0,3); h_usdg_mc_untag_zg[i]->Sumw2();
            h_c_mc_untag_zg[i] = new TH1D(Form("h_c_mc_untag_zg_%d",i),"",40,0,3); h_c_mc_untag_zg[i]->Sumw2();
            h_b_mc_untag_zg[i] = new TH1D(Form("h_b_mc_untag_zg_%d",i),"",40,0,3); h_b_mc_untag_zg[i]->Sumw2();

            h_usdg_mc_negTag_zg[i] = new TH1D(Form("h_usdg_mc_negTag_zg_%d",i),"",30,-10,3); h_usdg_mc_negTag_zg[i]->Sumw2();
            h_c_mc_negTag_zg[i] = new TH1D(Form("h_c_mc_negTag_zg_%d",i),"",30,-10,3); h_c_mc_negTag_zg[i]->Sumw2();
            h_b_mc_negTag_zg[i] = new TH1D(Form("h_b_mc_negTag_zg_%d",i),"",30,-10,3); h_b_mc_negTag_zg[i]->Sumw2();

            h_usdg_mc_negTagCut_zg[i] = new TH1D(Form("h_usdg_mc_negTagCut_zg_%d",i),"",30,-10,3); h_usdg_mc_negTagCut_zg[i]->Sumw2();
            h_c_mc_negTagCut_zg[i] = new TH1D(Form("h_c_mc_negTagCut_zg_%d",i),"",30,-10,3); h_c_mc_negTagCut_zg[i]->Sumw2();
            h_b_mc_negTagCut_zg[i] = new TH1D(Form("h_b_mc_negTagCut_zg_%d",i),"",30,-10,3); h_b_mc_negTagCut_zg[i]->Sumw2();
    
            h_data_negTag_zg[i] = new TH1D(Form("h_data_negTag_zg_%d",i),"",30,-10,3); h_data_negTag_zg[i]->Sumw2();
            h_data_negTagCut_zg[i] = new TH1D(Form("h_data_negTagCut_zg_%d",i),"",30,-10,3); h_data_negTagCut_zg[i]->Sumw2();
        }
	
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
	
	double jtpt, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, minJP, discr_csvV2, ndiscr_csvV2, discr_prob, svtxm, muPt;
	double jtSubJetJP1, jtSubJetJP2, jtSubJetCSV1, jtSubJetCSV2, jtSubJetNegCSV1, jtSubJetNegCSV2, jtSubJetSvtxm1, jtSubJetSvtxm2;
	int hiBin;
	double weight;
	
	//start with data
	//tdata->SetBranchAddress("jtpt",&jtpt);
        tdata->SetBranchAddress("matchedJtpt",&jtpt);
	tdata->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tdata->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tdata->SetBranchAddress("minCSV",&minCSV);
	tdata->SetBranchAddress("discr_prob", &discr_prob);
        tdata->SetBranchAddress("svtxm",&svtxm);
        tdata->SetBranchAddress("minJP",&minJP);
        tdata->SetBranchAddress("muPt",&muPt);
	tdata->SetBranchAddress("ndiscr_csvV2",&ndiscr_csvV2);
	tdata->SetBranchAddress("discr_csvV2",&discr_csvV2);
	tdata->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	tdata->SetBranchAddress("jtSubJetJP_1",&jtSubJetJP1);
	tdata->SetBranchAddress("jtSubJetJP_2",&jtSubJetJP2);
	tdata->SetBranchAddress("jtSubJetcsvV2_1",&jtSubJetCSV1);
	tdata->SetBranchAddress("jtSubJetcsvV2_2",&jtSubJetCSV2);
	tdata->SetBranchAddress("jtSubJetNegCsvV2_1",&jtSubJetNegCSV1);
	tdata->SetBranchAddress("jtSubJetNegCsvV2_2",&jtSubJetNegCSV2);
        tdata->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
        tdata->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);

	tdata->SetBranchAddress("weight",&weight);
	
	if(!ispp) tdata->SetBranchAddress("hiBin",&hiBin);
        for(int ievt=0; ievt<tdata->GetEntries(); ievt++){
		tdata->GetEntry(ievt);
		if((!doSubjets && jtpt>jtptCutMin) || (doSubjets && jtpt>jtptCutMin)){

                        if(jtpt>jtptCutMax) continue;
                        if(jtSubJetdR<0.1) continue;
                        if(muPt<4) continue;

                        //weight*=1e11;
			int ptBin=0;
			while(jtpt>xbins[ptBin+1] && ptBin<ptBins-1) ptBin++;
			
			int ptBin1=0, ptBin2=0;
			while(jtSubJetPt1>xbins[ptBin1+1] && ptBin1<ptBins-1) ptBin1++;
			while(jtSubJetPt2>xbins[ptBin2+1] && ptBin2<ptBins-1) ptBin2++;

                        int zgBin=0;
                        double zg;
                        if(!switchToDR) zg = getZg(jtSubJetPt1,jtSubJetPt2);
                        else zg = jtSubJetdR;
                        while(zg>xzgBins[zgBin+1]) zgBin++;

                        double minNegCSV = TMath::Min(jtSubJetNegCSV1,jtSubJetNegCSV2);

			if(!doSubjets){
				h_data[ptBin]->Fill(discr_prob, weight);
				if(discr_csvV2>taggerVal && svtxm>svtxCut) h_data_tag[ptBin]->Fill(discr_prob, weight);
				else h_data_untag[ptBin]->Fill(discr_prob, weight);
				
				if(ndiscr_csvV2>0) h_data_negTag[ptBin]->Fill(discr_prob, weight);
				if(ndiscr_csvV2>taggerVal && svtxm>svtxCut) h_data_negTagCut[ptBin]->Fill(discr_prob, weight); 
			}
			else{
                            if(jtSubJetdR<0.1) continue;
				h_data[ptBin1]->Fill(jtSubJetJP1, weight);
				h_data[ptBin2]->Fill(jtSubJetJP2, weight);
				if(jtSubJetCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_data_tag[ptBin1]->Fill(jtSubJetJP1, weight);
				else h_data_untag[ptBin1]->Fill(jtSubJetJP1, weight);
				if(jtSubJetCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_data_tag[ptBin2]->Fill(jtSubJetJP2, weight);
				else h_data_untag[ptBin2]->Fill(jtSubJetJP2, weight);
				
				if(jtSubJetNegCSV1>0 ) h_data_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
				if(jtSubJetNegCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_data_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);
				if(jtSubJetNegCSV2>0 ) h_data_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
				if(jtSubJetNegCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_data_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);
    
                                h_data_zg[zgBin]->Fill(minJP, weight);     
                                if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut ) h_data_tag_zg[zgBin]->Fill(minJP, weight);
                                else h_data_untag_zg[zgBin]->Fill(minJP, weight);

                                if(minNegCSV>0 ) h_data_negTag_zg[zgBin]->Fill(minJP, weight);
                                if(minNegCSV>taggerVal & jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_data_negTagCut_zg[zgBin]->Fill(minJP, weight);

			}
		}	
	}
	
	//then do MC
	int jtHadronFlavor, jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
	//tmc->SetBranchAddress("jtpt",&jtpt);
	tmc->SetBranchAddress("matchedJtpt",&jtpt);
        tmc->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tmc->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tmc->SetBranchAddress("minCSV",&minCSV);
	tmc->SetBranchAddress("minJP",&minJP);
	tmc->SetBranchAddress("svtxm",&svtxm);
        tmc->SetBranchAddress("muPt",&muPt);
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
        if(!ispp) tmc->SetBranchAddress("hiBin",&hiBin);
	for(int ievt=0; ievt<tmc->GetEntries(); ievt++){
		tmc->GetEntry(ievt);
		if((!doSubjets && jtpt>jtptCutMin) || (doSubjets && jtpt>jtptCutMin)){

                        if(jtpt>jtptCutMax) continue;
                        if(jtSubJetdR<0.1) continue;
                        if(muPt<5) continue;

                        //if(weight>1e-10) continue; //kill off pthat50 samples
			
			int ptBin=0;
			while(jtpt>xbins[ptBin+1] && ptBin<ptBins-1) ptBin++;

			int ptBin1=0, ptBin2=0;
			while(jtSubJetPt1>xbins[ptBin1+1] && ptBin1<ptBins-1) ptBin1++;
			while(jtSubJetPt2>xbins[ptBin2+1] && ptBin2<ptBins-1) ptBin2++;
    
                        int zgBin=0;
                        double zg=0;
                        if(!switchToDR) zg = getZg(jtSubJetPt1,jtSubJetPt2);
                        else zg = jtSubJetdR;
                        while(zg>xzgBins[zgBin+1]) zgBin++;

                        double minNegCSV = TMath::Min(jtSubJetNegCSV1,jtSubJetNegCSV2);

			if(!doSubjets){
                            if(abs(jtHadronFlavor)==4){
                                h_c_mc[ptBin]->Fill(discr_prob, weight);
                                if(discr_csvV2>taggerVal) h_c_mc_tag[ptBin]->Fill(discr_prob, weight);
                                else h_c_mc_untag[ptBin]->Fill(discr_prob, weight);

                                if(ndiscr_csvV2>0) h_c_mc_negTag[ptBin]->Fill(discr_prob, weight);
                                if(ndiscr_csvV2>taggerVal) h_c_mc_negTagCut[ptBin]->Fill(discr_prob, weight);
                            }
                            else if(abs(jtHadronFlavor)!=5){
					h_usdg_mc[ptBin]->Fill(discr_prob, weight);
					if(discr_csvV2>taggerVal && svtxm>svtxCut) h_usdg_mc_tag[ptBin]->Fill(discr_prob, weight);
					else h_usdg_mc_untag[ptBin]->Fill(discr_prob, weight);
					
					if(ndiscr_csvV2>0) h_usdg_mc_negTag[ptBin]->Fill(discr_prob, weight);
					if(ndiscr_csvV2>taggerVal && svtxm>svtxCut) h_usdg_mc_negTagCut[ptBin]->Fill(discr_prob, weight);
                                }
			}
			else{
                            if(jtSubJetdR<0.1) continue;
                            
                            if(abs(jtSubJetHadronFlavor1)==4){
                                h_c_mc[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetCSV1>taggerVal) h_c_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
                                //else h_c_mc_untag[ptBin1]->Fill(jtSubJetJP1, weight);

                                if(jtSubJetNegCSV1>0) h_c_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetNegCSV1>taggerVal) h_c_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);
                           
                                //CC + CB pairs
                                if(abs(jtSubJetHadronFlavor2)==4 || abs(jtSubJetHadronFlavor2)==5){
                                    h_c_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_c_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    //else h_c_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                    if(minNegCSV>0 ) h_c_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_c_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }
                                
                            }
                            /*else if(abs(jtSubJetHadronFlavor1)!=5){
                                        h_usdg_mc[ptBin1]->Fill(jtSubJetJP1, weight);
					if(jtSubJetCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_usdg_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
					else h_usdg_mc_untag[ptBin1]->Fill(jtSubJetJP1, weight);
					
					if(jtSubJetNegCSV1>0 ) h_usdg_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
					if(jtSubJetNegCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_usdg_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);
                               
                                        //LL + LC + LB pairs
                                        h_usdg_mc_zg[zgBin]->Fill(minJP, weight);
                                        if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                        else h_usdg_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                        if(minCSV>0 ) h_usdg_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                        if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);

                            }*/
                            if(abs(jtSubJetHadronFlavor2)==4){
                                h_c_mc[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_c_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
                                //else h_c_mc_untag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>0) h_c_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_c_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);
                            
                                // BC pairs
                                if(abs(jtSubJetHadronFlavor1)==5){
                                    h_c_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_c_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    //else h_c_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                    if(minNegCSV>0 ) h_c_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_c_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }
                            
                            }
                            /*else if(abs(jtSubJetHadronFlavor2)!=5){
                                h_usdg_mc[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_usdg_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
                                else h_usdg_mc_untag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>0 ) h_usdg_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_usdg_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);

                                //BL + CL pairs
                                if(abs(jtSubJetHadronFlavor1)==5 || abs(jtSubJetHadronFlavor1)==4){
                                    h_usdg_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    else h_usdg_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>0 ) h_usdg_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }

                            }*/
                            if(abs(jtSubJetHadronFlavor1)==5){
                                h_b_mc[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_b_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
                                //else h_b_mc_untag[ptBin1]->Fill(jtSubJetJP1, weight);

                                if(jtSubJetNegCSV1>0 ) h_b_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetNegCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_b_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);

                                //BB pairs
                                if(abs(jtSubJetHadronFlavor2)==5){
                                    h_b_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_b_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    //else h_b_mc_untag_zg[zgBin]->Fill(minJP, weight);

                                    if(minNegCSV>0 ) h_b_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_b_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }
                            }
                            if(abs(jtSubJetHadronFlavor2)==5){
                                h_b_mc[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_b_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
                                //else h_b_mc_untag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>0 ) h_b_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_b_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);
                            }

                        }	
		}
	}
	
	
	//then do lightMC
	tlightmc->SetBranchAddress("jtpt",&jtpt);
	tlightmc->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
	tlightmc->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
	tlightmc->SetBranchAddress("minCSV",&minCSV);
	tlightmc->SetBranchAddress("minJP",&minJP);
	tlightmc->SetBranchAddress("svtxm",&svtxm);
        tlightmc->SetBranchAddress("muPt",&muPt);
        tlightmc->SetBranchAddress("discr_prob", &discr_prob);
	tlightmc->SetBranchAddress("ndiscr_csvV2",&ndiscr_csvV2);
	tlightmc->SetBranchAddress("discr_csvV2",&discr_csvV2);
	tlightmc->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
	tlightmc->SetBranchAddress("jtHadronFlavor",&jtHadronFlavor);
	tlightmc->SetBranchAddress("weight",&weight);
	tlightmc->SetBranchAddress("jtSubJetJP_1",&jtSubJetJP1);
	tlightmc->SetBranchAddress("jtSubJetJP_2",&jtSubJetJP2);
	tlightmc->SetBranchAddress("jtSubJetcsvV2_1",&jtSubJetCSV1);
	tlightmc->SetBranchAddress("jtSubJetcsvV2_2",&jtSubJetCSV2);
	tlightmc->SetBranchAddress("jtSubJetNegCsvV2_1",&jtSubJetNegCSV1);
	tlightmc->SetBranchAddress("jtSubJetNegCsvV2_2",&jtSubJetNegCSV2);
	tlightmc->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
	tlightmc->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
        tlightmc->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
        tlightmc->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
        if(!ispp) tlightmc->SetBranchAddress("hiBin",&hiBin);
	for(int ievt=0; ievt<tlightmc->GetEntries(); ievt++){
		tlightmc->GetEntry(ievt);
		if((!doSubjets && jtpt>jtptCutMin) || (doSubjets && jtpt>jtptCutMin)){

                    if(jtpt>jtptCutMax) continue;
                    if(jtSubJetdR<0.1) continue;
                    if(muPt<5) continue;

			int ptBin=0;
			while(jtpt>xbins[ptBin+1] && ptBin<ptBins-1) ptBin++;

			int ptBin1=0, ptBin2=0;
			while(jtSubJetPt1>xbins[ptBin1+1] && ptBin1<ptBins-1) ptBin1++;
			while(jtSubJetPt2>xbins[ptBin2+1] && ptBin2<ptBins-1) ptBin2++;

                        int zgBin=0;
                        double zg=0;
                        if(!switchToDR) zg = getZg(jtSubJetPt1,jtSubJetPt2);
                        else zg = jtSubJetdR;
                        while(zg>xzgBins[zgBin+1]) zgBin++;

                        if(doSubjets){
                            if(jtSubJetdR<0.1) continue;
                        }
                        double minNegCSV = TMath::Min(jtSubJetNegCSV1,jtSubJetNegCSV2);

			if(!doSubjets){
				if(abs(jtHadronFlavor)!=5 && abs(jtHadronFlavor)!=4){
					h_usdg_mc[ptBin]->Fill(discr_prob, weight);
					if(discr_csvV2>taggerVal && svtxm>svtxCut) h_usdg_mc_tag[ptBin]->Fill(discr_prob, weight);
					else h_usdg_mc_untag[ptBin]->Fill(discr_prob, weight);
					
					if(ndiscr_csvV2>0) h_usdg_mc_negTag[ptBin]->Fill(discr_prob, weight);
					if(ndiscr_csvV2>taggerVal && svtxm>svtxCut) h_usdg_mc_negTagCut[ptBin]->Fill(discr_prob, weight);
				}
			}
                        else{
                            if(abs(jtSubJetHadronFlavor1)==4){
                                //h_c_mc[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetCSV1>taggerVal){}// h_c_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
                                else h_c_mc_untag[ptBin1]->Fill(jtSubJetJP1, weight);

                                //if(jtSubJetNegCSV1>0) h_c_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
                                //if(jtSubJetNegCSV1>taggerVal) h_c_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);

                                //CC + CB pairs
                                if(abs(jtSubJetHadronFlavor2)==4 || abs(jtSubJetHadronFlavor2)==5){
                                    //h_c_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut){}// h_c_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    else h_c_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                    //if(minNegCSV>0 ) h_c_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    //if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_c_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }

                            }
                            if(abs(jtSubJetHadronFlavor2)==4){
                                //h_c_mc[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetCSV2>taggerVal && jtSubJetSvtxm2>svtxCut){}// h_c_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
                                else h_c_mc_untag[ptBin2]->Fill(jtSubJetJP2, weight);
                                //if(jtSubJetNegCSV2>0) h_c_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
                                //if(jtSubJetNegCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_c_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);

                                // BC pairs
                                if(abs(jtSubJetHadronFlavor1)==5){
                                    //h_c_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut){}// h_c_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    else h_c_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                    //if(minNegCSV>0 ) h_c_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    //if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_c_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }

                            }
                            if(abs(jtSubJetHadronFlavor1)==5){
                                //h_b_mc[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetCSV1>taggerVal && jtSubJetSvtxm1>svtxCut){}//h_b_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
                                else h_b_mc_untag[ptBin1]->Fill(jtSubJetJP1, weight);

                                //if(jtSubJetNegCSV1>0 ) h_b_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
                                //if(jtSubJetNegCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_b_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);

                                //BB pairs
                                if(abs(jtSubJetHadronFlavor2)==5){
                                    //h_b_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut){}// h_b_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    else h_b_mc_untag_zg[zgBin]->Fill(minJP, weight);

                                    //if(minNegCSV>0 ) h_b_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    //if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_b_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }
                            }
                            if(abs(jtSubJetHadronFlavor2)==5){
                                //h_b_mc[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetCSV2>taggerVal && jtSubJetSvtxm2>svtxCut){}// h_b_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
                                else h_b_mc_untag[ptBin2]->Fill(jtSubJetJP2, weight);
                                //if(jtSubJetNegCSV2>0 ) h_b_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
                                //if(jtSubJetNegCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_b_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);
                            }

                            //original code starts here
                            if(abs(jtSubJetHadronFlavor1)!=4 && abs(jtSubJetHadronFlavor1)!=5){
                                h_usdg_mc[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_usdg_mc_tag[ptBin1]->Fill(jtSubJetJP1, weight);
                                else h_usdg_mc_untag[ptBin1]->Fill(jtSubJetJP1, weight);

                                if(jtSubJetNegCSV1>0 ) h_usdg_mc_negTag[ptBin1]->Fill(jtSubJetJP1, weight);
                                if(jtSubJetNegCSV1>taggerVal && jtSubJetSvtxm1>svtxCut) h_usdg_mc_negTagCut[ptBin1]->Fill(jtSubJetJP1, weight);

                                //LL + LC + LB pairs
                                h_usdg_mc_zg[zgBin]->Fill(minJP, weight);
                                if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                else h_usdg_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                if(minNegCSV>0 ) h_usdg_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);

                            }
                            if(abs(jtSubJetHadronFlavor2)!=4 && abs(jtSubJetHadronFlavor2)!=5){
                                h_usdg_mc[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_usdg_mc_tag[ptBin2]->Fill(jtSubJetJP2, weight);
                                else h_usdg_mc_untag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>0 ) h_usdg_mc_negTag[ptBin2]->Fill(jtSubJetJP2, weight);
                                if(jtSubJetNegCSV2>taggerVal && jtSubJetSvtxm2>svtxCut) h_usdg_mc_negTagCut[ptBin2]->Fill(jtSubJetJP2, weight);

                                //BL + CL pairs
                                if(abs(jtSubJetHadronFlavor1)==5 || abs(jtSubJetHadronFlavor1)==4){
                                    h_usdg_mc_zg[zgBin]->Fill(minJP, weight);
                                    if(minCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_tag_zg[zgBin]->Fill(minJP, weight);
                                    else h_usdg_mc_untag_zg[zgBin]->Fill(minJP, weight);
                                    if(minNegCSV>0 ) h_usdg_mc_negTag_zg[zgBin]->Fill(minJP, weight);
                                    if(minNegCSV>taggerVal && jtSubJetSvtxm1>svtxCut && jtSubJetSvtxm2>svtxCut) h_usdg_mc_negTagCut_zg[zgBin]->Fill(minJP, weight);
                                }
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
        for(int i=0; i<zgBins; i++){
                h_data_zg[i]->Write();
                h_data_tag_zg[i]->Write();
                h_data_untag_zg[i]->Write();

                h_usdg_mc_zg[i]->Write();
                h_c_mc_zg[i]->Write();
                h_b_mc_zg[i]->Write();

                h_usdg_mc_tag_zg[i]->Write();
                h_c_mc_tag_zg[i]->Write();
                h_b_mc_tag_zg[i]->Write();

                h_usdg_mc_untag_zg[i]->Write();
                h_c_mc_untag_zg[i]->Write();
                h_b_mc_untag_zg[i]->Write();

                h_usdg_mc_negTag_zg[i]->Write();
                h_c_mc_negTag_zg[i]->Write();
                h_b_mc_negTag_zg[i]->Write();

                h_usdg_mc_negTagCut_zg[i]->Write();
                h_c_mc_negTagCut_zg[i]->Write();
                h_b_mc_negTagCut_zg[i]->Write();

                h_data_negTag_zg[i]->Write();
                h_data_negTagCut_zg[i]->Write();
        }
	out->Close();
}
