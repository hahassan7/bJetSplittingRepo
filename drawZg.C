
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

void correctSubjetReco(TF1 *corrFactor, TH1D *input, bool doDivide=0){
	for(int ibin=1; ibin<=input->GetNbinsX(); ibin++){
		double binCenter = input->GetBinCenter(ibin);
		if(!doDivide){
			input->SetBinContent(ibin, input->GetBinContent(ibin)*corrFactor->Eval(binCenter));
			input->SetBinError(ibin, input->GetBinError(ibin)*corrFactor->Eval(binCenter));
		}
		else{
			if(corrFactor->Eval(binCenter)<0) continue;
			//cout << "for " << binCenter << " corr factor is " << corrFactor->Eval(binCenter) << endl;
			
			double offset = 0.;//abs(corrFactor->Eval(0.))+0.94;
			input->SetBinContent(ibin, input->GetBinContent(ibin)/(corrFactor->Eval(binCenter)+offset));
			input->SetBinError(ibin, input->GetBinError(ibin)/(corrFactor->Eval(binCenter)+offset));
			
		}
	}
	
}

double getYMean(TH1D *hist, TH1D *vals){
	double runTotal = 0;
	for(int ibin=1; ibin<=hist->GetNbinsX(); ibin++){
		runTotal += hist->GetBinContent(ibin)*vals->GetBinContent(ibin);
	}
	return runTotal/vals->Integral();
}

void drawZg(double cut=0.8, double jpCut=-9999, double svtxmCut=1.2){
	
	bool drawFromTree = true;
	bool writeOut = true;
	bool ispp = true;
	
	bool fitSFs = true;
	
	if(jpCut>0 && svtxmCut>0){
		cout << "WARNING! Both JP and SVTXM cuts are >0! You don't have enough stats for this!" << endl;
		exit(0);
	}
	
	TFile *fppdata = new TFile("ntuples/gspTreeOut_ppData_trgFilter_lowPtCuts.root");
	TFile *fppmc;
	if(!ispp) fppmc = new TFile("ntuples/gspTreeOut_PbPbMC_explicitJTA.root");
	//else fppmc = new TFile("ntuples/gspTreeOut_ppMC_muFilteredMC_Pythia8_fullSample_withFullJets.root");
	else fppmc = new TFile("ntuples/gspTreeOut_ppMC_QCDJet_noMuFilter_withRefInfo.root");
	//else fppmc = new TFile("ntuples/gspTreeOut_ppMC_muFilteredMC_Pythia8_fullSample.root");
	
	TFile *fppJES = new TFile("ntuples/gspTreeOut_ppMC_QCDJet_noMuFilter_drop3p9pctTracks.root");
	
	TFile *fpowheg = new TFile("ntuples/gspTreeOut_powhegSum.root");
	
	TH1D *mcLightFrac[2], *dataLightFrac[2], *mcPur[2], *dataPur[2], *mcEff[2], *dataEff[2], *scaleFactor[2];
	TH1D *jesEff[2];
	
	TFile *pursAndEffsFull = new TFile(Form("efficiencies/sf_fullJets_csv%g.root",cut));
	TFile *pursAndEffsSub[2];
	TFile *pursAndEffsJPCalibSys = new TFile(Form("efficiencies/sf_subjets_zg__muFiltered_csv%g_sv%g_dataJPCalibOnMC_ptBin1.root",cut, svtxmCut));
	TFile *pursAndEffsDMesonReweight = new TFile(Form("efficiencies/sf_subjets_zg__muFiltered_csv%g_sv%g_DMesonReweight_ptBin1.root",cut, svtxmCut));
	TFile *pursAndEffsJES[2]; 
	
	for(int ipt=0; ipt<2; ipt++){
		if(jpCut<-1000 && svtxmCut<-1000) pursAndEffsSub[ipt] = new TFile(Form("efficiencies/sf_subjets__muFiltered_csv%g_fullJetPtBin%d.root",cut,ipt+1));
		else if(jpCut==1.2) pursAndEffsSub[ipt] = new TFile(Form("efficiencies/sf_subjets_zg_csv%g_JP1p2Sel_fullJetPtBin%d.root",cut,ipt+1));
		else if(svtxmCut>0) pursAndEffsSub[ipt] = new TFile(Form("efficiencies/sf_subjets_zg__muFiltered_csv%g_sv%g_fullJetPtBin%d.root",cut,svtxmCut,ipt+1));
		else{ cout << "Can't find file with cut=" << cut << " svtxmCut=" << svtxmCut << endl;}
		
		mcLightFrac[ipt] = (TH1D*)pursAndEffsSub[ipt]->Get("mcLightFrac");
		dataLightFrac[ipt] = (TH1D*)pursAndEffsSub[ipt]->Get("dataLightFrac");
		mcPur[ipt] = (TH1D*)pursAndEffsSub[ipt]->Get("mcPur");
		dataPur[ipt] = (TH1D*)pursAndEffsSub[ipt]->Get("dataPur");
		mcEff[ipt] = (TH1D*)pursAndEffsSub[ipt]->Get("mcEff");
		dataEff[ipt] = (TH1D*)pursAndEffsSub[ipt]->Get("dataEff");
		scaleFactor[ipt] = (TH1D*)pursAndEffsSub[ipt]->Get("scaleFactors_btagSF");
		
		pursAndEffsJES[ipt] = new TFile(Form("efficiencies/sf_subjets_zg__muFiltered_csv%g_sv1.2_fullJetPtBin%d_shiftDown.root",cut,ipt+1));
		jesEff[ipt] = (TH1D*)pursAndEffsJES[ipt]->Get("dataEff");
	}
	
	TH1D *mcEffJPSys = (TH1D*)pursAndEffsJPCalibSys->Get("mcEff");
	TH1D *mcEffDMeson = (TH1D*)pursAndEffsDMesonReweight->Get("mcEff");
	
	const int ptBins = 2;
	const int jetPtBin[ptBins+1] = {120, 160, 400};
	const int jetPtBinShift[ptBins+1] = {116, 155, 400};
	//const int jetPtBinShift[ptBins+1] = {151, 171, 9999};
	TH1D *zgSpecIncl_data[ptBins];
	TH1D *zgSpecIncl_data_shift[ptBins];
	TH1D *zgSpecDR_data[ptBins];
	TH1D *zgSpecTag_data[ptBins];
	TH1D *zgSpecTag_data_shift[ptBins];
	TH1D *zgSpecBB_data[ptBins];
	TH1D *zgSpecBB_data_shift[ptBins];
	TH1D *zgSpecBB_data_v2[ptBins];
	
	TH1D *zgSpecIncl_mc_JES[ptBins];
	
	TH1D *drSpecIncl_data[ptBins];
	TH1D *drSpecGSP_data[ptBins];
	TH1D *drSpecIncl_data_shift[ptBins];
	TH1D *drSpecGSP_data_shift[ptBins];
	
	TH1D *zgSpecIncl_mc[ptBins];
	TH1D *zgSpecIncl_refAxis_mc[ptBins];
	TH1D *zgSpecTagCorr_mc[ptBins];
	TH1D *zgSpecTag_mc[ptBins];
	TH1D *zgSpecBB_mc[ptBins];
	TH1D *zgSpecBB_mc_JES[ptBins];
	TH1D *zgSpecBB_mc_jp0[ptBins];
	TH1D *zgSpecBB_mc_jp0p4[ptBins];
	TH1D *zgSpecBB_mc_jp1[ptBins];
	
	TH1D *drSpecIncl_mc[ptBins];
	TH1D *drSpecIncl_refAxis_mc[ptBins];
	TH1D *drSpecIncl_mc_JES[ptBins];
	TH1D *drSpecTag_mc[ptBins];
	TH1D *drSpecBB_mc[ptBins];
	
	TH1D *zgSpecPowheg[ptBins];
	TH1D *drSpecPowheg[ptBins];
	
	TH2D *zGdRCorr_mc[ptBins];
	TH2D *zGdRCorr_bb_mc[ptBins];
	TH2D *zGdRcorr_data[ptBins];
	TH2D *zGdRcorr_csv_data[ptBins];
	
	if(drawFromTree){
		for(int i=0; i<ptBins; i++){
			zgSpecIncl_data[i] = new TH1D(Form("zgSpecIncl_data_%d",i),"",10,0,0.5); zgSpecIncl_data[i]->Sumw2();
			zgSpecIncl_data_shift[i] = new TH1D(Form("zgSpecIncl_data_shift_%d",i),"",10,0,0.5); zgSpecIncl_data_shift[i]->Sumw2();
			zgSpecTag_data[i] = new TH1D(Form("zgSpecTag_data_%d",i),"",10,0,0.5); zgSpecTag_data[i]->Sumw2();
			zgSpecTag_data_shift[i] = new TH1D(Form("zgSpecTag_data_shift_%d",i),"",10,0,0.5); zgSpecTag_data_shift[i]->Sumw2();
			zgSpecBB_data[i] = new TH1D(Form("zgSpecBB_data_%d",i),"",10,0,0.5); zgSpecBB_data[i]->Sumw2();
			zgSpecBB_data_v2[i] = new TH1D(Form("zgSpecBB_data_v2_%d",i),"",10,0,0.5); zgSpecBB_data_v2[i]->Sumw2();
			zgSpecDR_data[i] = new TH1D(Form("zgSpecDR_data_%d",i),"",10,0,0.5); zgSpecDR_data[i]->Sumw2();
			
			drSpecIncl_data[i] = new TH1D(Form("drSpecIncl_data%d",i),"",20,0,0.4); drSpecIncl_data[i]->Sumw2();
			drSpecGSP_data[i] = new TH1D(Form("drSpecGSP_data%d",i),"",20,0,0.4); drSpecGSP_data[i]->Sumw2();
			drSpecIncl_data_shift[i] = new TH1D(Form("drSpecIncl_data_shift%d",i),"",20,0,0.4); drSpecIncl_data_shift[i]->Sumw2();
			drSpecGSP_data_shift[i] = new TH1D(Form("drSpecGSP_data_shift%d",i),"",20,0,0.4); drSpecGSP_data_shift[i]->Sumw2();
			
			zgSpecIncl_mc[i] = new TH1D(Form("zgSpecIncl_mc%d",i),"",10,0,0.5); zgSpecIncl_mc[i]->Sumw2();
			zgSpecIncl_mc_JES[i] = new TH1D(Form("zgSpecIncl_mc_JES%d",i),"",10,0,0.5); zgSpecIncl_mc_JES[i]->Sumw2();
			zgSpecIncl_refAxis_mc[i] = new TH1D(Form("zgSpecIncl_refAxis_mc%d",i),"",10,0,0.5); zgSpecIncl_refAxis_mc[i]->Sumw2();
			zgSpecTag_mc[i] = new TH1D(Form("zgSpecTag_mc%d",i),"",10,0,0.5); zgSpecTag_mc[i]->Sumw2();
			zgSpecTagCorr_mc[i] = new TH1D(Form("zgSpecCorr_mc%d",i),"",10,0,0.5); zgSpecTagCorr_mc[i]->Sumw2();
			zgSpecBB_mc[i] = new TH1D(Form("zgSpecBB_mc%d",i),"",10,0,0.5); zgSpecBB_mc[i]->Sumw2();
			zgSpecBB_mc_JES[i] = new TH1D(Form("zgSpecBB_mc_JES%d",i),"",10,0,0.5); zgSpecBB_mc_JES[i]->Sumw2();
			zgSpecBB_mc_jp0[i] = new TH1D(Form("zgSpecBB_mc_jp0_%d",i),"",10,0,0.5); zgSpecBB_mc_jp0[i]->Sumw2();
			zgSpecBB_mc_jp0p4[i] = new TH1D(Form("zgSpecBB_mc_jp0p4_%d",i),"",10,0,0.5); zgSpecBB_mc_jp0p4[i]->Sumw2();
			zgSpecBB_mc_jp1[i] = new TH1D(Form("zgSpecBB_mc_jp1_%d",i),"",10,0,0.5); zgSpecBB_mc_jp1[i]->Sumw2();
			
			drSpecIncl_mc[i] = new TH1D(Form("drSpecIncl_mc%d",i),"",20,0,0.4); drSpecIncl_mc[i]->Sumw2();
			drSpecIncl_mc_JES[i] = new TH1D(Form("drSpecIncl_mc_JES%d",i),"",20,0,0.4); drSpecIncl_mc_JES[i]->Sumw2();
			drSpecIncl_refAxis_mc[i] = new TH1D(Form("drSpecIncl_refAxis_mc%d",i),"",20,0,0.4); drSpecIncl_refAxis_mc[i]->Sumw2();
			drSpecTag_mc[i] = new TH1D(Form("drSpecTag_mc%d",i),"",20,0,0.4); drSpecTag_mc[i]->Sumw2();
			drSpecBB_mc[i] = new TH1D(Form("drSpecBB_mc%d",i),"",20,0,0.4); drSpecBB_mc[i]->Sumw2();
			
			zGdRCorr_mc[i] = new TH2D(Form("zGdRCorr_mc%d",i),"",10,0,0.5,40,0,0.4); zGdRCorr_mc[i]->Sumw2();
			zGdRCorr_bb_mc[i] = new TH2D(Form("zGdRCorr_bb_mc%d",i),"",10,0,0.5,40,0,0.4); zGdRCorr_bb_mc[i]->Sumw2();
			zGdRcorr_data[i] = new TH2D(Form("zGdRcorr_data%d",i),"",10,0,0.5,40,0,0.4); zGdRcorr_data[i]->Sumw2();
			zGdRcorr_csv_data[i] = new TH2D(Form("zGdRcorr_csv_data%d",i),"",10,0,0.5,40,0,0.4); zGdRcorr_csv_data[i]->Sumw2();
			
			zgSpecPowheg[i] = new TH1D(Form("zgSpecPowheg_%d",i),"",10,0,0.5); zgSpecPowheg[i]->Sumw2();
			drSpecPowheg[i] = new TH1D(Form("drSpecPowheg_%d",i),"",8,0,0.4); drSpecPowheg[i]->Sumw2();
			
		}
		
		//add some powheg
		TTree *powhegT = (TTree*)fpowheg->Get("gspTree");
		double jtpt, jteta, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, jtSubJetRefdR, minJP, jtSubJetSvtxm1, jtSubJetSvtxm2;
		int jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
		double jtSubJetHadronDR1, jtSubJetHadronDR2, jtSubJetHadronPt1, jtSubJetHadronPt2;
		double weight;
		powhegT->SetBranchAddress("matchedJtpt",&jtpt);
		powhegT->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
		powhegT->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
		powhegT->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
		powhegT->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
		powhegT->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
		powhegT->SetBranchAddress("weight",&weight);
		for(int ientry=0; ientry<powhegT->GetEntries(); ientry++){
			powhegT->GetEntry(ientry);
			
			if(jtpt<120) continue;
			if(jtpt>400) continue;
			
			int ibin=0;
			while(jtpt>jetPtBin[ibin+1]) ibin++;
			
			double zg = getzg(jtSubJetPt1, jtSubJetPt2);
			if(zg<0.1) continue;
			if(weight>10000) continue;
			//cout << "zg: "<< zg << endl;
			if(abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5 && jtSubJetdR>0.1){
				
				zgSpecPowheg[ibin]->Fill(zg, weight);
				drSpecPowheg[ibin]->Fill(jtSubJetdR, weight);
			}
		
		}
		for(int i=0; i<ptBins; i++){
			zgSpecPowheg[i]->Scale(1./zgSpecPowheg[i]->Integral());
			zgSpecPowheg[i]->Scale(1./zgSpecPowheg[i]->GetBinWidth(4));
			
			drSpecPowheg[i]->Scale(1./drSpecPowheg[i]->Integral());
			drSpecPowheg[i]->Scale(1./drSpecPowheg[i]->GetBinWidth(4));
		}
		
		
		//start with data
		TTree *dataT = (TTree*)fppdata->Get("gspTree");
		//double jtpt, jteta, jtSubJetPt1, jtSubJetPt2, minCSV, jtSubJetdR, minJP, jtSubJetSvtxm1, jtSubJetSvtxm2;
		int hiBin;
		dataT->SetBranchAddress("matchedJtpt",&jtpt);
		dataT->SetBranchAddress("jteta",&jteta);
		dataT->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
		dataT->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
		dataT->SetBranchAddress("minCSV",&minCSV);
		dataT->SetBranchAddress("minJP",&minJP);
		dataT->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
		dataT->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
		dataT->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
		if(!ispp) dataT->SetBranchAddress("hiBin",&hiBin);
		for(int ientry=0; ientry<dataT->GetEntries(); ientry++){
			dataT->GetEntry(ientry);
			
			if(jtpt<110) continue;
			if(jtpt>400) continue;
			//if(abs(jteta)<1.3) continue;
			//if(jtpt<160) continue;
			
			//if(!ispp){
			//	if(hiBin>20) continue;
			//}
			
			int ibin=0, ibinShift=0;
			while(jtpt>jetPtBin[ibin+1]) ibin++;
			
			while(jtpt>jetPtBinShift[ibinShift+1]) ibinShift++;
			
			double zg = getzg(jtSubJetPt1, jtSubJetPt2);
			if(zg<0.1) continue;
			if(jtpt>120 && jtSubJetdR>0.1) zgSpecIncl_data[ibin]->Fill(zg);
			if(jtpt>110 && jtSubJetdR>0.1) zgSpecIncl_data_shift[ibinShift]->Fill(zg);
			if(jtpt>120 && minCSV>cut && minJP>jpCut && jtSubJetSvtxm1>svtxmCut && jtSubJetSvtxm2>svtxmCut && jtSubJetdR>0.1) zgSpecTag_data[ibin]->Fill(zg);
			if(jtpt>110 && minCSV>cut && minJP>jpCut && jtSubJetSvtxm1>svtxmCut && jtSubJetSvtxm2>svtxmCut && jtSubJetdR>0.1) zgSpecTag_data_shift[ibinShift]->Fill(zg);
			if(jtSubJetdR>0.1) zgSpecDR_data[ibin]->Fill(zg);
			//if(jtSubJetdR>0.2 && minCSV>cut) zgSpecGSP_data[ibin]->Fill(zg);
			
			if(jtpt>120 && jtSubJetdR>0.1) zGdRcorr_data[ibin]->Fill(zg, jtSubJetdR);
			if(jtpt>120 && minCSV>cut && minJP>jpCut && jtSubJetSvtxm1>svtxmCut && jtSubJetSvtxm2>svtxmCut && jtSubJetdR>0.1) zGdRcorr_csv_data[ibin]->Fill(zg, jtSubJetdR);
			
			if(jtpt>120 && jtSubJetdR>0.1) drSpecIncl_data[ibin]->Fill(jtSubJetdR);
			if(jtpt>110 && jtSubJetdR>0.1) drSpecIncl_data_shift[ibinShift]->Fill(jtSubJetdR);
			if(jtpt>120 && minCSV>cut && minJP>jpCut && jtSubJetSvtxm1>svtxmCut && jtSubJetSvtxm2>svtxmCut && jtSubJetdR>0.1) drSpecGSP_data[ibin]->Fill(jtSubJetdR);
			if(jtpt>110 && minCSV>cut && minJP>jpCut && jtSubJetSvtxm1>svtxmCut && jtSubJetSvtxm2>svtxmCut && jtSubJetdR>0.1) drSpecGSP_data_shift[ibinShift]->Fill(jtSubJetdR);
		}
		for(int i=0; i<ptBins; i++){
			zgSpecIncl_data[i]->Scale(1./zgSpecIncl_data[i]->Integral());
			zgSpecIncl_data_shift[i]->Scale(1./zgSpecIncl_data_shift[i]->Integral());
			zgSpecDR_data[i]->Scale(1./zgSpecDR_data[i]->Integral());
			//zgSpecGSP_data[i]->Scale(1./zgSpecGSP_data[i]->Integral());
			zgSpecTag_data[i]->Scale(1./zgSpecTag_data[i]->Integral());
			
			zgSpecIncl_data[i]->Scale(1./zgSpecIncl_data[i]->GetBinWidth(4));
			zgSpecIncl_data_shift[i]->Scale(1./zgSpecIncl_data_shift[i]->GetBinWidth(4));
			zgSpecDR_data[i]->Scale(1./zgSpecDR_data[i]->GetBinWidth(4));
			//zgSpecGSP_data[i]->Scale(1./zgSpecGSP_data[i]->GetBinWidth(4));
			zgSpecTag_data[i]->Scale(1./zgSpecTag_data[i]->GetBinWidth(4));
			
			drSpecIncl_data[i]->Scale(1./drSpecIncl_data[i]->Integral());
			drSpecGSP_data[i]->Scale(1./drSpecGSP_data[i]->Integral());
			
			drSpecIncl_data[i]->Scale(1./drSpecIncl_data[i]->GetBinWidth(2));
			drSpecGSP_data[i]->Scale(1./drSpecGSP_data[i]->GetBinWidth(2));
			
			drSpecIncl_data_shift[i]->Scale(1./drSpecIncl_data_shift[i]->Integral());
			drSpecGSP_data_shift[i]->Scale(1./drSpecGSP_data_shift[i]->Integral());
			
			drSpecIncl_data_shift[i]->Scale(1./drSpecIncl_data_shift[i]->GetBinWidth(2));
			drSpecGSP_data_shift[i]->Scale(1./drSpecGSP_data_shift[i]->GetBinWidth(2));
			
		}
		
		TTree *mcjesT = (TTree*)fppJES->Get("gspTree");
		//int jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
		//double jtSubJetHadronDR1, jtSubJetHadronDR2, jtSubJetHadronPt1, jtSubJetHadronPt2;
		//double weight;
		mcjesT->SetBranchAddress("matchedJtpt",&jtpt);
		mcjesT->SetBranchAddress("jteta",&jteta);
		mcjesT->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
		mcjesT->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
		mcjesT->SetBranchAddress("minCSV",&minCSV);
		mcjesT->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
		mcjesT->SetBranchAddress("jtSubJetRefdR",&jtSubJetRefdR);
		mcjesT->SetBranchAddress("jtSubJetHadronFlavor1",&jtSubJetHadronFlavor1);
		mcjesT->SetBranchAddress("jtSubJetHadronFlavor2",&jtSubJetHadronFlavor2);
		mcjesT->SetBranchAddress("jtSubJetHadronDR1",&jtSubJetHadronDR1);
		mcjesT->SetBranchAddress("jtSubJetHadronDR2",&jtSubJetHadronDR2);
		mcjesT->SetBranchAddress("jtSubJetHadronPt1",&jtSubJetHadronPt1);
		mcjesT->SetBranchAddress("jtSubJetHadronPt2",&jtSubJetHadronPt2);
		mcjesT->SetBranchAddress("minJP",&minJP);
		mcjesT->SetBranchAddress("jtSubJetSvtxm1",&jtSubJetSvtxm1);
		mcjesT->SetBranchAddress("jtSubJetSvtxm2",&jtSubJetSvtxm2);
		mcjesT->SetBranchAddress("weight",&weight);
		if(!ispp) mcjesT->SetBranchAddress("hiBin",&hiBin);
		for(int ientry=0; ientry<mcjesT->GetEntries(); ientry++){
			mcjesT->GetEntry(ientry);
			
			if(jtpt<120) continue;
			if(jtpt>400) continue;
			
			int ibin=0;
			while(jtpt>jetPtBin[ibin+1]) ibin++;
			
			double zg = getzg(jtSubJetPt1, jtSubJetPt2);
			if(zg<0.1) continue;
			if(jtSubJetdR>0.1) zgSpecIncl_mc_JES[ibin]->Fill(zg, weight);
			if(jtSubJetdR>0.1) drSpecIncl_mc_JES[ibin]->Fill(jtSubJetdR, weight);
			
			if(abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5 && jtSubJetdR>0.1){
				zgSpecBB_mc_JES[ibin]->Fill(zg, weight);
			}
			
		}
		
		for(int i=0; i<ptBins; i++){
			zgSpecIncl_mc_JES[i]->Scale(1./zgSpecIncl_mc_JES[i]->Integral());
			zgSpecIncl_mc_JES[i]->Scale(1./zgSpecIncl_mc_JES[i]->GetBinWidth(4));
			
			drSpecIncl_mc_JES[i]->Scale(1./drSpecIncl_mc_JES[i]->Integral());
			drSpecIncl_mc_JES[i]->Scale(1./drSpecIncl_mc_JES[i]->GetBinWidth(4));
			
			zgSpecBB_mc_JES[i]->Scale(1./zgSpecBB_mc_JES[i]->Integral());
			zgSpecBB_mc_JES[i]->Scale(1./zgSpecBB_mc_JES[i]->GetBinWidth(4));
		}
		
		//finish with MC
		TTree *mcT = (TTree*)fppmc->Get("gspTree");
		//int jtSubJetHadronFlavor1, jtSubJetHadronFlavor2;
		//double jtSubJetHadronDR1, jtSubJetHadronDR2, jtSubJetHadronPt1, jtSubJetHadronPt2;
		//double weight;
		mcT->SetBranchAddress("matchedJtpt",&jtpt);
		mcT->SetBranchAddress("jteta",&jteta);
		mcT->SetBranchAddress("jtSubJetPt1",&jtSubJetPt1);
		mcT->SetBranchAddress("jtSubJetPt2",&jtSubJetPt2);
		mcT->SetBranchAddress("minCSV",&minCSV);
		mcT->SetBranchAddress("jtSubJetdR",&jtSubJetdR);
		mcT->SetBranchAddress("jtSubJetRefdR",&jtSubJetRefdR);
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
		if(!ispp) mcT->SetBranchAddress("hiBin",&hiBin);
		for(int ientry=0; ientry<mcT->GetEntries(); ientry++){
			mcT->GetEntry(ientry);
			
			if(jtpt<120) continue;
			if(jtpt>400) continue;
			//if(abs(jteta)>1.3) continue;
			
			if(!ispp){
				if(hiBin>20) continue;
			}
			
			//if(jtSubJetPt1<50 || jtSubJetPt2<10 || jtSubJetHadronDR2>0.1 || jtSubJetHadronDR2>0.1) continue;
			
			int ibin=0;
			while(jtpt>jetPtBin[ibin+1]) ibin++;
			
			double zg = getzg(jtSubJetPt1, jtSubJetPt2);
			if(zg<0.1) continue;
			if(jtSubJetdR>0.1) zgSpecIncl_mc[ibin]->Fill(zg, weight);
			if(abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5 && jtSubJetdR>0.1){
				
				zgSpecBB_mc[ibin]->Fill(zg, weight);
				if(jtSubJetHadronDR1<0.1 && jtSubJetHadronDR2<0.1) zgSpecBB_mc_jp0[ibin]->Fill(zg, weight);
				if(jtSubJetHadronDR1<0.1 && jtSubJetHadronDR2<0.1 && jtSubJetdR>0.1) zgSpecBB_mc_jp0p4[ibin]->Fill(zg, weight);
				if(jtSubJetHadronDR1<0.02 && jtSubJetHadronDR2<0.02) zgSpecBB_mc_jp1[ibin]->Fill(zg, weight);
			}
			
			zGdRCorr_mc[ibin]->Fill(zg, jtSubJetdR, weight);
			if(abs(jtSubJetHadronFlavor1)==5 && abs(jtSubJetHadronFlavor2)==5 && jtSubJetdR>0.1){
				zGdRCorr_bb_mc[ibin]->Fill(zg, jtSubJetdR, weight);
				drSpecBB_mc[ibin]->Fill(jtSubJetdR, weight);	
			}
			
			if(minCSV>cut && minJP>jpCut && jtSubJetSvtxm1>svtxmCut && jtSubJetSvtxm2>svtxmCut && jtSubJetdR>0.1) zgSpecTag_mc[ibin]->Fill(zg, weight);
		
			if(jtSubJetdR>0.1) drSpecIncl_mc[ibin]->Fill(jtSubJetdR, weight);
			if(jtSubJetRefdR>0.1) drSpecIncl_refAxis_mc[ibin]->Fill(jtSubJetRefdR, weight);
			if(minCSV>cut && minJP>jpCut && jtSubJetSvtxm1>svtxmCut && jtSubJetSvtxm2>svtxmCut && jtSubJetdR>0.1) drSpecTag_mc[ibin]->Fill(jtSubJetdR, weight);
		}
		for(int i=0; i<ptBins; i++){	
			drSpecIncl_mc[i]->Scale(1./drSpecIncl_mc[i]->Integral());
			drSpecIncl_refAxis_mc[i]->Scale(1./drSpecIncl_refAxis_mc[i]->Integral());
			drSpecTag_mc[i]->Scale(1./drSpecTag_mc[i]->Integral());
			drSpecBB_mc[i]->Scale(1./drSpecBB_mc[i]->Integral());
			
			drSpecIncl_mc[i]->Scale(1./drSpecIncl_mc[i]->GetBinWidth(4));
			drSpecIncl_refAxis_mc[i]->Scale(1./drSpecIncl_refAxis_mc[i]->GetBinWidth(4));
			drSpecTag_mc[i]->Scale(1./drSpecTag_mc[i]->GetBinWidth(4));
			drSpecBB_mc[i]->Scale(1./drSpecBB_mc[i]->GetBinWidth(4));
			
			zgSpecIncl_mc[i]->Scale(1./zgSpecIncl_mc[i]->Integral());
			zgSpecBB_mc[i]->Scale(1./zgSpecBB_mc[i]->Integral());
			//zgSpecTag_mc[i]->Scale(1./zgSpecTag_mc[i]->Integral());
			
			//zgSpecBB_mc_jp0[i]->Scale(1./zgSpecBB_mc_jp0[i]->Integral());
			zgSpecBB_mc_jp0p4[i]->Scale(1./zgSpecBB_mc_jp0p4[i]->Integral());
			zgSpecBB_mc_jp1[i]->Scale(1./zgSpecBB_mc_jp1[i]->Integral());
			
			zgSpecIncl_mc[i]->Scale(1./zgSpecIncl_mc[i]->GetBinWidth(4));
			zgSpecBB_mc[i]->Scale(1./zgSpecBB_mc[i]->GetBinWidth(4));
			//zgSpecTag_mc[i]->Scale(1./zgSpecTag_mc[i]->GetBinWidth(4));
			
			//zgSpecBB_mc_jp0[i]->Scale(1./zgSpecBB_mc_jp0[i]->GetBinWidth(4));
			zgSpecBB_mc_jp0p4[i]->Scale(1./zgSpecBB_mc_jp0p4[i]->GetBinWidth(4));
			zgSpecBB_mc_jp1[i]->Scale(1./zgSpecBB_mc_jp1[i]->GetBinWidth(4));
			
		}
	}
	
	
	if(!drawFromTree){
		for(int i=0; i<ptBins; i++){
			zgSpecIncl_data[i] = (TH1D*)fppdata->Get(Form("zgSpecIncl_data_%d",i))->Clone(Form("zgSpecIncl_data_%d",i)); //zgSpecIncl_data[i]->Sumw2();
			//zgSpecGSP_data[i] = (TH1D*)fppdata->Get(Form("zgSpecGSP_data_%d",i))->Clone(Form("zgSpecGSP_data_%d",i)); //zgSpecGSP_data[i]->Sumw2();
			zgSpecDR_data[i] = (TH1D*)fppdata->Get(Form("zgSpecDR_data%d",i))->Clone(Form("zgSpecDR_data%d",i)); //zgSpecDR_data[i]->Sumw2();
			//zgSpecTag_data[i] = (TH1D*)fppdata->Get(Form("zgSpecGSP_noDR_data%d",i))->Clone(Form("zgSpecGSP_noDR_data%d",i)); //zgSpecGSP_noDR_data[i]->Sumw2();
			
			zgSpecIncl_data[i]->Scale(1./zgSpecIncl_data[i]->Integral());
			//zgSpecGSP_data[i]->Scale(1./zgSpecGSP_data[i]->Integral());
			zgSpecDR_data[i]->Scale(1./zgSpecDR_data[i]->Integral());
			zgSpecTag_data[i]->Scale(1./zgSpecTag_data[i]->Integral());
			
			zgSpecIncl_data[i]->Scale(1./zgSpecIncl_data[i]->GetBinWidth(4));
			//zgSpecGSP_data[i]->Scale(1./zgSpecGSP_data[i]->GetBinWidth(4));
			zgSpecDR_data[i]->Scale(1./zgSpecDR_data[i]->GetBinWidth(4));
			zgSpecTag_data[i]->Scale(1./zgSpecTag_data[i]->GetBinWidth(4));
			
			zgSpecIncl_mc[i] = (TH1D*)fppmc->Get(Form("zgSpecIncl_mc%d",i))->Clone(Form("zgSpecIncl_mc%d",i));
			zgSpecBB_mc[i] = (TH1D*)fppmc->Get(Form("zgSpecGSP_mc%d",i))->Clone(Form("zgSpecBB_mc%d",i));
			
			zgSpecIncl_mc[i]->Scale(1./zgSpecIncl_mc[i]->Integral());
			//zgSpecBB_mc[i]->Scale(1./zgSpecBB_mc[i]->Integral());
			zgSpecIncl_mc[i]->Scale(1./zgSpecIncl_mc[i]->GetBinWidth(4));
			//zgSpecBB_mc[i]->Scale(1./zgSpecBB_mc[i]->GetBinWidth(4));
		}
	}
	
	//make correction factors to convert tagged -> true BB distributions
	
	TF1 *lightFrac[2];
	TF1 *mistagFrac[2];
	TF1 *lightFracData[2];
	TF1 *mistagFracData[2];
	
	for(int i=0; i<ptBins; i++){
	
		if(jpCut<0 && svtxmCut<0){
			lightFrac[i] = new TF1(Form("lightFrac_%d",i),"[0]+x*[1]",0.1,0.5);
			mistagFrac[i]  = new TF1(Form("mistagFrac_%d",i),"[0]+x*[1]",0.1,0.5);
			lightFracData[i] = new TF1(Form("lightFracData_%d",i),"[0]+x*[1]",0.1,0.5);
			mistagFracData[i] = new TF1(Form("mistagFracData_%d",i),"[0]+x*[1]",0.1,0.5);
			mcLightFrac[i]->Fit(lightFrac[i],"qN0","",0.1,0.5);
			mcPur[i]->Fit(mistagFrac[i],"qN0","",0.1,0.5);
			dataLightFrac[i]->Fit(lightFracData[i],"qN0","",0.1,0.5);
			dataPur[i]->Fit(mistagFracData[i],"qN0","",0.1,0.5);
		
			for(int ibin=1; ibin<=zgSpecTag_mc[i]->GetNbinsX(); ibin++){
				if(zgSpecTag_mc[i]->GetBinCenter(ibin)<0.1) continue;
			//do MC
				double lightF, mistagF;
				if(fitSFs){
					lightF = lightFrac[i]->Eval(zgSpecTag_mc[i]->GetBinCenter(ibin));
					mistagF = 1.-mistagFrac[i]->Eval(zgSpecTag_mc[i]->GetBinCenter(ibin));
				}  
				else{
					lightF = mcLightFrac[i]->GetBinContent(ibin);
					mistagF = mcPur[i]->GetBinContent(ibin);
				}
			      double corrBin = zgSpecTag_mc[i]->GetBinContent(ibin)*lightF/mistagF;// - zgSpecIncl_mc[i]->GetBinContent(ibin);
			      double corrErr = corrBin*sqrt(1./zgSpecTag_mc[i]->GetBinError(ibin) + 1./mcLightFrac[i]->GetBinError(ibin) + 1./mcPur[i]->GetBinError(ibin));
			      corrBin -= zgSpecIncl_mc[i]->GetBinContent(ibin);
			      corrErr = sqrt(pow(corrErr,2)+pow(zgSpecIncl_mc[i]->GetBinError(ibin),2));
			      zgSpecTagCorr_mc[i]->SetBinContent(ibin, corrBin);
			      zgSpecTagCorr_mc[i]->SetBinError(ibin, 0.01*corrErr);
			      
			//do Data
			      if(fitSFs){
			      	lightF = lightFracData[i]->Eval(zgSpecTag_mc[i]->GetBinCenter(ibin));
			      	mistagF = 1.-mistagFracData[i]->Eval(zgSpecTag_mc[i]->GetBinCenter(ibin));
			      }
			      else{
			      	lightF = dataLightFrac[i]->GetBinContent(ibin);
			      	mistagF = 1.-dataPur[i]->GetBinContent(ibin);
			      }
			      corrBin = zgSpecTag_data[i]->GetBinContent(ibin)*lightF/mistagF;
			      corrErr = corrBin*sqrt(1./zgSpecTag_data[i]->GetBinError(ibin) + 1./dataLightFrac[i]->GetBinError(ibin) + 1./dataPur[i]->GetBinError(ibin));
			      corrBin -= zgSpecIncl_data[i]->GetBinContent(ibin);
			      corrErr = sqrt(pow(corrErr,2)+pow(zgSpecIncl_data[i]->GetBinError(ibin),2));
			      zgSpecBB_data[i]->SetBinContent(ibin,corrBin);
			      zgSpecBB_data[i]->SetBinError(ibin,0.01*corrErr);
			}
			zgSpecTagCorr_mc[i]->Scale(1./zgSpecTagCorr_mc[i]->Integral());
			zgSpecBB_data[i]->Scale(1./zgSpecBB_data[i]->Integral());
			
			zgSpecTagCorr_mc[i]->Scale(1./zgSpecTagCorr_mc[i]->GetBinWidth(1));
			zgSpecBB_data[i]->Scale(1./zgSpecBB_data[i]->GetBinWidth(1));
		}
		else{
			for(int i=0; i<ptBins; i++){
				zgSpecBB_data[i] = (TH1D*)zgSpecTag_data[i]->Clone(Form("zgSpecBB_data_%d",i));
				zgSpecBB_data_shift[i] = (TH1D*)zgSpecTag_data_shift[i]->Clone(Form("zgSpecBB_data_shift_%d",i));
				zgSpecBB_data_v2[i] = (TH1D*)zgSpecTag_data[i]->Clone(Form("zgSpecBB_data_v2_%d",i));
				zgSpecTagCorr_mc[i] = (TH1D*)zgSpecTag_mc[i]->Clone(Form("zgSpecTagCorr_mc_%d",i));
			}
		}
	}
	
	
	//now correct for the tagging inefficiencies at low-zg
	TH1D *tagEffMC[ptBins];
	TH1D *tagEffData[ptBins];
	TH1D *tagEffShift[ptBins];
	TF1 *fTagCorrFactor[ptBins];
	TF1 *fTagCorrFactorSys[ptBins];
	TF1 *fTagCorrFactorShift[ptBins];
	TLatex *ll2[2];
	TCanvas *cctemp = new TCanvas("cctemp","",1200,600);
	cctemp->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		tagEffMC[i] = (TH1D*)mcEff[i]->Clone(Form("tagEffMC_%d",i));
		//tagEffMC[i] = (TH1D*)zgSpecTag_mc[i]->Clone(Form("tagEffMC_%d",i));
		//tagEffMC[i]->Divide(zgSpecBB_mc_jp0[i]);
		fTagCorrFactor[i] = new TF1(Form("fTagCorrFactor_pt%d",i),"[0]+[1]*TMath::Exp(-1*[2]*x)",0.1,0.5);
		//fTagCorrFactor[i] = new TF1(Form("fTagCorrFactor_pt%d",i),"pol2",0.2,0.5);
		fTagCorrFactor[i]->SetParLimits(1,-5,5);
		fTagCorrFactor[i]->SetParLimits(2,0.1,20);
		cctemp->cd(i+1);
		tagEffMC[i]->Fit(fTagCorrFactor[i],"","",0.1,0.5);
		fTagCorrFactorSys[i] = (TF1*)fTagCorrFactor[i]->Clone(Form("fTagCorrFactorSys_pt%d",i));
		//fTagCorrFactorSys[i]->SetParLimits(1,-5,5);
		//fTagCorrFactorSys[i]->SetParLimits(2,0.1,20);
		tagEffData[i] = (TH1D*)dataEff[i]->Clone(Form("tagEffData_%d",i));
		//tagEffData[i]->Multiply(scaleFactor);
		//tagEffData[i]->Scale(tagEffMC[i]->Integral()/tagEffData[i]->Integral());
		tagEffData[i]->Fit(fTagCorrFactorSys[i],"","",0.1,0.5);
		//fTagCorrFactorSys[i]->SetParameter(0, fTagCorrFactor[i]->GetParameter(0));
		tagEffMC[i]->SetMarkerStyle(20);
		tagEffMC[i]->SetMarkerColor(1);
		tagEffMC[i]->SetXTitle("zG");
		tagEffMC[i]->SetYTitle("Tagging Efficiency");
		tagEffMC[i]->SetMaximum(0.3);
		tagEffMC[i]->Draw();
		tagEffData[i]->SetMarkerStyle(21);
		tagEffData[i]->SetMarkerColor(2);
		tagEffData[i]->Draw("same");
		
		fTagCorrFactorShift[i] = (TF1*)fTagCorrFactor[i]->Clone(Form("fTagCorrFactorShift_pt%d",i));
		tagEffShift[i] = (TH1D*)jesEff[i]->Clone(Form("tagEffDataShift_%d",i));
		tagEffShift[i]->Fit(fTagCorrFactorShift[i],"qN0","",0.1,0.5);
		
		if(i<ptBins-1) ll2[i] = new TLatex(0.05,0.2,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		else ll2[i] = new TLatex(0.05,0.2,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		ll2[i]->Draw("Same");
		
		cout << "ptBin " << i+1 << endl;
		cout << "dataEff Integral: " << tagEffData[i]->Integral() << endl;
		cout << "MC Integral: " << tagEffMC[i]->Integral() << endl;
	}
	for(int i=0; i<ptBins; i++){
		
		if(fitSFs){
			correctSubjetReco(fTagCorrFactor[i], zgSpecBB_data[i], 1);
			correctSubjetReco(fTagCorrFactor[i], zgSpecBB_data_shift[i],1); //do we correct the shift properly? or just take the shape difference?
			correctSubjetReco(fTagCorrFactor[i], zgSpecTagCorr_mc[i], 1);
		}
		else{
			//zgSpecBB_data[i]->Divide(tagEffData[i]);
			//zgSpecTagCorr_mc[i]->Divide(tagEffData[i]);
		}
		
		zgSpecTag_mc[i]->Scale(1./zgSpecTag_mc[i]->Integral());
		zgSpecTag_mc[i]->Scale(1./zgSpecTag_mc[i]->GetBinWidth(1));
		
		zgSpecBB_data[i]->Scale(1./zgSpecBB_data[i]->Integral());
		zgSpecBB_data[i]->Scale(1./zgSpecBB_data[i]->GetBinWidth(1));
		
		zgSpecBB_data_shift[i]->Scale(1./zgSpecBB_data_shift[i]->Integral());
		zgSpecBB_data_shift[i]->Scale(1./zgSpecBB_data_shift[i]->GetBinWidth(1));
		
		zgSpecBB_mc_jp0[i]->Scale(1./zgSpecBB_mc_jp0[i]->Integral());
		zgSpecBB_mc_jp0[i]->Scale(1./zgSpecBB_mc_jp0[i]->GetBinWidth(1));
		
	}
	
	
	//now correct for the subjet reco inefficiencies at low-zg
	//cout << "correcting reco inefficiency" << endl;
	TFile *fIneff = new TFile("partonCheck_pythia8.root");
	TH1D *trueZg[ptBins];
	TH1D *subjetZg[ptBins];
	TF1 *fCorrFactor[ptBins];
	for(int i=0; i<ptBins; i++){
		trueZg[i] = (TH1D*)fIneff->Get("gspZg")->Clone(Form("trueZg_%d",i));
		subjetZg[i] = (TH1D*)fIneff->Get(Form("matchedPartonZg_%d",i))->Clone(Form("subjetZg_%d",i));
		trueZg[i]->Divide(subjetZg[i]);
		fCorrFactor[i] = new TF1(Form("fCorrFactor_pt%d",i),"[0]+[1]/x+[2]/(x*x)",0.1,0.5);
		trueZg[i]->Fit(fCorrFactor[i],"qN0","",0.1,0.5);
		
		subjetZg[i]->Scale(1./subjetZg[i]->Integral());
		subjetZg[i]->Scale(1./subjetZg[i]->GetBinWidth(1));
	}
	
	//Have to do the first dR bin as well...
	TF1 *fDRCorrFactor[ptBins];
	TF1 *fDRFlatCorrFactor[ptBins];
	TH1D *drEffHisto[ptBins];
	TFile *drTagEff = new TFile("drMCTagEff.root");
	for(int i=0; i<ptBins; i++){
		fDRCorrFactor[i] = new TF1(Form("fDRCorrFactor_%d",i),"pol1",0.1,0.4);
		fDRFlatCorrFactor[i] = new TF1(Form("fDRFlatCorrFactor_%d",i),"pol0",0.1,0.4);
		drEffHisto[i] = (TH1D*)drTagEff->Get(Form("tagged_%d",i));
		drEffHisto[i]->Fit(fDRCorrFactor[i],"qN0","",0.1,0.4);
		drEffHisto[i]->Fit(fDRFlatCorrFactor[i],"qN0","",0.1,0.4);
		
		correctSubjetReco(fDRCorrFactor[i], drSpecGSP_data[i], 1);
		correctSubjetReco(fDRCorrFactor[i], drSpecGSP_data_shift[i], 1);
		drSpecGSP_data[i]->Scale(1./drSpecGSP_data[i]->Integral());
		drSpecGSP_data[i]->Scale(1./drSpecGSP_data[i]->GetBinWidth(1));
		drSpecGSP_data_shift[i]->Scale(1./drSpecGSP_data_shift[i]->Integral());
		drSpecGSP_data_shift[i]->Scale(1./drSpecGSP_data_shift[i]->GetBinWidth(1));
	}
	
	for(int i=0; i<ptBins; i++){
		//correctSubjetReco(fCorrFactor[i], zgSpecBB_data[i]);
		//correctSubjetReco(fCorrFactor[i], zgSpecBB_data_v2[i]);
		//correctSubjetReco(fCorrFactor[i], zgSpecTagCorr_mc[i]);
		//correctSubjetReco(fCorrFactor[i], zgSpecBB_mc[i]);
		//correctSubjetReco(fCorrFactor[i], zgSpecBB_mc_jp0[i]);
		
		zgSpecTagCorr_mc[i]->Scale(1./zgSpecTagCorr_mc[i]->Integral());
		zgSpecBB_data[i]->Scale(1./zgSpecBB_data[i]->Integral());
		zgSpecBB_data_v2[i]->Scale(1./zgSpecBB_data_v2[i]->Integral());
		zgSpecBB_mc[i]->Scale(1./zgSpecBB_mc[i]->Integral());
		zgSpecBB_mc_jp0[i]->Scale(1./zgSpecBB_mc_jp0[i]->Integral());
		
		zgSpecTagCorr_mc[i]->Scale(1./zgSpecTagCorr_mc[i]->GetBinWidth(1));
		zgSpecBB_data[i]->Scale(1./zgSpecBB_data[i]->GetBinWidth(1));
		zgSpecBB_data_v2[i]->Scale(1./zgSpecBB_data_v2[i]->GetBinWidth(1));
		zgSpecBB_mc[i]->Scale(1./zgSpecBB_mc[i]->GetBinWidth(1));
		zgSpecBB_mc_jp0[i]->Scale(1./zgSpecBB_mc_jp0[i]->GetBinWidth(1));
		
	}
	cout << "done" << endl;
	
	TFile *fEffCorrs = new TFile("zgCorrFactors_gsp.root","recreate");
	fEffCorrs->cd();
	for(int i=0; i<ptBins; i++){
		fCorrFactor[i]->Write();
		fTagCorrFactor[i]->Write();
	}
	
	TFile *fparton = new TFile("partonCheck_fullQCDMC_eta1p5_fullJetPtBins.root");
	TH1D *partonZg = (TH1D*)fparton->Get("gspZg");
	partonZg->Scale(1./partonZg->Integral());
	partonZg->Scale(1./partonZg->GetBinWidth(1));

	// **************** START SYSTEMATIC UNCERTAINTIES HERE ******************//
	
	double purUnc[10][2], angRes[10][2], jesUnc[10][2], jerUnc[10][2], subjetTagEff[10][2];
	double jesUncIncl[10][2];
	double jesNonclosureIncl[10][2], jesNonclosureB[10][2];
	
	double purUncDR[20][2], jesUncDR[20][2], jerUncDR[20][2];
	double jesUncInclDR[20][2], jerUncInclDR[20][2];
	double jesNonclosureInclDR[20][2], jesNonclosureBDR[20][2];
	double subjetTagEffDR[20][2];
	double angResDR[20][2];
	
	TFile *JESRatios = new TFile("dR_qcdJet_JESRatios.root");
	TFile *JERRatios = new TFile("dR_qcdJet_JERRatios.root");
	TH1F *jesUncInclDRHistoUp[2], *jesUncInclDRHistoDown[2], *jerUncInclDRHisto[2], *jerUncIncl[2];
	jesUncInclDRHistoUp[0] = (TH1F*)JESRatios->Get("h_dR12_qcdpt100140_JESupRatio39");
	jesUncInclDRHistoDown[0] = (TH1F*)JESRatios->Get("h_dR12_qcdpt100140_JESdownRatio40");
	jesUncInclDRHistoUp[1] = (TH1F*)JESRatios->Get("h_dR12_qcdpt140up_JESupRatio41");
	jesUncInclDRHistoDown[1] = (TH1F*)JESRatios->Get("h_dR12_qcdpt140up_JESdownRatio42");
	jerUncInclDRHisto[0] = (TH1F*)JERRatios->Get("h_dR12_qcdpt100140_JERgausRatio57");
	jerUncInclDRHisto[1] = (TH1F*)JERRatios->Get("h_dR12_qcdpt140up_JERgausRatio58");
	
	TFile *inclJERRatios = new TFile("zg_qcdJER.root");
	jerUncIncl[0] = (TH1F*)inclJERRatios->Get("h_Zg_qcdpt120160_JERgausRatio25");
	jerUncIncl[1] = (TH1F*)inclJERRatios->Get("h_Zg_qcdpt160up_JERgausRatio26");
	
	TF1 *JESB[2], *JESIncl[2];
	for(int ptBin=0; ptBin<2; ptBin++){
		JESB[ptBin] = new TF1(Form("JESB_%d",ptBin),"pol1",0.1,0.5);
		JESIncl[ptBin] = new TF1(Form("JESIncl_%d",ptBin),"pol1",0.1,0.5);
		zgSpecIncl_mc_JES[ptBin]->Divide(zgSpecIncl_mc_JES[ptBin], zgSpecIncl_mc[ptBin],1,1,"B");
		zgSpecBB_mc_JES[ptBin]->Divide(zgSpecBB_mc_JES[ptBin], zgSpecBB_mc[ptBin],1,1,"B");
		
		zgSpecBB_mc_JES[ptBin]->Fit(JESB[ptBin],"qN0","",0.1,0.5);
		zgSpecIncl_mc_JES[ptBin]->Fit(JESIncl[ptBin],"qN0","",0.1,0.5);
		
		for(int zgBin=0; zgBin<10; zgBin++){

			//purUnc[zgBin][ptBin] = 0.03;
			angRes[zgBin][ptBin] = 0.04;
		
			//jesUncIncl[zgBin][ptBin] = abs(zgSpecIncl_mc_JES[ptBin]->GetBinContent(zgBin+1)-1.);
			//jesUnc[zgBin][ptBin] = abs(zgSpecBB_mc_JES[ptBin]->GetBinContent(zgBin+1)-1.);
			jesUncIncl[zgBin][ptBin] = abs(JESIncl[ptBin]->Eval(zgBin*0.05+0.025)-1.);
			jesUnc[zgBin][ptBin] = abs(JESB[ptBin]->Eval(zgBin*0.05+0.025)-1.);
			
			jerUncIncl[ptBin]->SetBinContent(zgBin+1,1.02);
		//}
		//inclusive jets...
		//jesUncIncl[zgBin][0] = 0.08;
		//jesUncIncl[zgBin][1] = 0.08;
			
		//gsp jets...
		//jesUnc[zgBin][0] = 0.05;
		//jesUnc[zgBin][1] = 0.10;
			
			if(ptBin==0) jerUnc[zgBin][ptBin] = 0.07;
			else jerUnc[zgBin][ptBin] = 0.05;
		}
	}
	
	for(int ptBin=0; ptBin<2; ptBin++){
		drSpecIncl_mc_JES[ptBin]->Divide(drSpecIncl_mc_JES[ptBin], drSpecIncl_mc[ptBin],1,1,"B");
		
		for(int idrbin=0; idrbin<20; idrbin++){
			jesUncDR[idrbin][ptBin] = 0.15;
			//jesUncDR[idrbin][ptBin] = abs(drSpecIncl_mc_JES[ptBin]->GetBinContent(idrbin+1)-1.);
			jesUncInclDR[idrbin][ptBin] = abs(drSpecIncl_mc_JES[ptBin]->GetBinContent(idrbin+1)-1.);
			
			jerUncDR[idrbin][ptBin] = 0.12;
			
			double binCenter = jerUncInclDRHisto[ptBin]->GetBinCenter(idrbin);
			
			//jesUncInclDR[idrbin][ptBin] = abs(TMath::Max(jesUncInclDRHistoUp[ptBin]->GetBinContent(idrbin+1), jesUncInclDRHistoDown[ptBin]->GetBinContent(idrbin+1))-1.);
			//jerUncInclDR[idrbin][ptBin] = abs(jerUncInclDRHisto[ptBin]->GetBinContent(idrbin+1)-1.);
			jerUncInclDR[idrbin][ptBin] = 0.02;
			
			subjetTagEffDR[idrbin][ptBin] = abs(fDRCorrFactor[ptBin]->Eval(binCenter)/fDRFlatCorrFactor[ptBin]->Eval(binCenter) - 1.);
			
		}
		
		/*jesUncInclDR[idrbin][0] = 0.05;
		jesUncInclDR[idrbin][1] = 0.05;
		
		jerUncInclDR[idrbin][0] = 0.05;
		jerUncInclDR[idrbin][1] = 0.05;*/
		
	}
	
	TF1 *JPsysFit[2];
	TF1 *DMesonSysFit[2];
	TF1 *MCNonclosure[2];
	TF1 *MCBNonclosure[2];
	TH1D *InclSpecPurSys[2], *InclSpecPurSysDR[2];
	double jpCalibSys[10][2];
	double DMesonSys[10][2];
	for(int i=0; i<ptBins; i++){
		
		zgSpecBB_data_shift[i]->Divide(zgSpecBB_data_shift[i],zgSpecBB_data[i],1,1,"B");
		zgSpecIncl_data_shift[i]->Divide(zgSpecIncl_data_shift[i],zgSpecIncl_data[i],1,1,"B");
		
		drSpecGSP_data_shift[i]->Divide(drSpecGSP_data_shift[i],drSpecGSP_data[i],1,1,"B");
		drSpecIncl_data_shift[i]->Divide(drSpecIncl_data_shift[i],drSpecIncl_data[i],1,1,"B");
		
		drSpecIncl_refAxis_mc[i]->Divide(drSpecIncl_mc[i],drSpecIncl_refAxis_mc[i],1,1,"B");
		
		InclSpecPurSys[i] = (TH1D*)zgSpecIncl_data[i]->Clone(Form("InclSpecPurSys_%d",i));
		InclSpecPurSys[i]->Scale(0.10); //using 10% nonpurity
		InclSpecPurSys[i]->Add(zgSpecBB_data[i]);
		InclSpecPurSys[i]->Scale(1./InclSpecPurSys[i]->Integral());
		InclSpecPurSys[i]->Scale(1./InclSpecPurSys[i]->GetBinWidth(1));
		InclSpecPurSys[i]->Divide(InclSpecPurSys[i],zgSpecBB_data[i],1,1,"B");
		
		InclSpecPurSysDR[i] = (TH1D*)drSpecIncl_data[i]->Clone(Form("InclSpecPurSysDR_%d",i));
		InclSpecPurSysDR[i]->Scale(0.10); //using 10% nonpurity
		InclSpecPurSysDR[i]->Add(drSpecGSP_data[i]);
		InclSpecPurSysDR[i]->Scale(1./InclSpecPurSysDR[i]->Integral());
		InclSpecPurSysDR[i]->Scale(1./InclSpecPurSysDR[i]->GetBinWidth(1));
		InclSpecPurSysDR[i]->Divide(InclSpecPurSysDR[i],drSpecGSP_data[i],1,1,"B");
		
		JPsysFit[i] = new TF1(Form("JPsysFit_%d",i),"pol1",0.1,0.5);
		DMesonSysFit[i] = new TF1(Form("DMesonSysFit_%d",i),"pol1",0.1,0.5);
		MCNonclosure[i] = new TF1(Form("MCNonclosure_%d",i),"pol1",0.1,0.5);
		MCBNonclosure[i] = new TF1(Form("MCBNonclosure_%d",i),"pol1",0.1,0.5);
	//mcEffJPsys->Divide(mcEff);
		if(i==0) mcEffJPSys->Add(mcEff[i],-1);
		mcEffJPSys->Fit(JPsysFit[i],"qN0","",0.1,0.5);
		if(i==0) mcEffDMeson->Add(mcEff[i],-1);
		mcEffDMeson->Fit(DMesonSysFit[i],"qN0","",0.1,0.5);
		zgSpecIncl_data_shift[i]->Fit(MCNonclosure[i],"qN0","",0.1,0.5);
		zgSpecBB_data_shift[i]->Fit(MCBNonclosure[i],"qN0","",0.1,0.5);
		for(int ibin=0; ibin<2; ibin++){ jpCalibSys[ibin][i]=0.; DMesonSys[ibin][i]=0.; jesNonclosureIncl[ibin][i]=0.; jesNonclosureB[ibin][i]=0.; purUnc[ibin][i]=0.;}
		for(int ibin=2; ibin<10; ibin++){
			jpCalibSys[ibin][i] = abs(JPsysFit[i]->Eval(0.05*ibin))/mcEff[i]->GetBinContent(mcEff[i]->FindBin(0.05*ibin));
			DMesonSys[ibin][i] = abs(DMesonSysFit[i]->Eval(0.05*ibin))/mcEff[i]->GetBinContent(mcEff[i]->FindBin(0.05*ibin));
			purUnc[ibin][i] = abs(InclSpecPurSys[i]->GetBinContent(ibin+1)-1.);
			//jesNonclosureIncl[ibin][i] = zgSpecIncl_data_shift[i]->GetBinContent(ibin+1)-1.;
			//jesNonclosureB[ibin][i] = zgSpecBB_data_shift[i]->GetBinContent(ibin+1)-1.;
			jesNonclosureIncl[ibin][i] =  MCNonclosure[i]->Eval(0.05*ibin+0.025)-1.;
			jesNonclosureB[ibin][i] =  MCBNonclosure[i]->Eval(0.05*ibin+0.025)-1.;
		}
		
		for(int ibin=0; ibin<20; ibin++){
			angResDR[ibin][i] = abs(drSpecIncl_refAxis_mc[i]->GetBinContent(ibin+1)-1.);
			purUncDR[ibin][i] = abs(InclSpecPurSysDR[i]->GetBinContent(ibin+1)-1.);
			jesNonclosureInclDR[ibin][i] = abs(drSpecIncl_data_shift[i]->GetBinContent(ibin+1)-1.);
			jesNonclosureBDR[ibin][i] = abs(drSpecGSP_data_shift[i]->GetBinContent(ibin+1)-1.);
		}
	}
	
	// Fill TGraphs of Systematic uncertainties
	
	TGraphErrors *bjetUnc[ptBins];
	TGraphErrors *drUnc[ptBins];
	TGraphErrors *inclJetUnc[ptBins];
	TGraphErrors *drInclUnc[ptBins];
	
	TGraphErrors *normUncert[ptBins];
	for(int i=0; i<ptBins; i++){
		bjetUnc[i] = new TGraphErrors(zgSpecBB_data[i]);
		inclJetUnc[i] = new TGraphErrors(zgSpecIncl_data[i]);
		drUnc[i] = new TGraphErrors(drSpecGSP_data[i]);
		drInclUnc[i] = new TGraphErrors(drSpecIncl_data[i]);
		inclJetUnc[i]->SetFillColor(kGray);
		inclJetUnc[i]->SetFillStyle(3002);
		inclJetUnc[i]->SetMarkerColor(1);
		inclJetUnc[i]->SetLineColor(1);
		bjetUnc[i]->SetFillColor(kMagenta+1);
		bjetUnc[i]->SetFillStyle(3002);
		bjetUnc[i]->SetMarkerColor(kMagenta+2);
		bjetUnc[i]->SetLineColor(kMagenta+2);
		drUnc[i]->SetFillColor(kMagenta+1);
		drUnc[i]->SetFillStyle(3002);
		drUnc[i]->SetMarkerColor(kMagenta+2);
		drUnc[i]->SetLineColor(kMagenta+2);
		drInclUnc[i]->SetFillColor(kGray);
		drInclUnc[i]->SetFillStyle(3002);
		drInclUnc[i]->SetMarkerColor(1);
		drInclUnc[i]->SetLineColor(1);
		for(int ibin=0; ibin<zgSpecBB_data[i]->GetNbinsX(); ibin++){
			
			//take full difference between data and MC?
			//subjetTagEff[ibin][i] = abs(fTagCorrFactorSys[i]->Eval(ibin*0.05+0.025)/fTagCorrFactor[i]->Eval(ibin*0.05+0.025) - 1.);
			
			//Instead lets just use the proper data-driven errors
			int translateBin = tagEffData[i]->FindBin(zgSpecBB_data[i]->GetBinCenter(ibin+1));
			if(tagEffData[i]->GetBinContent(translateBin)>0) subjetTagEff[ibin][i] = tagEffData[i]->GetBinError(translateBin)/tagEffData[i]->GetBinContent(translateBin);
			else subjetTagEff[ibin][i] = 0.;
			
			int jerUncInclBin = jerUncIncl[i]->FindBin(zgSpecBB_data[i]->GetBinCenter(ibin+1));
			
			double totErr = sqrt(pow(purUnc[ibin][i],2)+pow(angRes[ibin][i],2)+pow(jesUnc[ibin][i],2)+pow(jerUnc[ibin][i],2)+pow(subjetTagEff[ibin][i],2)+pow(jpCalibSys[ibin][i],2)+pow(DMesonSys[ibin][i],2)+pow(jesNonclosureB[ibin][i],2));
			double incErr = sqrt(pow(angRes[ibin][i],2)+pow(jesUncIncl[ibin][i],2)+pow(jerUncIncl[i]->GetBinContent(jerUncInclBin)-1.,2)+pow(jesNonclosureIncl[ibin][i],2));
			bjetUnc[i]->SetPointError(ibin, 0.025, totErr*zgSpecBB_data[i]->GetBinContent(ibin+1));
			inclJetUnc[i]->SetPointError(ibin, 0.025, incErr*zgSpecIncl_data[i]->GetBinContent(ibin+1));
		}
		
		//now do the normalization uncertainty
		double maxUncert = 0;
		for(int ibin=0; ibin<zgSpecBB_data[i]->GetNbinsX(); ibin++){
			//take the max difference 
			if(maxUncert < bjetUnc[i]->GetErrorY(ibin)) maxUncert = bjetUnc[i]->GetErrorY(ibin);
		}
		normUncert[i] = new TGraphErrors(1);
		normUncert[i]->SetPoint(1, 0.02, 2.2);  //not really sure where to put this guy... 2.2 seems reasonable
		normUncert[i]->SetPointError(1, 0.01, maxUncert/7./2.);
		
		for(int ibin=0; ibin<drSpecGSP_data[i]->GetNbinsX(); ibin++){
			double drTotErr = sqrt(pow(subjetTagEffDR[ibin][i],2)+pow(purUncDR[ibin][i],2)+pow(jesUncDR[ibin][i],2)+pow(jerUncDR[ibin][i],2)+pow(jesNonclosureBDR[ibin][i],2)+pow(angResDR[ibin][i],2));
			double drInclErr = sqrt(pow(jesUncInclDR[ibin][i],2)+pow(jerUncInclDR[ibin][i],2)+pow(jesNonclosureInclDR[ibin][i],2)+pow(angResDR[ibin][i],2));
			drUnc[i]->SetPointError(ibin, 0.01, drTotErr*drSpecGSP_data[i]->GetBinContent(ibin+1));
			drInclUnc[i]->SetPointError(ibin, 0.01, drInclErr*drSpecIncl_data[i]->GetBinContent(ibin+1));
		}
	}
	
	///*****  PLOT ZG RESULTS ******************//
	
	TLatex *labels[ptBins];
	TLegend *leg1 = new TLegend(0.5,0.62,0.88,0.77);
	TLegend *legNorm = new TLegend(0.22, 0.19, 0.67, 0.242);
	TCanvas *cc = new TCanvas("cc","",1200,600);
	cc->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cc->cd(i+1);
		zgSpecIncl_data[i]->SetYTitle("1/N dN/dz_{g}");
		zgSpecIncl_data[i]->SetXTitle("z_{g}");
		zgSpecIncl_data[i]->SetFillColor(kGray);
		zgSpecIncl_data[i]->SetFillStyle(3002);
		zgSpecIncl_data[i]->GetYaxis()->SetRangeUser(0,10);
		zgSpecIncl_data[i]->Draw();
		for(int ibin=1; ibin<=zgSpecIncl_data[i]->GetNbinsX(); ibin++){
			//cout << "zg " << 0.05*(ibin-1) << "-"<< 0.05*(ibin) << " " << zgSpecIncl_data[i]->GetBinContent(ibin) << " +/- " << zgSpecIncl_data[i]->GetBinError(ibin) << endl;
		}
		zgSpecTag_data[i]->SetLineColor(kOrange+1);
		zgSpecTag_data[i]->SetMarkerColor(kOrange+1);
		//zgSpecTag_data[i]->Draw("same");
		
		zgSpecBB_data[i]->SetYTitle("1/N dN/dz_{g}");
		zgSpecBB_data[i]->SetXTitle("z_{g}");
		zgSpecBB_data[i]->SetFillColor(kGray);
		zgSpecBB_data[i]->SetFillStyle(3002);
		zgSpecBB_data[i]->GetYaxis()->SetRangeUser(0,10);
		zgSpecBB_data[i]->SetLineColor(kMagenta+2);
		zgSpecBB_data[i]->SetMarkerColor(kMagenta+2);
		zgSpecBB_data[i]->Draw("");
		bjetUnc[i]->Draw("2,5,same");
		zgSpecBB_data[i]->Draw("same");
		inclJetUnc[i]->Draw("2,5,same");
		zgSpecIncl_data[i]->Draw("same");
		
		zgSpecTagCorr_mc[i]->SetYTitle("1/N dN/dz_{g}");
		zgSpecTagCorr_mc[i]->SetXTitle("z_{g}");
		zgSpecTagCorr_mc[i]->SetLineColor(4);
		zgSpecTagCorr_mc[i]->SetMarkerColor(4);
		//zgSpecTagCorr_mc[i]->Draw("same");
		zgSpecTagCorr_mc[i]->GetYaxis()->SetRangeUser(0,10);
		zgSpecBB_mc[i]->SetLineColor(4);
		zgSpecBB_mc[i]->SetMarkerColor(4);
		zgSpecBB_mc[i]->SetMarkerStyle(25);
		zgSpecBB_mc[i]->Draw("same");
		partonZg->SetLineColor(2);
		partonZg->SetMarkerColor(2);
		partonZg->SetMarkerStyle(25);
		//partonZg->Draw("same");
		subjetZg[i]->SetLineColor(4);
		subjetZg[i]->SetMarkerColor(4);
		subjetZg[i]->SetMarkerStyle(25);
		//subjetZg[i]->Draw("same");
		
		zgSpecPowheg[i]->SetLineColor(kGreen+1);
		zgSpecPowheg[i]->SetMarkerColor(kGreen+1);
		zgSpecPowheg[i]->Draw("same");
		for(int ibin=1; ibin<=zgSpecPowheg[i]->GetNbinsX(); ibin++){
			cout << "bin : " << ibin << " pwhg content: " << zgSpecPowheg[i]->GetBinContent(ibin) << " +/- " << zgSpecPowheg[i]->GetBinError(ibin) << endl;
		}
		
		normUncert[i]->SetFillColor(kGreen-9);
		normUncert[i]->SetLineColor(kGreen-9);
		normUncert[i]->SetFillStyle(3002);
		normUncert[i]->Draw("2,5,same");
		
		
		
		labels[i] = new TLatex(0.135,9.0,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		//else labels[i] = new TLatex(0.220,9.0,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		labels[i]->Draw("Same");
		
		if(i==0){
			leg1->AddEntry(zgSpecIncl_data[i], "Inclusive pp", "lpf");
			//leg1->AddEntry(zgSpecBB_data_v2[i], "BB (GSP) data, no tagEff Corr","lp");
			leg1->AddEntry(bjetUnc[i], "pp g->bb", "lpf");
			//leg1->AddEntry(zgSpecTagCorr_mc[i], "PYTHIA 8 g->bb", "lp"); //tagged and corrected
			leg1->AddEntry(zgSpecBB_mc[i], "PYTHIA 8 g->bb", "lp"); //flavor==5 USE THIS ONE!!
			//leg1->AddEntry(partonZg, "Parton-level", "lp");
			//leg1->AddEntry(zgSpecTag_mc[i], "CSV>0.8 + BB MC", "lp");
			//leg1->AddEntry(zgSpecIncl_mc[i],"Inclusive Pythia 8","lp");
			//leg1->AddEntry(subjetZg[i], "Parton-level Pythia 8 g->bb","lp"); //parton-level (not a fair comparison...)
			leg1->AddEntry(zgSpecPowheg[i], "POWHEG g->bb");
			leg1->Draw("same");
			
			legNorm->AddEntry(normUncert[i], "Normalization Uncert.", "f");
			legNorm->Draw("same");
		}
	}
	
	// ************** PLOT DR RESULTS *****************//
	
	TLegend *leg2 = new TLegend(0.217,0.643,0.613,0.803);
	TLatex *labelsdr[ptBins];
	TCanvas *ccdr = new TCanvas("ccdr","",1200,600);
	ccdr->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		ccdr->cd(i+1);
		drSpecIncl_data[i]->SetYTitle("1/N dN/d#DeltaR_{1,2}");
		drSpecIncl_data[i]->SetXTitle("#DeltaR_{1,2}");
		drSpecIncl_data[i]->GetYaxis()->SetRangeUser(0,10);
		drSpecIncl_data[i]->SetFillColor(kGray);
		drSpecIncl_data[i]->SetFillStyle(3002);
		drSpecIncl_data[i]->Draw();
		drSpecGSP_data[i]->SetLineColor(kMagenta+2);
		drSpecGSP_data[i]->SetMarkerColor(kMagenta+2);
		
		drSpecGSP_data[i]->SetYTitle("1/N dN/d#DeltaR_{1,2}");
		drSpecGSP_data[i]->SetXTitle("#DeltaR_{1,2}");
		drSpecGSP_data[i]->GetYaxis()->SetRangeUser(0,10);
		drSpecGSP_data[i]->SetFillColor(kGray);
		drSpecGSP_data[i]->SetFillStyle(3002);
		drSpecGSP_data[i]->Draw("");
		drUnc[i]->Draw("2,5,same");
		drSpecGSP_data[i]->Draw("same");
		drInclUnc[i]->Draw("2,5,same");
		drSpecIncl_data[i]->Draw("same");
		
		drSpecIncl_mc[i]->SetYTitle("1/N dN/d#DeltaR_{1,2}");
		drSpecIncl_mc[i]->SetXTitle("#DeltaR_{1,2}");
		drSpecIncl_mc[i]->GetYaxis()->SetRangeUser(0,10);
		drSpecIncl_mc[i]->SetLineColor(2);
		drSpecIncl_mc[i]->SetMarkerColor(2);
		//drSpecIncl_mc[i]->Draw("");
		
		drSpecBB_mc[i]->SetLineColor(4);
		drSpecBB_mc[i]->SetMarkerColor(4);
		drSpecBB_mc[i]->SetMarkerStyle(25);
		drSpecBB_mc[i]->Draw("same");
		
		drSpecPowheg[i]->SetLineColor(kGreen+1);
		drSpecPowheg[i]->SetMarkerColor(kGreen+1);
		drSpecPowheg[i]->Draw("same");
		
		labelsdr[i] = new TLatex(0.025,9.0,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		labelsdr[i]->Draw("Same");
		
		if(i==0){
			leg2->AddEntry(drSpecIncl_data[i], "Inclusive pp", "lpf");
			leg2->AddEntry(drUnc[i], "pp g->bb", "lpf");
			leg2->AddEntry(drSpecBB_mc[i], "PYTHIA 8 g->bb", "lp");
			//leg2->AddEntry(drSpecIncl_mc[i], "Inclusive Pythia 8", "lp");
			leg2->AddEntry(drSpecPowheg[i], "POWHEG g->bb", "lp");
			leg2->Draw("same");
		}
		
		//cout << "ptBin " << i << endl;
		//cout << "inclusive avg dr: "<< drSpecIncl_data[i]->GetMean() << " +/- " << drInclUnc[i]->GetRMS() << endl;
		//cout << "b-jet avg dr: " << drSpecGSP_data[i]->GetMean() << " +/- " << drUnc[i]->GetRMS() << endl;
	}
	
	// ****************** PLOT TAGGING EFFICIENCY CORRECTION *****************//
	
	/*TLatex *labels2[ptBins];
	TH1D *htmpForDiv[ptBins];
	TCanvas *cc2 = new TCanvas("cc2","",1200,600);
	cc2->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cc2->cd(i+1);
		htmpForDiv[i] = (TH1D*)zgSpecTag_data[i]->Clone(Form("htmpForDiv_%d",i));
		htmpForDiv[i]->Divide(zgSpecBB_data[i]);
		htmpForDiv[i]->SetYTitle("CSV-Tagged / Corrected z_{G}");
		htmpForDiv[i]->SetXTitle("z_{G}");
		htmpForDiv[i]->Draw();
		htmpForDiv[i]->GetYaxis()->SetRangeUser(0.5,2.0);
		
		if(i<ptBins-1) labels2[i] = new TLatex(0.25,0.8,Form("%d < Jet pT < %d",jetPtBin[i],jetPtBin[i+1]));
		else labels2[i] = new TLatex(0.25,0.8,Form("Jet pT > %d",jetPtBin[i]));
		labels2[i]->Draw("Same");

	}*/
	
	//***************** ZG SYSTEMATIC PLOT (GSP JETS) *******************//
	
	TH1D *purUncHisto[ptBins], *angResHisto[ptBins], *jetUncHisto[ptBins], *jerUncHisto[ptBins], *subjetTagUncHisto[ptBins], *jpCalibUncHisto[ptBins], *dMesonUncHisto[ptBins], *jesNonclosureHisto[ptBins], *totalUnc[ptBins];
	TLegend *uncLeg = new TLegend(0.19,0.64,0.9,0.93);
	TCanvas *cSys = new TCanvas("cSys","",1200,600);
	TLatex *labelsSys[2];
	cSys->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cSys->cd(i+1);
		purUncHisto[i] = new TH1D(Form("purUncHisto_%d",i),"",10,0,0.5);
		angResHisto[i] = new TH1D(Form("angResHisto_%d",i),"",10,0,0.5);
		jetUncHisto[i] = new TH1D(Form("jetUncHisto_%d",i),"",10,0,0.5);
		jerUncHisto[i] = new TH1D(Form("jerUncHisto_%d",i),"",10,0,0.5);
		subjetTagUncHisto[i] = new TH1D(Form("subjetTagUncHisto_%d",i),"",10,0,0.5);
		jpCalibUncHisto[i] = new TH1D(Form("jpCalibUncHisto_%d",i),"",10,0,0.5);
		dMesonUncHisto[i] = new TH1D(Form("dMesonUncHisto_%d",i),"",10,0,0.5);
		jesNonclosureHisto[i] = new TH1D(Form("jesNonclosureHisto_%d",i),"",10,0,0.5);
		totalUnc[i] = new TH1D(Form("totalUnc_%d",i),"",10,0,0.5);
		
		for(int ibin=1; ibin<=10; ibin++){
			purUncHisto[i]->SetBinContent(ibin, purUnc[ibin-1][i]);
			angResHisto[i]->SetBinContent(ibin, angRes[ibin-1][i]);
			jetUncHisto[i]->SetBinContent(ibin, jesUnc[ibin-1][i]);
			jerUncHisto[i]->SetBinContent(ibin, jerUnc[ibin-1][i]);
			subjetTagUncHisto[i]->SetBinContent(ibin, subjetTagEff[ibin-1][i]);
			jpCalibUncHisto[i]->SetBinContent(ibin, jpCalibSys[ibin-1][i]);
			dMesonUncHisto[i]->SetBinContent(ibin, DMesonSys[ibin-1][i]);
			jesNonclosureHisto[i]->SetBinContent(ibin, abs(jesNonclosureB[ibin-1][i]));
			totalUnc[i]->SetBinContent(ibin, sqrt(pow(purUnc[ibin-1][i],2)+pow(angRes[ibin-1][i],2)+pow(jesUnc[ibin-1][i],2)+pow(jerUnc[ibin-1][i],2)+pow(subjetTagEff[ibin-1][i],2)+pow(jpCalibSys[ibin-1][i],2)+pow(DMesonSys[ibin-1][i],2)+pow(jesNonclosureB[ibin-1][i],2)));
		}
		
		purUncHisto[i]->SetLineColor(kRed+2);
		angResHisto[i]->SetLineColor(kOrange+2);
		jetUncHisto[i]->SetLineColor(kGreen+2);
		jerUncHisto[i]->SetLineColor(kGreen-1);
		subjetTagUncHisto[i]->SetLineColor(kMagenta+2);
		jpCalibUncHisto[i]->SetLineColor(kCyan+2);
		dMesonUncHisto[i]->SetLineColor(kYellow+3);
		jesNonclosureHisto[i]->SetLineColor(kViolet+6);
		totalUnc[i]->SetLineColor(1);
		totalUnc[i]->SetFillColor(kGray);
		//totalUnc[i]->SetFillStyle(1);
		
		purUncHisto[i]->SetLineWidth(2);
		angResHisto[i]->SetLineWidth(2);
		jetUncHisto[i]->SetLineWidth(2);
		jerUncHisto[i]->SetLineWidth(2);
		jesNonclosureHisto[i]->SetLineWidth(2);
		subjetTagUncHisto[i]->SetLineWidth(2);
		jpCalibUncHisto[i]->SetLineWidth(2);
		dMesonUncHisto[i]->SetLineWidth(2);
		
		totalUnc[i]->SetMaximum(0.40);
		totalUnc[i]->SetMinimum(-0.01);
		totalUnc[i]->SetXTitle("zG");
		totalUnc[i]->SetYTitle("Uncertainty");
		totalUnc[i]->Draw();
		totalUnc[i]->GetXaxis()->SetRangeUser(0.1,0.5);
		purUncHisto[i]->Draw("same");
		angResHisto[i]->Draw("same");
		jetUncHisto[i]->Draw("same");
		jerUncHisto[i]->Draw("same");
		subjetTagUncHisto[i]->Draw("same");
		jesNonclosureHisto[i]->Draw("same");
		jpCalibUncHisto[i]->Draw("same");
		dMesonUncHisto[i]->Draw("same");
		
		if(i==1){
			uncLeg->AddEntry(purUncHisto[i],"Non-b#bar{b} Contamination","l");
			uncLeg->AddEntry(angResHisto[i],"Subjet Angular Resolution","l");
			uncLeg->AddEntry(jetUncHisto[i],"JES Uncertainty","l");
			uncLeg->AddEntry(jesNonclosureHisto[i], "JES MC Nonclosure","l");
			uncLeg->AddEntry(jerUncHisto[i],"JER Uncertainty","l");
			uncLeg->AddEntry(subjetTagUncHisto[i], "Subjet Tagging Efficiency","l");
			uncLeg->AddEntry(jpCalibUncHisto[i], "JP-Tagger Calibration","l");
			uncLeg->AddEntry(dMesonUncHisto[i], "c#rightarrowD#rightarrow#mu BR+Decay Reweighting","l");
			uncLeg->AddEntry(totalUnc[i], "Total","lf");
			uncLeg->SetTextSize(0.035);
			uncLeg->Draw("same");
		}
		
		if(i<ptBins-1) labelsSys[i] = new TLatex(0.15,0.35,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		else labelsSys[i] = new TLatex(0.15,0.35,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		labelsSys[i]->Draw("Same");
		
		/*cout << "ptbin " << i << endl;
		cout << "gsp jets............" << endl;
		cout << "purity: "<< getYMean(purUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "ang res: "<< getYMean(angResHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JES data/mc: "<< getYMean(jetUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JES mc: "<< getYMean(jesNonclosureHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JER: "<< getYMean(jerUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "tag eff: "<< getYMean(subjetTagUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "jp calib: "<< getYMean(jpCalibUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "d-meson BR: "<< getYMean(dMesonUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "TOTAL: "<< getYMean(totalUnc[i], zgSpecBB_data[i]) << endl;
		cout << endl;*/
	}
	
	
	//***************** ZG SYSTEMATIC PLOT (INCLUSIVE JETS) *******************//
	
	TH1D *angResHistoIncl[ptBins], *jetUncHistoIncl[ptBins], *jerUncHistoIncl[ptBins], *jesNonclosureHistoIncl[ptBins], *totalUncIncl[ptBins];
	TLegend *uncLegIncl = new TLegend(0.19,0.64,0.9,0.93);
	TCanvas *cSysIncl = new TCanvas("cSysIncl","",1200,600);
	TLatex *labelsSysIncl[2];
	cSysIncl->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cSysIncl->cd(i+1);
		angResHistoIncl[i] = new TH1D(Form("angResHistoIncl_%d",i),"",10,0,0.5);
		jetUncHistoIncl[i] = new TH1D(Form("jetUncHistoIncl_%d",i),"",10,0,0.5);
		jerUncHistoIncl[i] = new TH1D(Form("jerUncHistoIncl_%d",i),"",10,0,0.5);
		jesNonclosureHistoIncl[i] = new TH1D(Form("jesNonclosureHistoIncl_%d",i),"",10,0,0.5);
		totalUncIncl[i] = new TH1D(Form("totalUncIncl_%d",i),"",10,0,0.5);
		
		for(int ibin=1; ibin<=10; ibin++){
			
			int jerUncInclBin = jerUncIncl[i]->FindBin(zgSpecIncl_data[i]->GetBinCenter(ibin));
			
			angResHistoIncl[i]->SetBinContent(ibin, angRes[ibin-1][i]);
			jetUncHistoIncl[i]->SetBinContent(ibin, jesUncIncl[ibin-1][i]);
			jerUncHistoIncl[i]->SetBinContent(ibin, abs(jerUncIncl[i]->GetBinContent(jerUncInclBin)-1.));
			jesNonclosureHistoIncl[i]->SetBinContent(ibin, abs(jesNonclosureIncl[ibin-1][i]));
			totalUncIncl[i]->SetBinContent(ibin, sqrt(pow(angRes[ibin-1][i],2)+pow(jesUncIncl[ibin-1][i],2)+pow(jerUncHistoIncl[i]->GetBinContent(ibin),2)+pow(jesNonclosureIncl[ibin-1][i],2)));
		}
		
		angResHistoIncl[i]->SetLineColor(kOrange+2);
		jetUncHistoIncl[i]->SetLineColor(kGreen+2);
		jerUncHistoIncl[i]->SetLineColor(kGreen-1);
		jesNonclosureHistoIncl[i]->SetLineColor(kViolet+6);
		totalUncIncl[i]->SetLineColor(1);
		totalUncIncl[i]->SetFillColor(kGray);
		//totalUnc[i]->SetFillStyle(1);
		
		angResHistoIncl[i]->SetLineWidth(2);
		jetUncHistoIncl[i]->SetLineWidth(2);
		jerUncHistoIncl[i]->SetLineWidth(2);
		jesNonclosureHistoIncl[i]->SetLineWidth(2);
				
		totalUncIncl[i]->SetMaximum(0.40);
		totalUncIncl[i]->SetMinimum(-0.01);
		totalUncIncl[i]->SetXTitle("zG");
		totalUncIncl[i]->SetYTitle("Uncertainty");
		totalUncIncl[i]->Draw();
		totalUncIncl[i]->GetXaxis()->SetRangeUser(0.1,0.5);
		angResHistoIncl[i]->Draw("same");
		jetUncHistoIncl[i]->Draw("same");
		jerUncHistoIncl[i]->Draw("same");
		jesNonclosureHistoIncl[i]->Draw("same");
		
		if(i==1){
			uncLegIncl->AddEntry(angResHistoIncl[i],"Subjet Angular Resolution","l");
			uncLegIncl->AddEntry(jetUncHistoIncl[i],"JES Uncertainty","l");
			uncLegIncl->AddEntry(jesNonclosureHistoIncl[i], "JES MC Nonclosure","l");
			uncLegIncl->AddEntry(jerUncHistoIncl[i],"JER Uncertainty","l");
			uncLegIncl->AddEntry(totalUncIncl[i], "Total","lf");
			uncLegIncl->SetTextSize(0.035);
			uncLegIncl->Draw("same");
		}
		
		if(i<ptBins-1) labelsSysIncl[i] = new TLatex(0.15,0.35,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		else labelsSysIncl[i] = new TLatex(0.15,0.35,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		labelsSysIncl[i]->Draw("Same");
		
		/*cout << "ptbin " << i << endl;
		cout << "inclusive jets............" << endl;
		cout << "ang res: "<< getYMean(angResHistoIncl[i], zgSpecIncl_data[i]) << endl;
		cout << "JES data/mc: "<< getYMean(jetUncHistoIncl[i], zgSpecIncl_data[i]) << endl;
		cout << "JES mc: "<< getYMean(jesNonclosureHistoIncl[i], zgSpecIncl_data[i]) << endl;
		cout << "JER: "<< getYMean(jerUncHistoIncl[i], zgSpecIncl_data[i]) << endl;
		cout << "TOTAL: "<< getYMean(totalUncIncl[i], zgSpecIncl_data[i]) << endl;
		cout << endl;*/
	}
	
	///************* GSP jet DR SYSTEMATIC PLOT ******************//
	
	TH1D *purUncHistoDR[ptBins], *subjetTagEffHistoDR[ptBins], *jetUncHistoDR[ptBins], *jerUncHistoDR[ptBins], *jesNonclosureHistoDR[ptBins], *totalUncDR[ptBins], *angResHistoDR[ptBins];
	TLegend *uncLeg2 = new TLegend(0.19,0.64,0.9,0.93);
	TLatex *labelsDRSys[2];
	TCanvas *cSys2 = new TCanvas("cSys2","",1200,600);
	cSys2->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cSys2->cd(i+1);
		purUncHistoDR[i] = new TH1D(Form("purUncHistoDR_%d",i),"",20,0,0.4);
		jetUncHistoDR[i] = new TH1D(Form("jetUncHistoDR_%d",i),"",20,0,0.4);
		jerUncHistoDR[i] = new TH1D(Form("jerUncHistoDR_%d",i),"",20,0,0.4);
		angResHistoDR[i] = new TH1D(Form("angResHistoDR_%d",i),"",20,0,0.4);
		subjetTagEffHistoDR[i] = new TH1D(Form("subjetTagEffDR_%d",i),"",20,0,0.4);
		jesNonclosureHistoDR[i] = new TH1D(Form("jesNonclosureHistoDR_%d",i),"",20,0,0.4);
		totalUncDR[i] = new TH1D(Form("totalUncDR_%d",i),"",20,0,0.4);
		
		for(int ibin=1; ibin<=20; ibin++){
			purUncHistoDR[i]->SetBinContent(ibin, purUncDR[ibin-1][i]);
			jetUncHistoDR[i]->SetBinContent(ibin, jesUncDR[ibin-1][i]);
			jerUncHistoDR[i]->SetBinContent(ibin, jerUncDR[ibin-1][i]);
			angResHistoDR[i]->SetBinContent(ibin, angResDR[ibin-1][i]);
			subjetTagEffHistoDR[i]->SetBinContent(ibin, subjetTagEffDR[ibin-1][i]);
			jesNonclosureHistoDR[i]->SetBinContent(ibin, abs(jesNonclosureBDR[ibin-1][i]));
			totalUncDR[i]->SetBinContent(ibin, sqrt(pow(subjetTagEffDR[ibin-1][i],2)+pow(purUncDR[ibin-1][i],2)+pow(jesUncDR[ibin-1][i],2)+pow(jerUncDR[ibin-1][i],2)+pow(jesNonclosureBDR[ibin-1][i],2)+pow(angResDR[ibin-1][i],2)));
		}
		
		purUncHistoDR[i]->SetLineColor(kRed+2);
		jetUncHistoDR[i]->SetLineColor(kGreen+2);
		jerUncHistoDR[i]->SetLineColor(kGreen-1);
		angResHistoDR[i]->SetLineColor(kOrange+2);
		jesNonclosureHistoDR[i]->SetLineColor(kViolet+6);
		subjetTagEffHistoDR[i]->SetLineColor(kMagenta+2);
		totalUncDR[i]->SetLineColor(1);
		totalUncDR[i]->SetFillColor(kGray);
		
		purUncHistoDR[i]->SetLineWidth(2);
		jetUncHistoDR[i]->SetLineWidth(2);
		jerUncHistoDR[i]->SetLineWidth(2);
		angResHistoDR[i]->SetLineWidth(2);
		jesNonclosureHistoDR[i]->SetLineWidth(2);
		subjetTagEffHistoDR[i]->SetLineWidth(2);
				
		totalUncDR[i]->SetMaximum(0.60);
		totalUncDR[i]->SetMinimum(-0.01);
		totalUncDR[i]->SetXTitle("#Delta R");
		totalUncDR[i]->SetYTitle("Uncertainty");
		totalUncDR[i]->Draw();
		totalUncDR[i]->GetXaxis()->SetRangeUser(0.1,0.4);
		purUncHistoDR[i]->Draw("same");
		jetUncHistoDR[i]->Draw("same");
		jerUncHistoDR[i]->Draw("same");
		angResHistoDR[i]->Draw("same");
		jesNonclosureHistoDR[i]->Draw("same");
		subjetTagEffHistoDR[i]->Draw("same");
				
		if(i==1){
			uncLeg2->AddEntry(purUncHistoDR[i],"Non-b#bar{b} Contamination","l");
			uncLeg2->AddEntry(subjetTagEffHistoDR[i],"Subjet Tagging Efficiency","l");
			uncLeg2->AddEntry(jetUncHistoDR[i],"JES Uncertainty","l");
			uncLeg2->AddEntry(jesNonclosureHistoDR[i], "JES MC Nonclosure","l");
			uncLeg2->AddEntry(jerUncHistoDR[i],"JER Uncertainty","l");
			uncLeg2->AddEntry(angResHistoDR[i],"Subjet Angular Resolution", "l");
			uncLeg2->AddEntry(totalUncDR[i], "Total","lf");
			uncLeg2->SetTextSize(0.035);
			uncLeg2->Draw("same");
		}
		
		if(i<ptBins-1) labelsDRSys[i] = new TLatex(0.15,0.35,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		else labelsDRSys[i] = new TLatex(0.15,0.35,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		labelsDRSys[i]->Draw("Same");
		
		/*cout << "ptbin " << i << endl;
		cout << "purity: "<< getYMean(purUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JES data/mc: "<< getYMean(jetUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JES mc: "<< getYMean(jesNonclosureHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JER: "<< getYMean(jerUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "TOTAL: "<< getYMean(totalUnc[i], zgSpecBB_data[i]) << endl;*/
	}
	
	///************* Inclusive Jet DR SYSTEMATIC PLOT ******************//
		
	TH1D *jetUncHistoDRIncl[ptBins], *jerUncHistoDRIncl[ptBins], *jesNonclosureHistoDRIncl[ptBins], *angResHistoDRIncl[ptBins], *totalUncDRIncl[ptBins];
	TLegend *uncLeg2Incl = new TLegend(0.19,0.64,0.9,0.93);
	TLatex *labelsDRSysIncl[2];
	TCanvas *cSys2Incl = new TCanvas("cSys2Incl","",1200,600);
	cSys2Incl->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cSys2Incl->cd(i+1);
		jetUncHistoDRIncl[i] = new TH1D(Form("jetUncHistoDRIncl_%d",i),"",20,0,0.4);
		jerUncHistoDRIncl[i] = new TH1D(Form("jerUncHistoDRIncl_%d",i),"",20,0,0.4);
		angResHistoDRIncl[i] = new TH1D(Form("angResHistoDRIncl_%d",i),"",20,0,0.4);
		jesNonclosureHistoDRIncl[i] = new TH1D(Form("jesNonclosureHistoDRIncl_%d",i),"",20,0,0.4);
		totalUncDRIncl[i] = new TH1D(Form("totalUncDRIncl_%d",i),"",20,0,0.4);
		
		for(int ibin=1; ibin<=20; ibin++){
			if(jetUncHistoDRIncl[i]->GetBinCenter(ibin)<0.1){
				jetUncHistoDRIncl[i]->SetBinContent(ibin,0);
				jerUncHistoDRIncl[i]->SetBinContent(ibin,0);
				jesNonclosureHistoDRIncl[i]->SetBinContent(ibin,0);
			}
			else{
				jetUncHistoDRIncl[i]->SetBinContent(ibin, jesUncInclDR[ibin-1][i]);
				jerUncHistoDRIncl[i]->SetBinContent(ibin, jerUncInclDR[ibin-1][i]);
				angResHistoDRIncl[i]->SetBinContent(ibin, angResDR[ibin-1][i]);
				jesNonclosureHistoDRIncl[i]->SetBinContent(ibin, abs(jesNonclosureInclDR[ibin-1][i]));
				totalUncDRIncl[i]->SetBinContent(ibin, sqrt(pow(jesUncInclDR[ibin-1][i],2)+pow(jerUncInclDR[ibin-1][i],2)+pow(jesNonclosureInclDR[ibin-1][i],2)+pow(angResDR[ibin-1][i],2)));
			}
		}
		
		jetUncHistoDRIncl[i]->SetLineColor(kGreen+2);
		jerUncHistoDRIncl[i]->SetLineColor(kGreen-1);
		jesNonclosureHistoDRIncl[i]->SetLineColor(kViolet+6);
		angResHistoDRIncl[i]->SetLineColor(kOrange+2);
		totalUncDRIncl[i]->SetLineColor(1);
		totalUncDRIncl[i]->SetFillColor(kGray);
		
		jetUncHistoDRIncl[i]->SetLineWidth(2);
		jerUncHistoDRIncl[i]->SetLineWidth(2);
		jesNonclosureHistoDRIncl[i]->SetLineWidth(2);
		angResHistoDRIncl[i]->SetLineWidth(2);
				
		totalUncDRIncl[i]->SetMaximum(0.60);
		totalUncDRIncl[i]->SetMinimum(-0.01);
		totalUncDRIncl[i]->SetXTitle("#Delta R");
		totalUncDRIncl[i]->SetYTitle("Uncertainty");
		totalUncDRIncl[i]->Draw();
		totalUncDRIncl[i]->GetXaxis()->SetRangeUser(0.1,0.4);
		jetUncHistoDRIncl[i]->Draw("same");
		jerUncHistoDRIncl[i]->Draw("same");
		jesNonclosureHistoDRIncl[i]->Draw("same");
		angResHistoDRIncl[i]->Draw("same");
		
		if(i==1){
			uncLeg2Incl->AddEntry(jetUncHistoDRIncl[i],"JES Uncertainty","l");
			uncLeg2Incl->AddEntry(jesNonclosureHistoDRIncl[i], "JES MC Nonclosure","l");
			uncLeg2Incl->AddEntry(jerUncHistoDRIncl[i],"JER Uncertainty","l");
			uncLeg2Incl->AddEntry(angResHistoDR[i],"Subjet Angular Resolution", "l");
			uncLeg2Incl->AddEntry(totalUncDRIncl[i], "Total","lf");
			uncLeg2Incl->SetTextSize(0.035);
			uncLeg2Incl->Draw("same");
		}
		
		if(i<ptBins-1) labelsDRSysIncl[i] = new TLatex(0.15,0.35,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		else labelsDRSysIncl[i] = new TLatex(0.15,0.35,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		labelsDRSysIncl[i]->Draw("Same");
		
		/*cout << "ptbin " << i << endl;
		cout << "purity: "<< getYMean(purUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JES data/mc: "<< getYMean(jetUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JES mc: "<< getYMean(jesNonclosureHisto[i], zgSpecBB_data[i]) << endl;
		cout << "JER: "<< getYMean(jerUncHisto[i], zgSpecBB_data[i]) << endl;
		cout << "TOTAL: "<< getYMean(totalUnc[i], zgSpecBB_data[i]) << endl;*/
	}
	
	//****************** MC NONCLOSURE CHECK ********************//
	
	TLatex *labelsSmear[ptBins];
	TCanvas *cSmear = new TCanvas("cSmear","",1200,600);
	cSmear->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cSmear->cd(i+1);
		zgSpecBB_data_shift[i]->SetLineColor(kMagenta+2);
		zgSpecBB_data_shift[i]->SetYTitle("Relative z_{G} shape");
		zgSpecBB_data_shift[i]->SetXTitle("z_{G}");
		zgSpecBB_data_shift[i]->SetMarkerColor(kMagenta+2);
		zgSpecBB_data_shift[i]->Draw();
		zgSpecBB_data_shift[i]->GetYaxis()->SetRangeUser(0.5,1.5);
		zgSpecIncl_data_shift[i]->Draw("same");
		if(i<ptBins-1) labelsSmear[i] = new TLatex(0.05,1.3,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		else labelsSmear[i] = new TLatex(0.05,1.3,Form("Jet p_{T} > %d GeV",jetPtBin[i]));
		labelsSmear[i]->Draw("Same");
	}
	
	TLatex *labelsPur[ptBins];
	TCanvas *cPur = new TCanvas("cPur","",1200,600);
	cPur->Divide(ptBins,1);
	for(int i=0; i<ptBins; i++){
		cPur->cd(i+1);
		InclSpecPurSys[i]->SetYTitle("Purity Contamination");
		InclSpecPurSys[i]->SetXTitle("z_{g}");
		InclSpecPurSys[i]->Draw();
		InclSpecPurSys[i]->GetYaxis()->SetRangeUser(0.5,1.5);
		labelsPur[i] = new TLatex(0.05,1.3,Form("%d < Jet p_{T} < %d GeV",jetPtBin[i],jetPtBin[i+1]));
		labelsPur[i]->Draw("same");
	}
	
	if(drawFromTree){
		/*TCanvas *cc3 = new TCanvas("cc3","",1200,800);
		cc3->Divide(ptBins/2,ptBins/2);
		TH1D *py[4], *datapy[4];
		for(int i=0; i<ptBins; i++){
			cc3->cd(i+1);
			py[i] = zGdRCorr_mc[i]->ProjectionY();
			py[i]->SetMarkerColor(2);
			py[i]->SetLineColor(2);
			py[i]->Scale(1./py[i]->Integral());
			py[i]->Draw();
			datapy[i] = zGdRcorr_data[i]->ProjectionY();
			datapy[i]->Scale(1./datapy[i]->Integral());
			datapy[i]->Draw("Same");
		}
		
		TCanvas *cc4 = new TCanvas("cc4","",1200,800);
		cc4->Divide(ptBins/2,ptBins/2);
		TH1D *py2[4], *datapy2[4];
		for(int i=0; i<ptBins; i++){
			cc4->cd(i+1);
			py2[i] = zGdRCorr_bb_mc[i]->ProjectionY();
			py2[i]->SetMarkerColor(2);
			py2[i]->SetLineColor(2);
			py2[i]->Scale(1./py2[i]->Integral());
			py2[i]->Draw();
			datapy2[i] = zGdRcorr_csv_data[i]->ProjectionY();
			datapy2[i]->Scale(1./datapy2[i]->Integral());
			datapy2[i]->Draw("Same");
		}*/
	}
	
	if(writeOut){
		TFile *fout = new TFile("zgPlots.root","recreate");
		fout->cd();
		for(int i=0; i<ptBins; i++){
			zgSpecIncl_data[i]->Write();
			zgSpecBB_data[i]->Write();
			zgSpecTag_data[i]->Write();
			zgSpecDR_data[i]->Write();
			
			zgSpecIncl_mc[i]->Write();
			zgSpecBB_mc[i]->Write();
			
			zGdRCorr_mc[i]->Write();
			zGdRCorr_bb_mc[i]->Write();
			zGdRcorr_data[i]->Write();
			zGdRcorr_csv_data[i]->Write();
			
		}
		fout->Close();
	}
	
	
}