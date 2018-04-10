

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"

void compareDataMC(){
	
	bool saveHistos = false;
	bool doMuFilter = true;
	bool doPowheg = false;
	
	const bool ispp = 1;

	const int nVars=14;
	string vars[nVars] = {"jtpt", "jteta", "jtSubJetPt1", "jtSubJetPt2", "jtSubJetEta1","jtSubJetEta2","minCSV", "minJP","minJP","discr_csvV2","discr_prob","discr_prob","jtSubJetPt1","jtSubJetPt2"};
	string cuts[nVars] = {"matchedJtpt>120 && matchedJtpt<400 && abs(jteta)<2 && jtSubJetdR>0.1", 
	"matchedJtpt>120 && matchedJtpt<400 && abs(jteta)<2 && jtSubJetdR>0.1", 
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta1)<2 && jtSubJetdR>0.1", 
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta2)<2 && jtSubJetdR>0.1",
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta1)<2 && jtSubJetdR>0.1", 
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta2)<2 && jtSubJetdR>0.1",
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta1)<2 && abs(jtSubJetEta2)<2 && jtSubJetdR>0.1", 
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta1)<2 && abs(jtSubJetEta2)<2 && jtSubJetdR>0.1",
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta1)<2 && abs(jtSubJetEta2)<2 && jtSubJetdR>0.1 && jtSubJetSvtxm1>1.2 && jtSubJetSvtxm2>1.2 && minCSV>0.8",
	"matchedJtpt>120 && matchedJtpt<400 && abs(jteta)<2 && jtSubJetdR>0.1", 
	"matchedJtpt>120 && matchedJtpt<400 && abs(jteta)<2 && jtSubJetdR>0.1",
	"matchedJtpt>120 && matchedJtpt<400 && abs(jteta)<2 && jtSubJetdR>0.1 && jtSubJetSvtxm1>1.2 && jtSubJetSvtxm2>1.2 && minCSV>0.8",
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta1)<2.0 && jtSubJetcsvV2_1>0.8 && jtSubJetSvtxm1>1.2 && jtSubJetdR>0.1",
	"matchedJtpt>120 && matchedJtpt<400 && abs(jtSubJetEta2)<2.0 && jtSubJetcsvV2_2>0.8 && jtSubJetSvtxm2>1.2 && jtSubJetdR>0.1",};
	
	if(doMuFilter){
		for(int i=0; i<nVars; i++){
			cuts[i]+=" && muPt > 5";
		}
	}
	
	string flavorSelector[nVars] = {"jtHadronFlavor","jtHadronFlavor","jtSubJetHadronFlavor1","jtSubJetHadronFlavor2","jtSubJetHadronFlavor1","jtSubJetHadronFlavor2","jtSubJetHadronFlavor1","jtSubJetHadronFlavor1", "jtSubJetHadronFlavor1", "jtHadronFlavor", "jtHadronFlavor","jtHadronFlavor","jtSubJetHadronFlavor1","jtSubJetHadronFlavor2"};
	string evtCuts = "";
	
	string titles[nVars] = {"Jet pT", "Jet Eta","Leading Subjet pT", "Subleading Subjet pT", "Leading Subjet Eta", "Subleading Subjet Eta", "minCSV", "minJP", "minJP (Tagged)", "CSV Discr", "Jet Prob","Jet Prob (Tagged)","Tagged Leading Subjet pT", "Tagged Subleading Subjet pT"};
	
	bool isLogy[nVars] = {1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1};
	
	const int nBins[nVars] = {40, 40, 40, 40, 40, 40, 20, 20, 20, 20, 50, 50, 40, 40};
	const double lo[nVars] = {80, -2, 50, 50, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0};
	const double hi[nVars] = {200, 2, 200, 200, 2, 2, 1, 2.5, 2.5, 1, 3, 3, 150, 150};
	
	TH1D *varHistos[nVars];
	TH1D *varHistosLight[nVars];
	TH1D *varHistosCharm[nVars];
	TH1D *varHistosBottom[nVars];
	TH1D *hRatio[nVars];
	THStack *mcstack[nVars];
	TLegend *l1[nVars];
	
	TCanvas *cc[nVars];
	
	//WARNING!! If you don't have a muPt > 5 GeV cut, your JP distributions between data and MC will NOT match!!
	
	TFile *fdata, *fmc, *fmcLight, *fbmc;
	
	if(doMuFilter){
		fdata = new TFile("ntuples/gspTreeOut_ppData_muFilterPerJet_muPt4.root");
		fmc = new TFile("ntuples/gspTreeOut_ppMC_qcdPythia8_fullSample_withFullJets.root");
		fmcLight = new TFile("ntuples/gspTreeOut_ppMC_qcdPythia8_fullSample_withFullJets.root");
		fbmc = new TFile("ntuples/gspTreeOut_pt50BHadron_muPt4Filtered.root");
	}
	else if(doPowheg){
		fdata = new TFile("ntuples/gspTreeOut_ppData_trgFilter_lowPtCuts.root");
		fmc = new TFile("ntuples/gspTreeOut_Powheg.root");
		fmcLight = new TFile("ntuples/gspTreeOut_Powheg.root");
		fbmc = new TFile("ntuples/gspTreeOut_Powheg.root");
	}
	else{
		fdata = new TFile("ntuples/gspTreeOut_ppData_trgFilter_lowPtCuts.root");
		fmc = new TFile("ntuples/gspTreeOut_ppMC_QCDjets_noMuFilter.root");
		fmcLight = new TFile("ntuples/gspTreeOut_ppMC_QCDjets_noMuFilter.root");
		fbmc = new TFile("ntuples/gspTreeOut_ppMC_QCDjets_noMuFilter.root");
	}
	
	TTree *tdata = (TTree*)fdata->Get("gspTree");
	TTree *tmc = (TTree*)fmc->Get("gspTree");
	TTree *tmcLight = (TTree*)fmcLight->Get("gspTree");
	TTree *tbmc = (TTree*)fbmc->Get("gspTree");
	
	for(int i=0; i<nVars; i++){
		cc[i] = new TCanvas(Form("cc_%d",i),"",600,600);
		cc[i]->cd();
		TPad *p1 = new TPad(Form("p1_%d",i),"",0,0.2,1.0,1.0);
		p1->SetBottomMargin(0.01);
		TPad *p2 = new TPad(Form("p2_%d",i),"",0,0,1.0,0.2);
		p2->SetTopMargin(0.01);
		p2->SetBottomMargin(0.4);
		if(isLogy[i]) p1->SetLogy();
		p1->Draw();
		p2->Draw();
		string hName = vars[i]+Form("_histos_%d",i);
		varHistos[i] = new TH1D(hName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistos[i]->Sumw2();
		string lightName = vars[i]+Form("_light_%d",i);
		varHistosLight[i] = new TH1D(lightName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistosLight[i]->Sumw2();
		varHistosLight[i]->SetFillColor(kBlue-3);
		string charmName = vars[i]+Form("_charm_%d",i);
		varHistosCharm[i] = new TH1D(charmName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistosCharm[i]->Sumw2();
		varHistosCharm[i]->SetFillColor(kGreen-3);
		string botName = vars[i]+Form("_bottom_%d",i);
		varHistosBottom[i] = new TH1D(botName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistosBottom[i]->Sumw2();
		varHistosBottom[i]->SetFillColor(kRed-3);
		tdata->Project(hName.c_str(),vars[i].c_str(),Form("weight*(%s)", cuts[i].c_str()));
		tmcLight->Project(lightName.c_str(),vars[i].c_str(),Form("weight*(%s && abs(%s)!=4 && abs(%s)!=5)", cuts[i].c_str(), flavorSelector[i].c_str(), flavorSelector[i].c_str()));
		tmc->Project(charmName.c_str(),vars[i].c_str(),Form("weight*(%s && abs(%s)==4)", cuts[i].c_str(), flavorSelector[i].c_str()));
		tbmc->Project(botName.c_str(),vars[i].c_str(),Form("weight*2.5*(%s && abs(%s)==5)", cuts[i].c_str(), flavorSelector[i].c_str()));
		
		TH1D *htemp = (TH1D*)varHistosLight[i]->Clone("htemp");
		htemp->Add(varHistosCharm[i]);
		htemp->Add(varHistosBottom[i]);
		varHistosLight[i]->Scale(varHistos[i]->Integral(10,20)/htemp->Integral(10,20));
		varHistosCharm[i]->Scale(varHistos[i]->Integral(10,20)/htemp->Integral(10,20));
		varHistosBottom[i]->Scale(varHistos[i]->Integral(10,20)/htemp->Integral(10,20));
		TH1D *htemp2 = (TH1D*)varHistosLight[i]->Clone("htemp2");
		htemp2->Add(varHistosCharm[i]);
		htemp2->Add(varHistosBottom[i]);

		//varHistos[i]->SetXTitle(vars[i].c_str());
		varHistos[i]->SetYTitle("Counts");
		varHistos[i]->SetTitleSize(0,"X");
		varHistos[i]->SetLabelSize(0,"X");

		mcstack[i] = new THStack(Form("mcstack_%d",i),"");
		mcstack[i]->Add(varHistosBottom[i]);
		mcstack[i]->Add(varHistosCharm[i]);
		mcstack[i]->Add(varHistosLight[i]);

		p1->cd();
		varHistos[i]->SetMinimum(1);
		varHistos[i]->Draw();
		mcstack[i]->Draw("same,hist");
		varHistos[i]->Draw("same");
		p1->RedrawAxis();

		l1[i] = new TLegend(0.5,0.6,0.9,0.8);
		if(!ispp) l1[i]->AddEntry(varHistos[i],"pPb data");
		else l1[i]->AddEntry(varHistos[i],"pp data");
		l1[i]->AddEntry(varHistosLight[i],"usdg jets","f");
		l1[i]->AddEntry(varHistosCharm[i],"charm jets","f");
		l1[i]->AddEntry(varHistosBottom[i],"bottom jets","f");
		l1[i]->Draw("same");

		p2->cd();
		hRatio[i] = (TH1D*)varHistos[i]->Clone(Form("hRatio_%d",i));
		hRatio[i]->Divide(htemp2);
		hRatio[i]->SetXTitle(titles[i].c_str());
		hRatio[i]->SetYTitle("Data/MC");
		hRatio[i]->SetMaximum(1.7);
		hRatio[i]->SetMinimum(0.3);
		hRatio[i]->SetLabelSize(0.2,"X");
		hRatio[i]->SetLabelSize(0.2,"Y");
		hRatio[i]->SetNdivisions(505,"Y");
		hRatio[i]->SetTitleSize(0.2,"X");
		hRatio[i]->SetTitleSize(0.2,"Y");
		hRatio[i]->SetTitleOffset(0.3,"Y");
		hRatio[i]->SetTitleOffset(0.85,"X");
		hRatio[i]->SetTickLength(0.1,"X");
		hRatio[i]->Draw();
		
		string saveTitle;
		for(unsigned int ipos=0; ipos<(unsigned int)titles[i].size(); ipos++){
			if(titles[i].at(ipos) != ' ') saveTitle.push_back(titles[i].at(ipos));
		}
		if(saveHistos){
			if(doMuFilter) cc[i]->SaveAs(Form("dataMCComp_withMuFilter/%s.pdf",saveTitle.c_str()));
			else cc[i]->SaveAs(Form("dataMCComp_noMuFilter/%s.pdf",saveTitle.c_str()));
		} 

	}
	
}