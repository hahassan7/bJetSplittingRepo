

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"

void compareDataMC(){
	
	const bool ispp = 1;

	const int nVars=9;
	string vars[nVars] = {"jtpt", "jteta", "jtSubJetPt1", "jtSubJetPt2", "minCSV", "minJP","discr_csvV2","ndiscr_csvV2","discr_prob"};
	string cuts[nVars] = {"jtpt>80 && jtpt<200 && abs(jteta)<2", "jtpt>80 && jtpt<200 && abs(jteta)<2", "jtSubJetPt1>60 && jtSubJetPt1<200 && abs(jtSubJetEta1)<2", "jtSubJetPt1>60 && jtSubJetPt1<200 && abs(jtSubJetEta1)<2", "jtSubJetPt1>60 && jtSubJetPt1<200 && jtSubJetPt2>60 && jtSubJetPt2<200 && abs(jtSubJetEta1)<2 && abs(jtSubJetEta2)<2", "jtSubJetPt1>60 && jtSubJetPt1<200 && jtSubJetPt2>60 && jtSubJetPt2<200 && abs(jtSubJetEta1)<2 && abs(jtSubJetEta2)<2", "jtpt>80 && jtpt<200 && abs(jteta)<2", "jtpt>80 && jtpt<200 && abs(jteta)<2","jtpt>80 && jtpt<200 && abs(jteta)<2"};
	string flavorSelector[nVars] = {"jtHadronFlavor","jtHadronFlavor","jtSubJetHadronFlavor1","jtSubJetHadronFlavor2","jtSubJetHadronFlavor1","jtSubJetHadronFlavor1", "jtHadronFlavor", "jtHadronFlavor","jtHadronFlavor"};
	string evtCuts = "";
	
	string titles[nVars] = {"Jet pT", "Jet Eta", "Leading Subjet pT", "Subleading Subjet pT", "minCSV", "minJP", "CSV Discr.", "Neg CSV Discr.","Jet Prob"};
	
	const int nBins[nVars] = {40, 40, 40, 40, 20, 20, 20, 20, 50};
	const double lo[nVars] = {80, -2, 50, 50, 0, 0, 0, 0, 0};
	const double hi[nVars] = {200, 2, 200, 200, 1, 2.5, 1, 1, 3};
	
	TH1D *varHistos[nVars];
	TH1D *varHistosLight[nVars];
	TH1D *varHistosCharm[nVars];
	TH1D *varHistosBottom[nVars];
	TH1D *hRatio[nVars];
	THStack *mcstack[nVars];
	TLegend *l1[nVars];
	
	TCanvas *cc[nVars];
	
	TFile *fdata = new TFile("gspTreeOut_ppData_explicitJTA_recalibJP_pt60.root");
	TFile *fmc = new TFile("gspTreeOut_ppMC_Pythia8_explicitJTA_recalibJP_withNegCSV.root");
	
	TTree *tdata = (TTree*)fdata->Get("gspTree");
	TTree *tmc = (TTree*)fmc->Get("gspTree");
	
	for(int i=0; i<nVars; i++){
		cc[i] = new TCanvas(Form("cc_%d",i),"",600,600);
		cc[i]->cd();
		TPad *p1 = new TPad(Form("p1_%d",i),"",0,0.2,1.0,1.0);
		p1->SetBottomMargin(0.01);
		TPad *p2 = new TPad(Form("p2_%d",i),"",0,0,1.0,0.2);
		p2->SetTopMargin(0.01);
		p2->SetBottomMargin(0.4);
		p1->Draw();
		p2->Draw();
		string hName = vars[i]+"_histos";
		varHistos[i] = new TH1D(hName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistos[i]->Sumw2();
		string lightName = vars[i]+"_light";
		varHistosLight[i] = new TH1D(lightName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistosLight[i]->Sumw2();
		varHistosLight[i]->SetFillColor(kBlue-3);
		string charmName = vars[i]+"_charm";
		varHistosCharm[i] = new TH1D(charmName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistosCharm[i]->Sumw2();
		varHistosCharm[i]->SetFillColor(kGreen-3);
		string botName = vars[i]+"_bottom";
		varHistosBottom[i] = new TH1D(botName.c_str(),"",nBins[i],lo[i],hi[i]);
		varHistosBottom[i]->Sumw2();
		varHistosBottom[i]->SetFillColor(kRed-3);
		tdata->Project(hName.c_str(),vars[i].c_str(),Form("weight*(%s)", cuts[i].c_str()));
		tmc->Project(lightName.c_str(),vars[i].c_str(),Form("weight*(%s && abs(%s)!=4 && abs(%s)!=5)", cuts[i].c_str(), flavorSelector[i].c_str(), flavorSelector[i].c_str()));
		tmc->Project(charmName.c_str(),vars[i].c_str(),Form("weight*(%s && abs(%s)==4)", cuts[i].c_str(), flavorSelector[i].c_str()));
		tmc->Project(botName.c_str(),vars[i].c_str(),Form("weight*(%s && abs(%s)==5)", cuts[i].c_str(), flavorSelector[i].c_str()));
		
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

	}
	
}
