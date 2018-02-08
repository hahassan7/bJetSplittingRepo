
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


using std::endl;
using std::cout;
using std::vector;

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<vector<float> > >+;
#pragma link C++ class vector<vector<int> >+;
#endif

const int MAXJETS=200;
const double PI=3.1415926;

double getDR(double eta1, double phi1, double eta2, double phi2){
	//project phi's to 0->pi
	phi1 = TMath::ACos(TMath::Cos(phi1));
	phi2 = TMath::ACos(TMath::Cos(phi2));
	return(sqrt(pow(eta1-eta2,2)+pow(phi1-phi2,2)));
}

double getzg(float j1, float j2){
	return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}

int getFlavorCategory(int flav1, int flav2){
	if(abs(flav1)==5 && abs(flav2)==5) return 0;
	else if((abs(flav1)==4 && abs(flav2)==5) || (abs(flav1)==5 && abs(flav2)==4)) return 1;
	else if(abs(flav1)==4 && abs(flav2)==4) return 2;
	else if(abs(flav1)==5 || abs(flav2)==5) return 3;
	else if(abs(flav1)==4 || abs(flav2)==4) return 4;
	else return 5;
}

void formatHistos(TH1D **histos){
	TH1D *hAll = (TH1D*)histos[0]->Clone("hAll");
	for(int i=1; i<6; i++){
		hAll->Add(histos[i]);
	}
	for(int i=0; i<6; i++){
		histos[i]->Divide(hAll);
	}
	histos[0]->SetFillColor(kOrange+2); histos[0]->SetTitle("Double-b subjets");
	histos[1]->SetFillColor(kCyan+1); histos[1]->SetTitle("b+c subjets");
	histos[2]->SetFillColor(kGreen+2); histos[2]->SetTitle("Double-c subjets");
	histos[3]->SetFillColor(2); histos[3]->SetTitle("b+light subjets");
	histos[4]->SetFillColor(kMagenta+2); histos[4]->SetTitle("c+light subjets");
	histos[5]->SetFillColor(4); histos[5]->SetTitle("light subjets");
}

void scanGSPevents(){
	
	const bool ispp = true;
	const bool isMC = true;
	
	const int nPt = 4;
	double ptBins[nPt+1] = {100, 120, 160, 250, 9999};
	
	const int nCent = ispp ? 1 : 3;
	int *centBins;
	
	int PbPbCentBins[4] = {0, 20, 60, 200};
	int ppCentBins[2] = {0, 200};
	ispp ? centBins = &ppCentBins[0] : centBins = &PbPbCentBins[0];
	
	TH1D *zgSpec_ppIncl[nPt];
	TH1D *zgSpec_ppGSP[nPt];
	for(int i=0; i<nPt; i++){
		zgSpec_ppIncl[i] = new TH1D(Form("zgSpec_ppIncl_%d",i),"",10,0,0.5);
		zgSpec_ppGSP[i] = new TH1D(Form("zgSpec_ppGSP_%d",i),"",10,0,0.5);
	}
	TH1D *flavorFrac[6]; //bb, bc, cc, bl, cl, ll
	TH1D *flavorFrac_withCuts[6];
	for(int i=0; i<6; i++){
		flavorFrac[i] = new TH1D(Form("flavorFrac_%d",i),"",20,0,1);
		flavorFrac_withCuts[i] = new TH1D(Form("flavorFrac_withCuts_%d",i),"",20,0,1);
	}
	TH1D *allDiBjetsCut[4][nCent], *allJetsCut[4][nCent];
	TH1D *allDiBjets[nCent];
	for(int j=0; j<nCent; j++){
		allDiBjets[j] = new TH1D(Form("allDiBjets_%d",j),"",40,0,1); allDiBjets[j]->Sumw2();
		for(int i=0; i<4; i++){
			allDiBjetsCut[i][j] = new TH1D(Form("allDiBjetsCut_%d_%d",i,j),"",40,0,1); allDiBjetsCut[i][j]->Sumw2();
			allJetsCut[i][j] = new TH1D(Form("allJetsCut_%d_%d",i,j),"",40,0,1); allJetsCut[i][j]->Sumw2();
		}
	}
	
	TH1D *allDiBjets_dR = new TH1D("allDiBjets_dR","",10,0,0.4); allDiBjets_dR->Sumw2();
	TH1D *allDBjetsCut_dR = new TH1D("dRTagEff_csv09","",10,0,0.4); allDBjetsCut_dR->Sumw2();
	
	TH1D *allDiBjets_zg = new TH1D("allDiBjets_zg","",10,0,0.5); allDiBjets_zg->Sumw2();
	TH1D *allDBjetsCut_zg = new TH1D("zgTagEff_csv09","",10,0,0.5); allDBjetsCut_zg->Sumw2();
	
	TH2D *allDiBjets_corr = new TH2D("allDiBjets_corr","",10,0,0.5,10,0,0.4); allDiBjets_corr->Sumw2();
	TH2D *allDiBjetsCut_corr = new TH2D("CorrTagEff_csv09","",10,0,0.5,10,0,0.4); allDiBjetsCut_corr->Sumw2();
	
	//****************************************************************//
	
	int nref;
	float jtpt[MAXJETS], jteta[MAXJETS], jtphi[MAXJETS], refpt[MAXJETS];
	int refparton_flavorForB[MAXJETS], nsvtx[MAXJETS];
	float discr_csvV2[MAXJETS], discr_prob[MAXJETS];
	vector<vector<float> > *jtSubJetPt=0, *jtSubJetEta=0, *jtSubJetPhi=0;
	vector<vector<float> > *jtSubJetPartonFlavor=0, *jtSubJetHadronFlavor=0;
	vector<vector<float> > *jtSubJetCSV=0, *jtSubJetJP=0;
	int refparton_flavorProcess[MAXJETS];
	float jtHadronFlavor[MAXJETS];
	float jtm[MAXJETS];
	vector<vector<float> >*svtxm=0, *svtxpt=0, *svtxdl=0, *svtxdls=0;
	vector<vector<int> >*svtxntrk=0, *jtSubJetVtxType=0;
	
	vector<vector<vector<float> > > *jtSubJetHadronDR=0, *jtSubJetPartonDR=0, *jtSubJetPartonPt=0, *jtSubJetPartonPdg=0;
	vector<vector<vector<float> > > *jtSubJetHadronPt=0, *jtSubJetHadronEta=0, *jtSubJetHadronPhi=0, *jtSubJetHadronPdg=0;
	
	vector<vector<vector<float> > > *jtSubJetSvtxm=0, *jtSubJetSvtxpt=0, *jtSubJetSvtxeta=0, *jtSubJetSvtxphi=0;
	vector<vector<vector<float> > > *jtSubJetSvtxNtrk=0, *jtSubJetSvtxdl=0, *jtSubJetSvtxdls=0;
	
	///*************************************************************//
	
	double t_jtpt, t_jteta, t_jtphi, t_discr_csvV2, t_jtSubJetPt1, t_jtSubJetPt2, t_jtSubJetEta1, t_jtSubJetEta2, t_refpt;
	double t_jtSubJetPhi1, t_jtSubJetPhi2, t_jtSubJetCSV1, t_jtSubJetCSV2, t_jtSubJetdR, t_discr_prob;
	double t_jtSubJetVtxType1, t_jtSubJetVtxType2, t_minCSV, t_jtSubJetJP1, t_jtSubJetJP2, t_minJP;
        double t_jtSubJetHadronDR1, t_jtSubJetHadronDR2;
	int t_jtHadronFlavor, t_jtSubJetHadronFlavor1, t_jtSubJetHadronFlavor2;
	int hiBin;
	
	TFile *fout = new TFile("gspTreeOut.root","recreate");
	TTree *gspTree = new TTree("gspTree","");
	gspTree->Branch("jtpt",&t_jtpt);
	gspTree->Branch("jteta",&t_jteta);
	gspTree->Branch("jtphi",&t_jtphi);
	gspTree->Branch("discr_csvV2",&t_discr_csvV2);
	gspTree->Branch("discr_prob", &t_discr_prob);
	gspTree->Branch("jtSubJetPt1",&t_jtSubJetPt1);
	gspTree->Branch("jtSubJetPt2",&t_jtSubJetPt2);
	gspTree->Branch("jtSubJetEta1", &t_jtSubJetEta1);
	gspTree->Branch("jtSubJetPhi1", &t_jtSubJetPhi1);
	gspTree->Branch("jtSubJetEta2", &t_jtSubJetEta2);
	gspTree->Branch("jtSubJetPhi2", &t_jtSubJetPhi2);
	gspTree->Branch("jtSubJetcsvV1_1",&t_jtSubJetCSV1);
	gspTree->Branch("jtSubJetcsvV1_2",&t_jtSubJetCSV2);
	gspTree->Branch("jtSubJetJP_1",&t_jtSubJetJP1);
	gspTree->Branch("jtSubJetJP_2",&t_jtSubJetJP2);
	gspTree->Branch("minCSV",&t_minCSV);
	gspTree->Branch("minJP",&t_minJP);
	gspTree->Branch("jtSubJetdR",&t_jtSubJetdR);
	gspTree->Branch("jtSubJetVtxType1",&t_jtSubJetVtxType1);
	gspTree->Branch("jtSubJetVtxType2",&t_jtSubJetVtxType2);
	if(isMC){
		gspTree->Branch("refpt",&t_refpt);
		gspTree->Branch("jtHadronFlavor",&t_jtHadronFlavor);
		gspTree->Branch("jtSubJetHadronFlavor1",&t_jtSubJetHadronFlavor1);
		gspTree->Branch("jtSubJetHadronFlavor2",&t_jtSubJetHadronFlavor2);
                gspTree->Branch("jtSubJetHadronDR1",&t_jtSubJetHadronDR1);
                gspTree->Branch("jtSubJetHadronDR2",&t_jtSubJetHadronDR2);
        }
	if(!ispp){
		gspTree->Branch("hiBin",&hiBin);
	}
	
	
	std::string filelist;
	if(ispp && isMC) filelist = "ppMC_allPthat_withJPtags.txt";
	else if(ispp) filelist = "ppData_Jet80Trg.txt";
      else if(isMC) filelist = "PbPbMC_qcdJetList_pthat120.txt";
      else filelist = "PbPbData_jet80TrgList_part2.txt";
      //else{ cout << "Not doing PbPb data yet! " << endl; exit(0); }
      //std::string filename = "HiForestAOD.root";
	std::string filename;
	ifstream instr(filelist.c_str());
	TFile *fin;
	int ifile=0;
	//for(int ifile=0; ifile<1; ifile++){
	while(instr>>filename && ifile <5000){
		
		ifile++;
		fin = TFile::Open(filename.c_str());
		cout << "starting file " << filename << "... " << endl;
		
		if(!fin){ cout << "WARNING! Can't find file" << endl; continue; }
		
		int pvFilter, beamScrapeFilter, noiseFilter;
		int collEvtSel, hfCoincFilter;
		float pthat;
		TTree *skimTree = (TTree*)fin->Get("skimanalysis/HltTree");
		if(ispp){
			skimTree->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter);
			skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapeFilter);
		}
		else{
			skimTree->SetBranchAddress("pcollisionEventSelection", &collEvtSel);
			skimTree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
			skimTree->SetBranchAddress("phfCoincFilter3",&hfCoincFilter);
			
		}
		skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&noiseFilter);
		
		TTree *evtTree = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
		evtTree->SetBranchAddress("hiBin",&hiBin);
		
		TTree *subjetTree;
		if(!ispp) subjetTree = (TTree*)fin->Get("akCsSoftDrop4PFJetAnalyzer/t");
		else if(isMC) subjetTree = (TTree*)fin->Get("akSoftDrop4PFJetAnalyzer/t");
		else subjetTree = (TTree*)fin->Get("akSoftDrop4PFJetAnalyzer/t");
		
		subjetTree->SetBranchAddress("nref", &nref);
		subjetTree->SetBranchAddress("nsvtx", nsvtx);
		subjetTree->SetBranchAddress("jtpt", jtpt);
		subjetTree->SetBranchAddress("jteta", jteta);
		subjetTree->SetBranchAddress("jtphi", jtphi);
		if(isMC){
			subjetTree->SetBranchAddress("refpt",refpt);
			subjetTree->SetBranchAddress("pthat",&pthat);
			subjetTree->SetBranchAddress("jtHadronFlavor",jtHadronFlavor);
			subjetTree->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
			//subjetTree->SetBranchAddress("refparton_flavorProcess",refparton_flavorProcess);
		}
		subjetTree->SetBranchAddress("discr_csvV2",discr_csvV2);
		subjetTree->SetBranchAddress("discr_prob",discr_prob);
		
		subjetTree->SetBranchAddress("jtSubJetPt", &jtSubJetPt);
		subjetTree->SetBranchAddress("jtSubJetEta", &jtSubJetEta);
		subjetTree->SetBranchAddress("jtSubJetPhi", &jtSubJetPhi);
		subjetTree->SetBranchAddress("jtSubJetcsvV1",&jtSubJetCSV);
		subjetTree->SetBranchAddress("jtSubJetJP",&jtSubJetJP);
		if(isMC){ 
			subjetTree->SetBranchAddress("jtSubJetPartonFlavor", &jtSubJetPartonFlavor);
			subjetTree->SetBranchAddress("jtSubJetHadronFlavor", &jtSubJetHadronFlavor);
			subjetTree->SetBranchAddress("jtSubJetHadronDR",&jtSubJetHadronDR);
			subjetTree->SetBranchAddress("jtSubJetPartonDR",&jtSubJetPartonDR);
			subjetTree->SetBranchAddress("jtSubJetPartonPt",&jtSubJetPartonPt);
			subjetTree->SetBranchAddress("jtSubJetPartonPdg",&jtSubJetPartonPdg);
			subjetTree->SetBranchAddress("jtSubJetHadronPt",&jtSubJetHadronPt);
			subjetTree->SetBranchAddress("jtSubJetHadronEta",&jtSubJetHadronEta);
			subjetTree->SetBranchAddress("jtSubJetHadronPhi",&jtSubJetHadronPhi);
			subjetTree->SetBranchAddress("jtSubJetHadronPdg",&jtSubJetHadronPdg);
		}
		subjetTree->SetBranchAddress("jtSubJetVtxType",&jtSubJetVtxType);
		
		for(int ievt=0; ievt<subjetTree->GetEntries(); ievt++){
						
			if(ievt%10000==0) cout << "at entry " << ievt << " of " << subjetTree->GetEntries() << endl;
			subjetTree->GetEntry(ievt);
			skimTree->GetEntry(ievt);
			evtTree->GetEntry(ievt);
			
			if(ispp && (!pvFilter || !beamScrapeFilter || !noiseFilter)) continue;
			if(!ispp && (!collEvtSel || !pvFilter || !hfCoincFilter || !noiseFilter)) continue;
			
			if(nref>MAXJETS) cout << "OH NO!!!" << endl;
			for(int ijet=0; ijet<nref; ijet++){
				
				if(jtpt[ijet]<30 || abs(jteta[ijet])>2.0) continue;
				if(jtSubJetPt->at(ijet).size()<2) continue;
				//if(jtSubJetPt->at(ijet).at(0) < 20 || jtSubJetPt->at(ijet).at(1) < 20) continue;
				
				t_jtpt = jtpt[ijet];
				t_jteta = jteta[ijet];
				t_jtphi = jtphi[ijet];
				t_refpt = refpt[ijet];
				t_discr_csvV2 = discr_csvV2[ijet];
				t_jtSubJetPt1 = jtSubJetPt->at(ijet).at(0);
				t_jtSubJetPt2 = jtSubJetPt->at(ijet).at(1);
				t_jtSubJetEta1 = jtSubJetEta->at(ijet).at(0);
				t_jtSubJetEta2 = jtSubJetEta->at(ijet).at(1);
				t_jtSubJetPhi1 = jtSubJetPhi->at(ijet).at(0);
				t_jtSubJetPhi2 = jtSubJetPhi->at(ijet).at(1);
				t_jtSubJetCSV1 = jtSubJetCSV->at(ijet).at(0);
				t_jtSubJetCSV2 = jtSubJetCSV->at(ijet).at(1);
				t_jtSubJetJP1 = jtSubJetJP->at(ijet).at(0);
				t_jtSubJetJP2 = jtSubJetJP->at(ijet).at(1);
				t_minCSV = t_jtSubJetCSV1 > t_jtSubJetCSV2 ? t_jtSubJetCSV2 : t_jtSubJetCSV1;
				t_jtSubJetdR = getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1));
				t_jtSubJetVtxType1 = jtSubJetVtxType->at(ijet).at(0);
				t_jtSubJetVtxType2 = jtSubJetVtxType->at(ijet).at(1);
				if(isMC){
					t_jtHadronFlavor = jtHadronFlavor[ijet];
					t_jtSubJetHadronFlavor1 = jtSubJetHadronFlavor->at(ijet).at(0);
					t_jtSubJetHadronFlavor2 = jtSubJetHadronFlavor->at(ijet).at(1);
				
                                        for(int isubjet=0; isubjet<2; isubjet++){
                                            hadronDR=999;
                                            for(int ihadron=0; ihadron<jtSubJetHadronDR->at(ijet).at(isubjet).size(); ihadron++){
                                                int hpdg = abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron));
                                                if(jtSubJetHadronFlavor->at(ijet).at(isubjet)==5){
                                                    if(hpdg>500 && hpdg<600 && jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron)<t_jtSubJetHadronDR1) t_jtSubJetHadronDR1 = jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron);
                                                }
                                                else if(jtSubJetHadronFlavor->at(ijet).at(isubjet)==4){
                                                    if(hpdg>400 && hpdg<500 && jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron)<t_jtSubJetHadronDR1) t_jtSubJetHadronDR1 = jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron);
                                                }
                                            }
                                            if(isubjet==0) t_jtSubJetHadronDR1 = hadronDR;
                                            else t_jtSubJetHadronDR2 = hadronDR;
                                        }
                                }
				gspTree->Fill();
				
				double zg = getzg(jtSubJetPt->at(ijet).at(0), jtSubJetPt->at(ijet).at(1));
				
				int ibin=0;
				while(t_jtpt>ptBins[ibin+1]) ibin++;
				zgSpec_ppIncl[ibin]->Fill(zg);
				if(t_minCSV > 0.9) zgSpec_ppGSP[ibin]->Fill(zg);
				
				if(isMC){
					
					int icentBin=0;
					if(!ispp){
						while(hiBin>centBins[icentBin+1]) icentBin++;
					}
					
					int flavCat = getFlavorCategory(t_jtSubJetHadronFlavor1, t_jtSubJetHadronFlavor2);
					int bin = flavorFrac[flavCat]->FindBin(t_minCSV);
					for(int ibin=1; ibin<=bin; ibin++){ 
						if(1) flavorFrac[flavCat]->AddBinContent(ibin);
						if(t_jtSubJetdR > 0.1) flavorFrac_withCuts[flavCat]->AddBinContent(ibin);
					}
					
					//build ROC curves vs CSV
					if(flavCat==0) allDiBjetsCut[0][icentBin]->Fill(t_minCSV);
					if(1) allJetsCut[0][icentBin]->Fill(t_minCSV);
					
					if(flavCat==0 && t_jtSubJetdR>0.1) allDiBjetsCut[1][icentBin]->Fill(t_minCSV);
					if(t_jtSubJetdR>0.1) allJetsCut[1][icentBin]->Fill(t_minCSV);
					
					if(flavCat==0 && t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allDiBjetsCut[2][icentBin]->Fill(t_minCSV);
					if(t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allJetsCut[2][icentBin]->Fill(t_minCSV);
					
					if(flavCat==0 && t_jtSubJetdR>0.1 && t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allDiBjetsCut[3][icentBin]->Fill(t_minCSV);
					if(t_jtSubJetdR>0.1 && t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allJetsCut[3][icentBin]->Fill(t_minCSV);
					
					if(flavCat==0) allDiBjets[icentBin]->Fill(t_minCSV);
					
					//build efficiency vs dR (csv>0.9)
					if(flavCat==0) allDiBjets_dR->Fill(t_jtSubJetdR);
					if(flavCat==0 && t_minCSV>0.9) allDBjetsCut_dR->Fill(t_jtSubJetdR);
					
					//build efficiency vs zg (csv>0.9)
					if(flavCat==0) allDiBjets_zg->Fill(zg);
					if(flavCat==0 && t_minCSV>0.9) allDBjetsCut_zg->Fill(zg);
					
					//correlation
					if(flavCat==0) allDiBjets_corr->Fill(zg, t_jtSubJetdR);
					if(flavCat==0 && t_minCSV>0.9) allDiBjetsCut_corr->Fill(zg, t_jtSubJetdR);
				}
			}
		}
		fin->Close();
	}
	double efficiency[4][nCent][40], purity[4][nCent][40];
	if(isMC){
		for(int icent=0; icent<nCent; icent++){
			double dib= allDiBjets[icent]->Integral();
			for(int hist=0; hist<4; hist++){
				for(int i=0; i<allJetsCut[hist][icent]->GetNbinsX(); i++){
					double allj=0, dibCut=0;
					for(int j=i; j<allJetsCut[hist][icent]->GetNbinsX(); j++){
						allj+= allJetsCut[hist][icent]->GetBinContent(j+1);
						dibCut+= allDiBjetsCut[hist][icent]->GetBinContent(j+1);
					}
					cout << "icent: "<< icent << " hist: "<< hist << " CSV>" << (double)i/40. << " allj: "<< allj << " dib: "<< dib << " dibCut: "<< dibCut << endl;
					efficiency[hist][icent][i] = dibCut/dib;
					purity[hist][icent][i] = dibCut/allj;
				}
			}
		}
	}
	
	cout << " purity, efficiency, CSV > X" << endl;
	for(int i=0; i<40; i++){
		cout << purity[0][0][i] << ", " << efficiency[0][0][i] << ", " << (double)i/40. << endl;
	}
	
	TGraph *purEffCurve_noCuts[nCent], *purEffCurve_drCut[nCent], *purEffCurve_vtxTypeCut[nCent], *purEffCurve_vtxType_drCut[nCent];
	for(int icent=0; icent<nCent; icent++){
		purEffCurve_noCuts[icent] = new TGraph(40, purity[0][icent], efficiency[0][icent]);
		purEffCurve_noCuts[icent]->SetNameTitle(Form("purEffCurve_noCuts_cent%d",icent),Form("purEffCurve_noCuts_cent%d",icent));
		purEffCurve_noCuts[icent]->SetLineColor(kGreen);
		purEffCurve_drCut[icent] = new TGraph(40, purity[1][icent], efficiency[1][icent]);
		purEffCurve_drCut[icent]->SetNameTitle(Form("purEffCurve_drCut_cent%d",icent),Form("purEffCurve_drCut_cent%d",icent));
		purEffCurve_drCut[icent]->SetLineColor(kGreen+2);
		purEffCurve_vtxTypeCut[icent] = new TGraph(40, purity[2][icent], efficiency[2][icent]);
		purEffCurve_vtxTypeCut[icent]->SetNameTitle(Form("purEffCurve_vtxTypeCut_cent%d",icent),Form("purEffCurve_vtxTypeCut_cent%d",icent));
		purEffCurve_vtxTypeCut[icent]->SetLineColor(kBlue-1);
		purEffCurve_vtxType_drCut[icent] = new TGraph(40, purity[3][icent], efficiency[3][icent]);
		purEffCurve_vtxType_drCut[icent]->SetNameTitle(Form("purEffCurve_vtxType_drCut_cent%d",icent),Form("purEffCurve_vtxType_drCut_cent%d",icent));
		purEffCurve_vtxType_drCut[icent]->SetLineColor(kBlue+2);
	}
	
	allDBjetsCut_dR->Divide(allDiBjets_dR);
	allDBjetsCut_zg->Divide(allDiBjets_zg);
	
	allDiBjetsCut_corr->Divide(allDiBjets_corr);
	
	fout->cd();
	gspTree->Write();
	for(int i=0; i<nPt; i++){
		zgSpec_ppIncl[i]->Write();
		zgSpec_ppGSP[i]->Write();
	}
	THStack *ths = new THStack("flavorFrac","flavorFrac");
	THStack *thsCuts = new THStack("flavorFrac_withCuts","flavorFrac_withCuts");
	formatHistos(flavorFrac);
	formatHistos(flavorFrac_withCuts);
	if(isMC){
		for(int i=5; i>=0; i--){
			ths->Add(flavorFrac[i]);
			thsCuts->Add(flavorFrac_withCuts[i]);
			flavorFrac[i]->Write();
			flavorFrac_withCuts[i]->Write();
		}
		ths->Write();
		thsCuts->Write();
		for(int icent=0; icent<nCent; icent++){
			purEffCurve_noCuts[icent]->Write();
			purEffCurve_drCut[icent]->Write();
			purEffCurve_vtxTypeCut[icent]->Write();
			purEffCurve_vtxType_drCut[icent]->Write();
		
			for(int i=0; i<4; i++){
				allDiBjetsCut[i][icent]->Write();
				allJetsCut[i][icent]->Write();
			}
			allDiBjets[icent]->Write();
		}
		allDBjetsCut_dR->Write();
		allDBjetsCut_zg->Write();
		allDiBjetsCut_corr->Write();
		allDiBjets_corr->Write();
	}
	fout->Close();
	
	
}
