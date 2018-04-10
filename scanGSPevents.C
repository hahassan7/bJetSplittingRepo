
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

//const int pthatBins[8] = {50,80,100,120,170,220,280,99999};
//const double xsecs[8] = {3.778E-03, 4.412E-04, 1.511E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 0};
//const int pthatEntries[7] = {176770, 195948, 193735, 191889, 179361, 194037, 265958}; //all pthats
//const int pthatEntries[7] = {740218, 693350, 610162, 1054000, 729448, 182013, 61566}; //big qcd-jet production (no full jets)
//const int pthatEntries[7] = {574183, 480754, 526452, 1067910, 820217, 204831, 69108}; //big qcd-jet production with full jets
//const int pthatEntries[7] = {271792, 558570, 500232, 985696, 795315, 198557, 67224}; //big qcd-jet production with subjet PF cand id's
//const int pthatEntries[7] = {171494, 492437, 560897, 1054710, 730146, 182688, 61740}; //big qcd-jet prod with 3.9% dropped tracks

//const int pthatEntries[7] = {439251, 366255, 411245, 651466, 475354, 265898, 87628}; //all pthat, pythia 6
//const int pthatEntries[8] = {0, 439251, 366255, 411245, 651466, 475354, 235375, 77559}; //pythia 6, 3 pthat220 jobs failed

//const int pthatEntries[7] = {122407, 157171, 211615, 410730, 366838, 374814, 145740}; //mu filtered samples
//const double mujetFracs[7] = {0.0247, 0.0375, 0.0480, 0.0530, 0.0691, 0.0691, 0.0691}; 
const double mujetFracs[7] = {0.133, 0.137, 0.136, 0.137, 0.139, 0.139, 0.139};

const double bjetFracs[7] = {0.058, 0.064, 0.075, 0.084, 0.084, 0.084, 0.084};
//const int bjetPthatEntries[7] = {470484, 360543, 116334, 478550, 452606, 112199, 36883}; //pythia 6 bjets without a pthat100 sample
const int bjetPthatEntries[7] = {209384, 144644, 161036, 264249, 193729, 195918, 68702}; //pythia 8 b-jets

const int pthatBins[2] = {50,99999}; //these three are for the extra bbbar prod
const double xsecs[2] = {1.697E-05, 0};
const int pthatEntries[1] = {81522};

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

void scanGSPevents(int nFiles=2000){
	
	const bool ispp = true;
	const bool isPowheg = false;
        const bool isMC = true;
        const bool doCaloJets = false;
        const bool isbJet = false;
        const bool doMuFilter = true;
        const bool doDMesonWeight = false;
        const bool useRefAxis = false;

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
	float discr_csvV2[MAXJETS], ndiscr_csvV2[MAXJETS], discr_prob[MAXJETS];
        float pthat;
	float vz;
        vector<vector<float> > *jtSubJetPt=0, *jtSubJetEta=0, *jtSubJetPhi=0;
        vector<vector<float> > *refSubJetPt=0, *refSubJetEta=0, *refSubJetPhi=0;
	vector<vector<float> > *jtSubJetPartonFlavor=0, *jtSubJetHadronFlavor=0;
	vector<vector<float> > *jtSubJetCSV=0, *jtSubJetJP=0, *jtSubJetNegCSV=0;
        vector<vector<float> > *jtSubJetRefPt=0, *jtSubJetRefEta=0, *jtSubJetRefPhi=0;
	int refparton_flavorProcess[MAXJETS];
	float jtHadronFlavor[MAXJETS];
	float jtm[MAXJETS];
	vector<vector<float> >*svtxm=0, *svtxpt=0, *svtxdl=0, *svtxdls=0;
	vector<vector<int> >*svtxntrk=0, *jtSubJetVtxType=0;

        vector<vector<int> > *jtSubJetConstitID=0;

	vector<vector<vector<float> > > *jtSubJetHadronDR=0, *jtSubJetPartonDR=0, *jtSubJetPartonPt=0, *jtSubJetPartonPdg=0;
	vector<vector<vector<float> > > *jtSubJetHadronPt=0, *jtSubJetHadronEta=0, *jtSubJetHadronPhi=0, *jtSubJetHadronPdg=0;
	
	vector<vector<vector<float> > > *jtSubJetSvtxm=0, *jtSubJetSvtxpt=0, *jtSubJetSvtxeta=0, *jtSubJetSvtxphi=0;
	vector<vector<vector<float> > > *jtSubJetSvtxNtrk=0, *jtSubJetSvtxdl=0, *jtSubJetSvtxdls=0;
	
        vector<float> *eta=0, *phi=0;
        vector<int> *pdg=0;
        vector<vector<int> > *motherIdx=0;
        int nparts;

	///*************************************************************//
	
	double t_jtpt, t_jteta, t_jtphi, t_matchedJtpt, t_discr_csvV2, t_ndiscr_csvV2, t_jtSubJetPt1, t_jtSubJetPt2, t_jtSubJetEta1, t_jtSubJetEta2, t_refpt;
	double t_jtSubJetPhi1, t_jtSubJetPhi2, t_jtSubJetCSV1, t_jtSubJetCSV2, t_jtSubJetdR, t_discr_prob;
	double t_jtSubJetVtxType1, t_jtSubJetVtxType2, t_minCSV, t_jtSubJetJP1, t_jtSubJetJP2, t_minJP;
        double t_jtSubJetNegCSV1, t_jtSubJetNegCSV2;
        double t_jtSubJetSvtxm1, t_jtSubJetSvtxm2, t_svtxm;
	double t_jtSubJetHadronDR1, t_jtSubJetHadronDR2, t_jtSubJetHadronPt1, t_jtSubJetHadronPt2;
        int t_jtHadronFlavor, t_jtSubJetHadronFlavor1, t_jtSubJetHadronFlavor2;
        double t_jtSubJetRefPt1, t_jtSubJetRefPt2, t_jtSubJetRefEta1, t_jtSubJetRefEta2, t_jtSubJetRefPhi1, t_jtSubJetRefPhi2;
	int t_jtSubJetPFType1, t_jtSubJetPFType2;
        int hiBin, t_dMesonChannel;
        double weight;
        float weightF;
        double t_muPt;
        double t_jtSubJetRefdR;

	TFile *fout = new TFile("gspTreeOut.root","recreate");
	TTree *gspTree = new TTree("gspTree","");
	gspTree->Branch("jtpt",&t_jtpt);
	gspTree->Branch("jteta",&t_jteta);
	gspTree->Branch("jtphi",&t_jtphi);
	gspTree->Branch("matchedJtpt",&t_matchedJtpt);
        gspTree->Branch("discr_csvV2",&t_discr_csvV2);
	gspTree->Branch("ndiscr_csvV2",&t_ndiscr_csvV2);
        gspTree->Branch("discr_prob", &t_discr_prob);
        gspTree->Branch("svtxm",&t_svtxm);
	gspTree->Branch("jtSubJetPt1",&t_jtSubJetPt1);
	gspTree->Branch("jtSubJetPt2",&t_jtSubJetPt2);
	gspTree->Branch("jtSubJetEta1", &t_jtSubJetEta1);
	gspTree->Branch("jtSubJetPhi1", &t_jtSubJetPhi1);
	gspTree->Branch("jtSubJetEta2", &t_jtSubJetEta2);
	gspTree->Branch("jtSubJetPhi2", &t_jtSubJetPhi2);
	gspTree->Branch("jtSubJetcsvV2_1",&t_jtSubJetCSV1);
	gspTree->Branch("jtSubJetcsvV2_2",&t_jtSubJetCSV2);
        gspTree->Branch("jtSubJetNegCsvV2_1",&t_jtSubJetNegCSV1);
        gspTree->Branch("jtSubJetNegCsvV2_2",&t_jtSubJetNegCSV2);
	gspTree->Branch("jtSubJetJP_1",&t_jtSubJetJP1);
	gspTree->Branch("jtSubJetJP_2",&t_jtSubJetJP2);
	gspTree->Branch("weight",&weight);
        gspTree->Branch("minCSV",&t_minCSV);
	gspTree->Branch("minJP",&t_minJP);
	gspTree->Branch("jtSubJetdR",&t_jtSubJetdR);
	gspTree->Branch("jtSubJetRefdR",&t_jtSubJetRefdR);
        gspTree->Branch("jtSubJetSvtxm1",&t_jtSubJetSvtxm1);
        gspTree->Branch("jtSubJetSvtxm2",&t_jtSubJetSvtxm2);
        gspTree->Branch("jtSubJetVtxType1",&t_jtSubJetVtxType1);
	gspTree->Branch("jtSubJetVtxType2",&t_jtSubJetVtxType2);
        gspTree->Branch("jtSubJetPFType1",&t_jtSubJetPFType1);
        gspTree->Branch("jtSubJetPFType2",&t_jtSubJetPFType2);
        gspTree->Branch("muPt",&t_muPt);
        gspTree->Branch("DMesonChannel",&t_dMesonChannel);
        if(isMC){
		gspTree->Branch("refpt",&t_refpt);
		if(!doCaloJets) gspTree->Branch("jtHadronFlavor",&t_jtHadronFlavor);
                else gspTree->Branch("refparton_flavorForB",&t_jtHadronFlavor);
		gspTree->Branch("jtSubJetHadronFlavor1",&t_jtSubJetHadronFlavor1);
		gspTree->Branch("jtSubJetHadronFlavor2",&t_jtSubJetHadronFlavor2);
                gspTree->Branch("jtSubJetHadronDR1",&t_jtSubJetHadronDR1);
                gspTree->Branch("jtSubJetHadronDR2",&t_jtSubJetHadronDR2);
                gspTree->Branch("jtSubJetHadronPt1",&t_jtSubJetHadronPt1);
                gspTree->Branch("jtSubJetHadronPt2",&t_jtSubJetHadronPt2);
                gspTree->Branch("jtSubJetRefPt1",&t_jtSubJetRefPt1);
                gspTree->Branch("jtSubJetRefPt2",&t_jtSubJetRefPt2);
                gspTree->Branch("jtSubJetRefEta1",&t_jtSubJetRefEta1);
                gspTree->Branch("jtSubJetRefEta2",&t_jtSubJetRefEta2);
                gspTree->Branch("jtSubJetRefPhi1",&t_jtSubJetRefPhi1);
                gspTree->Branch("jtSubJetRefPhi2",&t_jtSubJetRefPhi2);
        }
	if(!ispp){
		gspTree->Branch("hiBin",&hiBin);
	}
	
	
        std::string filelist;
        //if(ispp && isMC && !isbJet) filelist = "ppMC_pythia8_fullSet_muFiltered_withFullJets.txt";
        //if(isPowheg) filelist = "powheg_pp5TeV_List.txt";
        if(isPowheg) filelist = "powhegProdKurt.txt";
        //else if(ispp && isMC && !isbJet) filelist = "Pythia8_QCDjetMC_PFcandID.txt";
        //else if(ispp && isMC && !isbJet) filelist = "ppMC_pythia8_fullSet_muFiltered_withFullJets.txt";
        //if(ispp && isMC && !isbJet) filelist = "ppMC_pythia6_officialQCDSample.txt";
        else if(ispp && isMC && !isbJet) filelist = "pp_pythia8_BHadronPt50List.txt";
        else if(ispp && isMC) filelist = "Pythia8_bFilteredList.txt";
        else if(isMC) filelist = "PbPbMC_pthat120List_caloJet.txt";
        else if(ispp) filelist = "ppData_Jet80Trg_withMuons.txt";
        //else if(isMC) filelist = "ppMC_bjet_allPthat_withJPtags.txt";
        else filelist = "PbPbData_jet80TrgList_part2.txt";
        //else{ cout << "Not doing PbPb data yet! " << endl; exit(0); }
        //std::string filename = "HiForestAOD.root";
        std::string filename;
	ifstream instr(filelist.c_str());
	TFile *fin;
	int ifile=0;
	//for(int ifile=0; ifile<1; ifile++){
	cout << "warning! WEIGHT is turned off!" << endl;
        while(instr>>filename && ifile <nFiles){
		
		ifile++;
		fin = TFile::Open(filename.c_str());
		cout << "starting file " << filename << "... " << endl;
		
		if(!fin){ cout << "WARNING! Can't find file" << endl; continue; }
		
		int pvFilter, beamScrapeFilter, noiseFilter, trigger;
		int collEvtSel, hfCoincFilter;
		float pthat;
		TTree *hltTree = (TTree*)fin->Get("hltanalysis/HltTree");
                if(isMC) hltTree->SetBranchAddress("HLT_AK4PFJet80_Eta5p1ForPPRef_v1",&trigger);
                else hltTree->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1",&trigger);
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
                evtTree->SetBranchAddress("vz",&vz);

                TTree *genTree;
                if(isMC){
                    genTree = (TTree*)fin->Get("HiGenParticleAna/hi");
                    genTree->SetBranchAddress("mult",&nparts);
                    genTree->SetBranchAddress("eta",&eta);
                    genTree->SetBranchAddress("phi",&phi);
                    genTree->SetBranchAddress("pdg",&pdg);
                    genTree->SetBranchAddress("motherIdx",&motherIdx);
                }

		TTree *subjetTree;
		if(!ispp && doCaloJets) subjetTree = (TTree*)fin->Get("akPu4CaloJetAnalyzer/t");
                else if(!ispp) subjetTree = (TTree*)fin->Get("akCsSoftDrop4PFJetAnalyzer/t");
		else if(isMC) subjetTree = (TTree*)fin->Get("akSoftDrop4PFJetAnalyzer/t");
		else subjetTree = (TTree*)fin->Get("akSoftDrop4PFJetAnalyzer/t");
                TTree *fullJetTree;
                if(!ispp && doCaloJets) fullJetTree = (TTree*)fin->Get("akPu4CaloJetAnalyzer/t");
                else if(!ispp) fullJetTree = (TTree*)fin->Get("akCs4PFJetAnalyzer/t");
                else fullJetTree = (TTree*)fin->Get("ak4PFJetAnalyzer/t");

                int nMu;
                float muPt[100], muEta[100], muPhi[100];
                int trkLayerWMeas[100], pixLayerWMeas[100], isArbitrated[100];
                float trkDxy[100], trkDz[100];
                TTree *muTree;
                if(doMuFilter){
                    muTree = (TTree*)fin->Get("hltMuTree/HLTMuTree");
                    muTree->SetBranchAddress("Glb_nptl",&nMu);
                    muTree->SetBranchAddress("Glb_pt",muPt);
                    muTree->SetBranchAddress("Glb_eta",muEta);
                    muTree->SetBranchAddress("Glb_phi",muPhi);
                    muTree->SetBranchAddress("Glb_isArbitrated",isArbitrated);
                    muTree->SetBranchAddress("Glb_trkLayerWMeas",&trkLayerWMeas);
                    muTree->SetBranchAddress("Glb_pixLayerWMeas",&pixLayerWMeas);
                    muTree->SetBranchAddress("Glb_trkDxy",trkDxy);
                    muTree->SetBranchAddress("Glb_trkDz",trkDz);
                }
    
                float fulljtpt[MAXJETS], fulljteta[MAXJETS], fulljtphi[MAXJETS];
                int fullnref;
                fullJetTree->SetBranchAddress("jtpt",fulljtpt);
                fullJetTree->SetBranchAddress("jteta",fulljteta);
                fullJetTree->SetBranchAddress("jtphi",fulljtphi);
                fullJetTree->SetBranchAddress("nref",&fullnref);

		subjetTree->SetBranchAddress("nref", &nref);
		subjetTree->SetBranchAddress("nsvtx", nsvtx);
		subjetTree->SetBranchAddress("jtpt", jtpt);
		subjetTree->SetBranchAddress("jteta", jteta);
		subjetTree->SetBranchAddress("jtphi", jtphi);
                subjetTree->SetBranchAddress("svtxm", &svtxm);
		if(isMC){
			subjetTree->SetBranchAddress("refpt",refpt);
			subjetTree->SetBranchAddress("pthat",&pthat);
			if(!doCaloJets) subjetTree->SetBranchAddress("jtHadronFlavor",jtHadronFlavor);
                        else subjetTree->SetBranchAddress("refparton_flavorForB",jtHadronFlavor);
                        if(!doCaloJets) subjetTree->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
			//subjetTree->SetBranchAddress("refparton_flavorProcess",refparton_flavorProcess);
		}
		subjetTree->SetBranchAddress("discr_csvV2",discr_csvV2);
		subjetTree->SetBranchAddress("ndiscr_csvV2",ndiscr_csvV2);
                subjetTree->SetBranchAddress("discr_prob",discr_prob);
	
                if(isMC) evtTree->SetBranchAddress("pthat",&pthat);
                if(isPowheg) subjetTree->SetBranchAddress("weight",&weightF);

		subjetTree->SetBranchAddress("jtSubJetPt", &jtSubJetPt);
                subjetTree->SetBranchAddress("jtSubJetEta", &jtSubJetEta);
                subjetTree->SetBranchAddress("jtSubJetPhi", &jtSubJetPhi);
                if(isMC && useRefAxis){    
                    subjetTree->SetBranchAddress("refSubJetPt", &refSubJetPt);
                    subjetTree->SetBranchAddress("refSubJetEta", &refSubJetEta);
                    subjetTree->SetBranchAddress("refSubJetPhi", &refSubJetPhi);
                }
		subjetTree->SetBranchAddress("jtSubJetcsvV2",&jtSubJetCSV);
                subjetTree->SetBranchAddress("jtSubJetNegCsvV2",&jtSubJetNegCSV);
		subjetTree->SetBranchAddress("jtSubJetJP",&jtSubJetJP);
                subjetTree->SetBranchAddress("jtSubJetSvtxm",&jtSubJetSvtxm);
		//subjetTree->SetBranchAddress("jtSubJetConstituentsId",&jtSubJetConstitID);
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
			hltTree->GetEntry(ievt);
                        subjetTree->GetEntry(ievt);
			skimTree->GetEntry(ievt);
			evtTree->GetEntry(ievt);
                        if(doMuFilter) muTree->GetEntry(ievt);
	                if(doDMesonWeight && isMC) genTree->GetEntry(ievt);
                        fullJetTree->GetEntry(ievt);

			if(ispp && (!pvFilter || !beamScrapeFilter || !noiseFilter)) continue;
			if(!ispp && (!collEvtSel || !pvFilter || !hfCoincFilter || !noiseFilter)) continue;
		
                        if(abs(vz)>15) continue;
                        if(!trigger) continue;

                        if(isPowheg){
                            if(weightF>1e6) continue;
                        }
			if(nref>MAXJETS) cout << "OH NO!!!" << endl;
	
                        //implement bijective matching to coincide with 16-006
                        int matchedGroomToFull[MAXJETS];
                        int matchedFullToGroom[MAXJETS];
                        for(int ijet=0; ijet<nref; ijet++){
                            if(jtpt[ijet]<80 || abs(jteta[ijet])>2.0) continue;
                            double maxDR = 0.4;
                            for(int matchedJet=0; matchedJet<fullnref; matchedJet++){
                                if(fulljtpt[matchedJet]<50 || abs(fulljteta[matchedJet])>2.2) continue; 
                                double dr = getDR(jteta[ijet],jtphi[ijet],fulljteta[matchedJet],fulljtphi[matchedJet]);
                                if(dr < maxDR){
                                    matchedGroomToFull[ijet] = matchedJet;
                                    maxDR = dr;
                                }
                            }
                        }
                        for(int matchedJet=0; matchedJet<fullnref; matchedJet++){
                            if(fulljtpt[matchedJet]<50 || abs(fulljteta[matchedJet])>2.2) continue;
                            double maxDR = 0.4;
                            for(int ijet=0; ijet<nref; ijet++){
                                if(jtpt[ijet]<80 || abs(jteta[ijet])>2.0) continue;
                                double dr = getDR(fulljteta[matchedJet],fulljtphi[matchedJet],jteta[ijet],jtphi[ijet]);
                                if(dr < maxDR){
                                    matchedFullToGroom[matchedJet] = ijet;
                                    maxDR = dr;
                                }
                            }
                        }
                        
                        for(int ijet=0; ijet<nref; ijet++){
				
				//if(jtpt[ijet]<80 || abs(jteta[ijet])>2.0) continue;
                                if(jtpt[ijet]<80 || abs(jteta[ijet])>2.0) continue;

                                int ibin=0;
                                if(isMC && 0 && !isPowheg){
                                    if(pthat<=80.1) continue;
                                    if(isbJet){
                                        while(pthat>pthatBins[ibin+1]) ibin++;
                                        weight = bjetFracs[ibin]*(xsecs[ibin]-xsecs[ibin+1])/bjetPthatEntries[ibin];
                                    }
                                    else{
                                        while(pthat>pthatBins[ibin+1]) ibin++;
                                        if(doMuFilter) weight = mujetFracs[ibin]*(xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
                                        else weight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
                                    }
                                }
                                else if(!isPowheg) weight=2.5*7.581e-06/445955.;
                                else weight = weightF;

                                int muSel=0;
                                if(doMuFilter){
                                    double dr=999;
                                    if(nMu<1) continue;
                                    int nPassMu=0;
                                    for(int imu=0; imu<nMu; imu++){
                                        //use "soft" muon definition
                                        if(muPt[imu]<4 || muPt[imu]>400) continue;
                                        if(!isArbitrated[imu] || trkLayerWMeas[imu] <=5 || pixLayerWMeas[imu] <= 0 || trkDxy[imu] >= 0.3 || trkDz[imu] >= 20) continue;
                                        nPassMu++;
                                        double dphi = TMath::ACos(TMath::Cos(jtphi[ijet]-muPhi[imu]));
                                        double deta = jteta[ijet] - muEta[imu];
                                        double drTemp = TMath::Sqrt(pow(dphi,2)+pow(deta,2));
                                        if(drTemp < dr){ dr = drTemp; muSel = imu; }
                                    }
                                    if(!nPassMu) continue;
                                    if(dr>0.4) continue;
                                }

                                t_dMesonChannel=-1;
                                double partonDR=999;
                                int partonIdx=-1;
                                if(doDMesonWeight && isMC){
                                    for(int ipart=0; ipart<nparts; ipart++){
                                        double dphi = TMath::ACos(TMath::Cos(jtphi[ijet]-phi->at(ipart)));
                                        double deta = jteta[ijet] - eta->at(ipart);
                                        double drTemp = TMath::Sqrt(pow(dphi,2)+pow(deta,2));
                                        if(pdg->at(ipart)==13 && drTemp<partonDR){
                                            int nextMother = ipart;
                                            unsigned int loop=0;
                                            while(motherIdx->at(nextMother).size()>0 && loop < 20){
                                                if(motherIdx->at(nextMother).at(0)<2){ break; }
                                                if(abs(pdg->at(motherIdx->at(nextMother).at(0)))>400 && abs(pdg->at(motherIdx->at(nextMother).at(0)))<500){
                                                    //cout << "mother Index: "<< motherIdx->at(nextMother).at(0) << " pdg: "<< pdg->at(motherIdx->at(nextMother).at(0)) << endl;
                                                    partonIdx=ipart;
                                                    partonDR=drTemp;
                                                    break;
                                                }
                                                else{
                                                    nextMother = motherIdx->at(nextMother).at(0);
                                                }
                                                loop++;
                                            }
                                        }
                                    }
                                    if(partonIdx>-1){
                                        if(abs(pdg->at(motherIdx->at(partonIdx).at(0)))==421) t_dMesonChannel=0;
                                        else if(abs(pdg->at(motherIdx->at(partonIdx).at(0)))==411) t_dMesonChannel=1;
                                        else if(abs(pdg->at(motherIdx->at(partonIdx).at(0)))==431) t_dMesonChannel=2;
                                        if(t_dMesonChannel==0) weight*=1.37*1.13;
                                        if(t_dMesonChannel==1) weight*=0.91*0.98;
                                        if(t_dMesonChannel==2) weight*=0.67*1.16;
                                    }
                                }

                                if(ijet == matchedFullToGroom[matchedGroomToFull[ijet]]) t_matchedJtpt = fulljtpt[matchedGroomToFull[ijet]];
                                else t_matchedJtpt = -1;
				t_jtpt = jtpt[ijet];
				t_jteta = jteta[ijet];
				t_jtphi = jtphi[ijet];
				t_refpt = refpt[ijet];
				t_discr_csvV2 = discr_csvV2[ijet];
                                t_ndiscr_csvV2 = ndiscr_csvV2[ijet];
 				t_discr_prob = discr_prob[ijet];
                                t_jtHadronFlavor = jtHadronFlavor[ijet];
                                if(svtxm->at(ijet).size()>0) t_svtxm = svtxm->at(ijet).at(0);
                                else t_svtxm = -999;
                                if(!doCaloJets){
                                    if(jtSubJetPt->at(ijet).size()>1 && (!isMC || !useRefAxis || refSubJetPt->at(ijet).size()>1)){
                                        t_jtSubJetPt1 = jtSubJetPt->at(ijet).at(0);
                                        t_jtSubJetPt2 = jtSubJetPt->at(ijet).at(1);
                                        t_jtSubJetEta1 = jtSubJetEta->at(ijet).at(0);
                                        t_jtSubJetEta2 = jtSubJetEta->at(ijet).at(1);
                                        t_jtSubJetPhi1 = jtSubJetPhi->at(ijet).at(0);
                                        t_jtSubJetPhi2 = jtSubJetPhi->at(ijet).at(1);
                                        if(isMC && useRefAxis){
                                            t_jtSubJetRefPt1 = refSubJetPt->at(ijet).at(0);
                                            t_jtSubJetRefPt2 = refSubJetPt->at(ijet).at(1);
                                            t_jtSubJetRefEta1 = refSubJetEta->at(ijet).at(0);
                                            t_jtSubJetRefEta2 = refSubJetEta->at(ijet).at(1);
                                            t_jtSubJetRefPhi1 = refSubJetPhi->at(ijet).at(0);
                                            t_jtSubJetRefPhi2 = refSubJetPhi->at(ijet).at(1);
                                        }
                                        t_jtSubJetCSV1 = jtSubJetCSV->at(ijet).at(0);
                                        t_jtSubJetCSV2 = jtSubJetCSV->at(ijet).at(1);
                                        t_jtSubJetNegCSV1 = jtSubJetNegCSV->at(ijet).at(0);
                                        t_jtSubJetNegCSV2 = jtSubJetNegCSV->at(ijet).at(1);
                                        t_jtSubJetJP1 = jtSubJetJP->at(ijet).at(0);
                                        t_jtSubJetJP2 = jtSubJetJP->at(ijet).at(1);
                                        t_minCSV = t_jtSubJetCSV1 > t_jtSubJetCSV2 ? t_jtSubJetCSV2 : t_jtSubJetCSV1;
                                        t_minJP = t_jtSubJetJP1 > t_jtSubJetJP2 ? t_jtSubJetJP2 : t_jtSubJetJP1;
                                        t_jtSubJetdR = getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1));
                                        if(isMC && useRefAxis) t_jtSubJetRefdR = getDR(refSubJetEta->at(ijet).at(0), refSubJetPhi->at(ijet).at(0), refSubJetEta->at(ijet).at(1), refSubJetPhi->at(ijet).at(1));
                                        else t_jtSubJetRefdR = -1;
                                        t_jtSubJetVtxType1 = jtSubJetVtxType->at(ijet).at(0);
                                        t_jtSubJetVtxType2 = jtSubJetVtxType->at(ijet).at(1);
                                        //t_jtSubJetPFType1 = jtSubJetConstitID->at(ijet).at(0);
                                        //t_jtSubJetPFType2 = jtSubJetConstitID->at(ijet).at(1);
                                        if(doMuFilter) t_muPt = muPt[muSel];
                                        if(jtSubJetSvtxm->at(ijet).at(0).size()>0) t_jtSubJetSvtxm1 = jtSubJetSvtxm->at(ijet).at(0).at(0); //get first svtx
                                        else t_jtSubJetSvtxm1 = -999;
                                        if(jtSubJetSvtxm->at(ijet).at(1).size()>0) t_jtSubJetSvtxm2 = jtSubJetSvtxm->at(ijet).at(1).at(0);
                                        else t_jtSubJetSvtxm2 = -999;
                                        if(isMC){
                                            /*if((int)jtSubJetRefPt->size()>ijet){
                                                if(jtSubJetRefPt->at(ijet).size()>0){
                                                    t_jtSubJetRefPt1 = jtSubJetRefPt->at(ijet).at(0);
                                                    t_jtSubJetRefEta1 = jtSubJetRefEta->at(ijet).at(0);
                                                    t_jtSubJetRefPhi1 = jtSubJetRefPhi->at(ijet).at(0);
                                                }
                                                else{
                                                    t_jtSubJetRefPt1 = -999;
                                                    t_jtSubJetRefEta1 = -999;
                                                    t_jtSubJetRefPhi1 = -999;
                                                    t_jtSubJetRefPt2 = -999;
                                                    t_jtSubJetRefEta2 = -999;
                                                    t_jtSubJetRefPhi2 = -999;
                                                }
                                                if(jtSubJetRefPt->at(ijet).size()>1){
                                                    t_jtSubJetRefPt2 = jtSubJetRefPt->at(ijet).at(1);
                                                    t_jtSubJetRefEta2 = jtSubJetRefEta->at(ijet).at(1);
                                                    t_jtSubJetRefPhi2 = jtSubJetRefPhi->at(ijet).at(1);
                                                }
                                                else{
                                                    t_jtSubJetRefPt2 = -999;
                                                    t_jtSubJetRefEta2 = -999;
                                                    t_jtSubJetRefPhi2 = -999;
                                                }
                                            }
                                            else{ 
                                                cout << "WARNING! Ref size < reco size!" << endl; 
                                                t_jtSubJetRefPt1 = -999;
                                                t_jtSubJetRefEta1 = -999;
                                                t_jtSubJetRefPhi1 = -999;
                                                t_jtSubJetRefPt2 = -999;
                                                t_jtSubJetRefEta2 = -999;
                                                t_jtSubJetRefPhi2 = -999;
                                            }*/
                                            t_jtSubJetHadronFlavor1 = jtSubJetHadronFlavor->at(ijet).at(0);
                                            t_jtSubJetHadronFlavor2 = jtSubJetHadronFlavor->at(ijet).at(1);

                                            //get subjet-to-hadron dR
                                            for(int isubjet=0; isubjet<2; isubjet++){
                                                double hadronDR=999, hadronPt=-1;
                                                for(unsigned int ihadron=0; ihadron<jtSubJetHadronDR->at(ijet).at(isubjet).size(); ihadron++){
                                                    int hpdg = abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron));
                                                    if(jtSubJetHadronFlavor->at(ijet).at(isubjet)==5){
                                                        if(hpdg>500 && hpdg<600 && jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron)<hadronDR){
                                                            hadronDR = jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron);
                                                            hadronPt = jtSubJetHadronPt->at(ijet).at(isubjet).at(ihadron);
                                                        }
                                                    }
                                                    else if(jtSubJetHadronFlavor->at(ijet).at(isubjet)==4){
                                                        if(hpdg>400 && hpdg<500 && jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron)<hadronDR){
                                                            hadronDR = jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron);
                                                            hadronPt = jtSubJetHadronPt->at(ijet).at(isubjet).at(ihadron);
                                                        }
                                                    }
                                                }
                                                if(isubjet==0){
                                                    t_jtSubJetHadronDR1 = hadronDR;
                                                    t_jtSubJetHadronPt1 = hadronPt;
                                                }
                                                else{
                                                    t_jtSubJetHadronDR2 = hadronDR;
                                                    t_jtSubJetHadronPt2 = hadronPt;
                                                }
                                            }

                                            if(t_jtSubJetdR>0.3 && t_minJP<0.01 && t_jtSubJetHadronFlavor1==5 && t_jtSubJetHadronFlavor2==5 && t_jtSubJetPt2/(t_jtSubJetPt1+t_jtSubJetPt2)<0.2){
                                                //test for pathological b-ass'n cases
                                                //std::cout << "found pathological case - investigating..." << endl;
                                                //cout << "subjet 1 pt: "<< t_jtSubJetPt1 << " " << t_jtSubJetEta1 << " " << t_jtSubJetPhi1 << " " << t_jtSubJetJP1 << " " << t_jtSubJetCSV1 << endl;
                                                // cout << "hadron pdg, pt, eta, phi" << endl;
                                                for(unsigned int iparton=0; iparton<jtSubJetHadronPt->at(ijet).at(0).size(); iparton++){
                                                    //cout << jtSubJetHadronPdg->at(ijet).at(0).at(iparton) << " " << jtSubJetHadronPt->at(ijet).at(0).at(iparton) << " " << jtSubJetHadronEta->at(ijet).at(0).at(iparton) << " " << jtSubJetHadronPhi->at(ijet).at(0).at(iparton) << endl;
                                                }
                                                //cout << "subjet 2 pt, eta, phi, jp, csv: "<< t_jtSubJetPt2 << " " << t_jtSubJetEta2 << " " << t_jtSubJetPhi2 << " " << t_jtSubJetJP2 << " " << t_jtSubJetCSV2 << endl;
                                                //cout << "hadron pdg, pt, eta, phi" << endl;
                                                for(unsigned int iparton=0; iparton<jtSubJetHadronPt->at(ijet).at(1).size(); iparton++){
                                                    //cout << jtSubJetHadronPdg->at(ijet).at(1).at(iparton) << " " << jtSubJetHadronPt->at(ijet).at(1).at(iparton) << " " << jtSubJetHadronEta->at(ijet).at(1).at(iparton) << " " << jtSubJetHadronPhi->at(ijet).at(1).at(iparton) << endl;
                                                }
                                            }
                                        }
                                    }
                                    else{
                                        t_jtSubJetPt1 = -999;
                                        t_jtSubJetPt2 = -999; 
                                        t_jtSubJetEta1 = -999;
                                        t_jtSubJetEta2 = -999;
                                        t_jtSubJetPhi1 = -999;
                                        t_jtSubJetPhi2 = -999;
                                        t_jtSubJetCSV1 = -999;
                                        t_jtSubJetCSV2 = -999;
                                        t_jtSubJetNegCSV1 = -999;
                                        t_jtSubJetNegCSV2 = -999;
                                        t_jtSubJetJP1 = -999;
                                        t_jtSubJetJP2 = -999;
                                        t_minCSV = -999;
                                        t_minJP = -999;
                                        t_jtSubJetdR = -999;
                                        t_jtSubJetVtxType1 = -999; 
                                        t_jtSubJetVtxType2 = -999;
                                        t_jtSubJetSvtxm1 = -999;
                                        t_jtSubJetSvtxm2 = -999;
                                        if(isMC){
                                            t_jtHadronFlavor = -999;
                                            t_jtSubJetHadronFlavor1 = -999; 
                                            t_jtSubJetHadronFlavor2 = -999;
                                            t_jtSubJetHadronDR1 = 999;
                                            t_jtSubJetHadronDR2 = 999;
                                            t_jtSubJetRefPt1 = -999;
                                            t_jtSubJetRefPt2 = -999;
                                            t_jtSubJetRefEta1 = -999;
                                            t_jtSubJetRefEta2 = -999;
                                            t_jtSubJetRefPhi1 = -999;
                                            t_jtSubJetRefPhi2 = -999;
                                        }
                                    }
                                }

				gspTree->Fill();
			
                                if(doCaloJets) continue;
                                if(jtSubJetPt->at(ijet).size()<2) continue;

				double zg = getzg(jtSubJetPt->at(ijet).at(0), jtSubJetPt->at(ijet).at(1));
				
				ibin=0;
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
					if(flavCat==0) allDiBjetsCut[0][icentBin]->Fill(t_discr_csvV2);
					if(1) allJetsCut[0][icentBin]->Fill(t_discr_csvV2);
					
					if(flavCat==0 && t_jtSubJetdR>0.1) allDiBjetsCut[1][icentBin]->Fill(t_discr_csvV2);
					if(t_jtSubJetdR>0.1) allJetsCut[1][icentBin]->Fill(t_discr_csvV2);
					
					if(flavCat==0 && t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allDiBjetsCut[2][icentBin]->Fill(t_discr_csvV2);
					if(t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allJetsCut[2][icentBin]->Fill(t_discr_csvV2);
					
					if(flavCat==0 && t_jtSubJetdR>0.1 && t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allDiBjetsCut[3][icentBin]->Fill(t_discr_csvV2);
					if(t_jtSubJetdR>0.1 && t_jtSubJetVtxType1==0 && t_jtSubJetVtxType2==0) allJetsCut[3][icentBin]->Fill(t_discr_csvV2);
					
					if(flavCat==0) allDiBjets[icentBin]->Fill(t_discr_csvV2);
					
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
                cout << " purity, efficiency, CSV > X" << endl;
                for(int i=0; i<40; i++){
                    cout << purity[0][0][i] << ", " << efficiency[0][0][i] << ", " << (double)i/40. << endl;
                }
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
