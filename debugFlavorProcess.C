#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "THStack.h"

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
	return(sqrt(pow(eta1-eta2,2)+pow(phi1-phi2,2)));
}

double getdPhi(double phi1, double phi2){
	double dphi = abs(phi2-phi1);
	while(dphi>2.*PI) dphi-=(2.*PI);
	return dphi;
}

double getzg(float j1, float j2){
	return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}

double getBinomialError(int bin, TH1D* first, TH1D *second){

	TH1D *htemp = (TH1D*)first->Clone("htemp");
	htemp->Divide(first,second,1,1,"B");
	return htemp->GetBinError(bin-1);

}

void debugFlavorProcess(){

	bool ispp = true;
	double subjetDRCut = 0.2;
	double subjetVtxDRCut = 999;
	bool rejectPseudovertex = false;
      
	TH1D *statusType = new TH1D("statusType","",100,0,100);
	TH1D *statusJet = new TH1D("statusJet","",100,0,100);
	TH1D *doublebStatusCorr = new TH1D("doublebStatusCorr","",100,0,100);
	
	TH1D *allJets = new TH1D("allJets","",30,0,1);
	TH1D *lightFrac = new TH1D("lightFrac","",30,0,1);
	TH1D *charmFrac = new TH1D("charmFrac","",30,0,1);
	TH1D *singleBFrac = new TH1D("singleBFrac","",30,0,1);
	TH1D *bcFrac = new TH1D("bcFrac","",30,0,1);
	TH1D *doubleBFrac = new TH1D("doubleBFrac","",30,0,1);
	TH1D *unknownFrac = new TH1D("unknownFrac","",30,0,1);
	
	TH1D *allJets_dR = new TH1D("allJets_dR","",30,0,1);
	TH1D *lightFrac_dR = new TH1D("lightFrac_dR","",30,0,1);
	TH1D *charmFrac_dR = new TH1D("charmFrac_dR","",30,0,1);
	TH1D *singleBFrac_dR = new TH1D("singleBFrac_dR","",30,0,1);
	TH1D *bcFrac_dR = new TH1D("bcFrac_dR","",30,0,1);
	TH1D *doubleBFrac_dR = new TH1D("doubleBFrac_dR","",30,0,1);
	TH1D *unknownFrac_dR = new TH1D("unknownFrac_dR","",30,0,1);
	
	TH1D *allJets_dR_subVtxOne = new TH1D("allJets_dR_subVtxOne","",30,0,1);
	TH1D *lightFrac_dR_subVtxOne = new TH1D("lightFrac_dR_subVtxOne","",30,0,1);
	TH1D *charmFrac_dR_subVtxOne = new TH1D("charmFrac_dR_subVtxOne","",30,0,1);
	TH1D *singleBFrac_dR_subVtxOne = new TH1D("singleBFrac_dR_subVtxOne","",30,0,1);
	TH1D *bcFrac_dR_subVtxOne = new TH1D("bcFrac_dR_subVtxOne","",30,0,1);
	TH1D *doubleBFrac_dR_subVtxOne = new TH1D("doubleBFrac_dR_subVtxOne","",30,0,1);
	TH1D *unknownFrac_dR_subVtxOne = new TH1D("unknownFrac_dR_subVtxOne","",30,0,1);
	
	TH1D *allJets_dR_subVtxBoth = new TH1D("allJets_dR_subVtxBoth","",30,0,1);
	TH1D *lightFrac_dR_subVtxBoth = new TH1D("lightFrac_dR_subVtxBoth","",30,0,1);
	TH1D *charmFrac_dR_subVtxBoth = new TH1D("charmFrac_dR_subVtxBoth","",30,0,1);
	TH1D *singleBFrac_dR_subVtxBoth = new TH1D("singleBFrac_dR_subVtxBoth","",30,0,1);
	TH1D *bcFrac_dR_subVtxBoth = new TH1D("bcFrac_dR_subVtxBoth","",30,0,1);
	TH1D *doubleBFrac_dR_subVtxBoth = new TH1D("doubleBFrac_dR_subVtxBoth","",30,0,1);
	TH1D *unknownFrac_dR_subVtxBoth = new TH1D("unknownFrac_dR_subVtxBoth","",30,0,1);
	
	TH2D *dRCSVCorrelation_singleB = new TH2D("dRCSVCorrelation_singleB","",30,0,1,50,0,0.5);
	TH2D *dRCSVCorrelation_doubleB = new TH2D("dRCSVCorrelation_doubleB","",30,0,1,50,0,0.5);
	TH2D *dRCSVCorrelation_doubleB_doubleType0 = new TH2D("dRCSVCorrelation_doubleB_doubleType0","",30,0,1,50,0,0.5);
	TH2D *dRCSVCorrelation_doubleB_singleType0 = new TH2D("dRCSVCorrelation_doubleB_singleType0","",30,0,1,50,0,0.5);
	TH2D *dRCSVCorrelation_doubleB_noType0 = new TH2D("dRCSVCorrelation_doubleB_noType0","",30,0,1,50,0,0.5);
	
	TH2D *csv12Correlation = new TH2D("csv12Correlation","",50,0,1,50,0,1);
	TH2D *csv12Correlation_doubleB = new TH2D("csv12Correlation_doubleB","",50,0,1,50,0,1);
	TH2D *csv12Correlation_singleB = new TH2D("csv12Correlation_singleB","",50,0,1,50,0,1);
	
	TH1D *csvDistr = new TH1D("csvDistr","",50,0,1);
	TH1D *csvDistr_bsubjet = new TH1D("csvDistr_bsubjet","",50,0,1);
	
	TH1D *subjetVtxToDR = new TH1D("subjetVtxToDR","",50,0,0.5);
	TH1D *subjetVtxToDR_csv09 = new TH1D("subjetVtxToDR_csv09","",50,0,0.5);
	TH1D *subjetVtxToDR_bFlav = new TH1D("subjetVtxToDR_bFlav","",50,0,0.5);
	
	TH1D *dRtoMatchedSubjet = new TH1D("dRtoMatchedSubjet","",50,0,0.5);
	TH1D *dRtoUnmatchedSubjet = new TH1D("dRtoUnmatchedSubjet","",50,0,0.5);
	
	TH2D *singleDoubleNVtx = new TH2D("singleDoubleNVtx","",10,0,10,2,1,3);
	
	double nDiBSubjets[40];
	double totalSubjets[40];	
	
	double nDiBVtxDRSubjets[40];
	double totalVtxDRSubjets[40];
	for(int i=0; i<40; i++){ totalSubjets[i]=0; nDiBSubjets[i]=0; totalVtxDRSubjets[i]=0; nDiBVtxDRSubjets[i]=0; }
	double totalDiBSubjets=0;
	
	//really only one b-hadron (one b-parton), one b-hadron (multiple b-partons), multiple b-hadrons in single subjet,
	//one b-hadron + one c+hadron
	
	TH1D *failureModes = new TH1D("failureModes","",4,1,5);
	
	TH1D *jetSVMass_sB = new TH1D("jetSVMass_sB","",100,0,40); jetSVMass_sB->Sumw2();
	TH1D *jetSVMass_dB = new TH1D("jetSVMass_dB","",100,0,40); jetSVMass_dB->Sumw2();
	
	TH1D *subjetSVsBMass = new TH1D("subjetSVsBMass","",30,0,6); subjetSVsBMass->Sumw2();
	TH1D *subjetSVdBMass = new TH1D("subjetSVdBMass","",30,0,6); subjetSVdBMass->Sumw2();
	
	//****************************************************************//
	
	int nref;
	float jtpt[MAXJETS], jteta[MAXJETS], jtphi[MAXJETS];
	int refparton_flavorForB[MAXJETS], nsvtx[MAXJETS];
	float discr_csvV2[MAXJETS];
	vector<vector<float> > *jtSubJetPt=0, *jtSubJetEta=0, *jtSubJetPhi=0;
	vector<vector<float> > *jtSubJetPartonFlavor=0, *jtSubJetHadronFlavor=0;
	vector<vector<float> > *jtSubJetCSV=0;
	int refparton_flavorProcess[MAXJETS];
	float jtHadronFlavor[MAXJETS];
	float jtm[MAXJETS];
	vector<vector<float> >*svtxm=0, *svtxpt=0, *svtxdl=0, *svtxdls=0;
	vector<vector<int> >*svtxntrk=0, *jtSubJetVtxType=0;
	
	vector<vector<vector<float> > > *jtSubJetHadronDR=0, *jtSubJetPartonDR=0, *jtSubJetPartonPt=0, *jtSubJetPartonPdg=0;
	vector<vector<vector<float> > > *jtSubJetHadronPt=0, *jtSubJetHadronEta=0, *jtSubJetHadronPhi=0, *jtSubJetHadronPdg=0;
	
	vector<vector<vector<float> > > *jtSubJetSvtxm=0, *jtSubJetSvtxpt=0, *jtSubJetSvtxeta=0, *jtSubJetSvtxphi=0;
	vector<vector<vector<float> > > *jtSubJetSvtxNtrk=0, *jtSubJetSvtxdl=0, *jtSubJetSvtxdls=0;
	
	vector<vector<int> > *motherIdx=0, *daughterIdx=0;
	
	vector<float> *pt=0, *eta=0, *phi=0;
	vector<int> *pdg=0, *sta=0;
	
	int nSingleJets=0;
	int nSingleJetsPSVtx=0, nSingleJetsPSVtxSingleB=0;
	int nSingleJetsDoublePSVtx=0;
	int nSingleJetsNoVtx=0;
	int nSingleJetsDoubleNoVtx=0;
	int nSingleJetsPSNoVtx=0;
	
	///*************************************************************//
	
	
	std::string filelist;
	if(ispp) filelist = "ppMC_qcdJetList.txt";
      else filelist = "PbPbMC_QCDjetList.txt";
      //std::string filename = "HiForestAOD.root";
	std::string filename;
	ifstream instr(filelist.c_str());
	TFile *fin;
	int ifile=0;
	//for(int ifile=0; ifile<1; ifile++){
	while(instr>>filename && ifile <5000){
		
		ifile++;
		cout << "top of file loop" << endl;
		cout << "infile loads correctly" << endl;
		fin = TFile::Open(filename.c_str());
		cout << "starting file " << filename << "... " << endl;
		
		TTree *subjetTree;
                if(!ispp) subjetTree = (TTree*)fin->Get("akCsSoftDrop4PFJetAnalyzer/t");
                else subjetTree = (TTree*)fin->Get("akSoftDrop4PFJetAnalyzer/t");
                TTree *partonTree = (TTree*)fin->Get("HiGenParticleAna/hi");
		
		subjetTree->SetBranchAddress("nref", &nref);
		subjetTree->SetBranchAddress("nsvtx", nsvtx);
		subjetTree->SetBranchAddress("jtpt", jtpt);
		subjetTree->SetBranchAddress("jteta", jteta);
		subjetTree->SetBranchAddress("jtphi", jtphi);
		subjetTree->SetBranchAddress("svtxm", &svtxm);
		subjetTree->SetBranchAddress("svtxpt", &svtxpt);
		subjetTree->SetBranchAddress("svtxdl", &svtxdl);
		subjetTree->SetBranchAddress("svtxdls", &svtxdls);
		subjetTree->SetBranchAddress("svtxntrk", &svtxntrk);
		subjetTree->SetBranchAddress("jtHadronFlavor",jtHadronFlavor);
		subjetTree->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
		subjetTree->SetBranchAddress("refparton_flavorProcess",refparton_flavorProcess);
		subjetTree->SetBranchAddress("discr_csvV2",discr_csvV2);
		
		subjetTree->SetBranchAddress("jtSubJetPt", &jtSubJetPt);
		subjetTree->SetBranchAddress("jtSubJetEta", &jtSubJetEta);
		subjetTree->SetBranchAddress("jtSubJetPhi", &jtSubJetPhi);
		subjetTree->SetBranchAddress("jtSubJetcsvV1",&jtSubJetCSV);
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
		
		subjetTree->SetBranchAddress("jtSubJetVtxType",&jtSubJetVtxType);
		subjetTree->SetBranchAddress("jtSubJetSvtxm",&jtSubJetSvtxm);
		subjetTree->SetBranchAddress("jtSubJetSvtxpt",&jtSubJetSvtxpt);
		subjetTree->SetBranchAddress("jtSubJetSvtxeta",&jtSubJetSvtxeta);
		subjetTree->SetBranchAddress("jtSubJetSvtxphi",&jtSubJetSvtxphi);
		subjetTree->SetBranchAddress("jtSubJetSvtxNtrk",&jtSubJetSvtxNtrk);
		subjetTree->SetBranchAddress("jtSubJetSvtxdl",&jtSubJetSvtxdl);
		subjetTree->SetBranchAddress("jtSubJetSvtxdls",&jtSubJetSvtxdls);
		
		partonTree->SetBranchAddress("pt",&pt);
		partonTree->SetBranchAddress("eta",&eta);
		partonTree->SetBranchAddress("phi",&phi);
		partonTree->SetBranchAddress("pdg",&pdg);
		partonTree->SetBranchAddress("motherIdx",&motherIdx);
		partonTree->SetBranchAddress("daughterIdx",&daughterIdx);
		partonTree->SetBranchAddress("sta",&sta);
		
		for(int ievt=0; ievt<subjetTree->GetEntries(); ievt++){
						
			if(ievt%1000==0) cout << "at entry " << ievt << " of " << subjetTree->GetEntries() << endl;
			subjetTree->GetEntry(ievt);
			partonTree->GetEntry(ievt);
	
			if(nref>MAXJETS) cout << "OH NO!!!" << endl;
			//cout << "jet size: "<< nref << endl;
			//cout << "parton size: "<< pt->size() << endl;
			for(int ijet=0; ijet<nref; ijet++){
				
				if(jtpt[ijet]<80 || abs(jteta[ijet])>2.0) continue;
				if(jtSubJetPt->at(ijet).size()<2) continue;
				if(jtSubJetPt->at(ijet).at(0) < 20 || jtSubJetPt->at(ijet).at(1) < 20) continue;
				
				bool dijetFlag=false;
				for(int ijet2=0; ijet2<nref && ijet!=ijet2; ijet2++){
					if(jtpt[ijet2]<80 || abs(jteta[ijet2])>2.0) continue;
					double dijetdPhi = getdPhi(jtphi[ijet], jtphi[ijet2]);
					if(dijetdPhi > 2*PI/3. && dijetdPhi < 4*PI/3. && discr_csvV2[ijet2]<0.5){ dijetFlag = true; break; }
				}
				if(!dijetFlag) continue;
				
				bool foundPseudovtx=false;
				if(rejectPseudovertex){
					if(jtSubJetVtxType->at(ijet).at(0)==1 || jtSubJetVtxType->at(ijet).at(1)==1) foundPseudovtx=true;
				}
				//cout << "at jet : "<< ijet << " pt: "<< jtpt[ijet] << " eta: "<< jteta[ijet] << endl;
				
				double csvArr1 = (jtSubJetCSV->at(ijet).at(0));
				double csvArr2 = (jtSubJetCSV->at(ijet).at(1));
				double avgCSV = (csvArr2+csvArr1)/2.;
				if(avgCSV<0) avgCSV=0;
				double minCSV = csvArr1 < csvArr2 ? csvArr1 : csvArr2;
				
				double subjetVtxDR1 = 999;
				double subjetVtxDR2 = 999;
				if(jtSubJetSvtxm->at(ijet).at(0).size()>0) subjetVtxDR1 = getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetSvtxeta->at(ijet).at(0).at(0), jtSubJetSvtxphi->at(ijet).at(0).at(0));
				if(jtSubJetSvtxm->at(ijet).at(1).size()>0) subjetVtxDR2 = getDR(jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1), jtSubJetSvtxeta->at(ijet).at(1).at(0), jtSubJetSvtxphi->at(ijet).at(1).at(0));
				
				double subjetDR = getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1));
				double subjetZg = getzg(jtSubJetPt->at(ijet).at(0), jtSubJetPt->at(ijet).at(0));
				
				csv12Correlation->Fill(csvArr1, csvArr2);
				if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 || abs(jtSubJetHadronFlavor->at(ijet).at(0))==5) csv12Correlation_singleB->Fill(csvArr1, csvArr2);
				if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 && abs(jtSubJetHadronFlavor->at(ijet).at(0))==5) csv12Correlation_doubleB->Fill(csvArr1, csvArr2);
				
				for(unsigned int isubjet=0; isubjet<jtSubJetPt->at(ijet).size(); isubjet++){
					for(unsigned int ivtx=0; ivtx<jtSubJetSvtxm->at(ijet).at(isubjet).size(); ivtx++){
						double vtxDR = getDR(jtSubJetEta->at(ijet).at(isubjet), jtSubJetPhi->at(ijet).at(isubjet), jtSubJetSvtxeta->at(ijet).at(isubjet).at(ivtx), jtSubJetSvtxphi->at(ijet).at(isubjet).at(ivtx));
						subjetVtxToDR->Fill(vtxDR);
						if(jtSubJetCSV->at(ijet).at(isubjet)>0.9) subjetVtxToDR_csv09->Fill(vtxDR);
						if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))==5) subjetVtxToDR_bFlav->Fill(vtxDR);
					}
					
					csvDistr->Fill(jtSubJetCSV->at(ijet).at(isubjet));
					if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))==5) csvDistr_bsubjet->Fill(jtSubJetCSV->at(ijet).at(isubjet));
				}
				
				allJets->Fill(minCSV);
				if(subjetDR>subjetDRCut){
					allJets_dR->Fill(minCSV);
					if(subjetVtxDR1<subjetVtxDRCut || subjetVtxDR2<subjetVtxDRCut) allJets_dR_subVtxOne->Fill(minCSV);
					if(subjetVtxDR1<subjetVtxDRCut && subjetVtxDR2<subjetVtxDRCut) allJets_dR_subVtxBoth->Fill(minCSV);
				}
				
				if(jtHadronFlavor[ijet]==4){
					charmFrac->Fill(minCSV);
					if(subjetDR>subjetDRCut){
						charmFrac_dR->Fill(minCSV);
						if(subjetVtxDR1<subjetVtxDRCut || subjetVtxDR2<subjetVtxDRCut) charmFrac_dR_subVtxOne->Fill(minCSV);
						if(subjetVtxDR1<subjetVtxDRCut && subjetVtxDR2<subjetVtxDRCut) charmFrac_dR_subVtxBoth->Fill(minCSV);
					}
				}
				else if(jtHadronFlavor[ijet]!=5){
					lightFrac->Fill(minCSV);
					if(subjetDR>subjetDRCut){
						lightFrac_dR->Fill(minCSV);
						if(subjetVtxDR1<subjetVtxDRCut || subjetVtxDR2<subjetVtxDRCut) lightFrac_dR_subVtxOne->Fill(minCSV);
						if(subjetVtxDR1<subjetVtxDRCut && subjetVtxDR2<subjetVtxDRCut) lightFrac_dR_subVtxBoth->Fill(minCSV);
					}
				}
				
				int csvBin = 0;
				if(minCSV>0){
					csvBin = minCSV*40;
				}
				
				double minSubjetVtxDR = subjetVtxDR1 > subjetVtxDR2 ? subjetVtxDR2 : subjetVtxDR1;
				int vtxDRBin=39;
				if(minSubjetVtxDR>0){
					vtxDRBin = minSubjetVtxDR*80.;	
				}
				if(vtxDRBin>=40 || vtxDRBin<0) vtxDRBin=39;
				
				for(int ibin=0; ibin<csvBin; ibin++){ if(subjetDR>subjetDRCut && minSubjetVtxDR<subjetVtxDRCut && (!rejectPseudovertex || !foundPseudovtx)) totalSubjets[ibin]++; }
				for(int ibin=vtxDRBin; ibin>=0; ibin--){ if(subjetDR>subjetDRCut && minSubjetVtxDR<subjetVtxDRCut && (!rejectPseudovertex || !foundPseudovtx)) totalVtxDRSubjets[ibin]++; }
				
				//cout << "cp1" << endl;
								
				if(jtHadronFlavor[ijet]!=5) continue;
					
				//cout << "cp2" << endl;
				
				int origQuarkStatus=-1;
				bool breakFlag=false;
				for(unsigned int iparton=0; iparton<pt->size(); iparton++){
					if(abs(pdg->at(iparton))==5){
						if(pt->at(iparton)<0.01) continue;
						if(getDR(jteta[ijet], jtphi[ijet], eta->at(iparton), phi->at(iparton)) < 0.4){
							if(motherIdx->at(iparton).at(0)==-999){
								origQuarkStatus = sta->at(iparton);
								break;
							}
							while(pt->at(motherIdx->at(iparton).at(0)) > 0.01 && abs(pdg->at(motherIdx->at(iparton).at(0)))==5){
								iparton = motherIdx->at(iparton).at(0);
								if(iparton > pt->size()){ breakFlag = true; break; }
								if(iparton == (unsigned int)motherIdx->at(motherIdx->at(iparton).at(0)).at(0)) { breakFlag = true; break; }
							}
							if(!breakFlag){
								//cout << "parton id: "<< iparton << " pdg: "<< pdg->at(iparton) << " pt: "<< pt->at(iparton) << " motherIdx: " << motherIdx->at(iparton).at(0) <<  " motherPdg: "<< pdg->at(motherIdx->at(iparton).at(0)) << endl;
								origQuarkStatus = sta->at(iparton);
								break;
							}
							else break;
						}
					}
				}
				statusType->Fill(origQuarkStatus);
				
				//cout << "cp3" << endl;
				
				
				if(refparton_flavorProcess[ijet]==5 || refparton_flavorProcess[ijet]==6){
					statusJet->Fill(origQuarkStatus);
				}
				
			//check out the subjets...
				unsigned int flavorCount=0;
				unsigned int cCount=0;
				for(unsigned int isubjet=0; isubjet<jtSubJetHadronFlavor->at(ijet).size(); isubjet++){
					if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))==5) flavorCount++;
					else if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))==4) cCount++;
				}
				
				//cout << "cp4" << endl;
				//cout << "flavorCount: "<< flavorCount << endl;
								
				if(flavorCount>1){					
					doublebStatusCorr->Fill(origQuarkStatus);
					doubleBFrac->Fill(minCSV);
					if(subjetDR>subjetDRCut){
						doubleBFrac_dR->Fill(minCSV);
						if(subjetVtxDR1<subjetVtxDRCut || subjetVtxDR2<subjetVtxDRCut) doubleBFrac_dR_subVtxOne->Fill(minCSV);
						if(subjetVtxDR1<subjetVtxDRCut && subjetVtxDR2<subjetVtxDRCut) doubleBFrac_dR_subVtxBoth->Fill(minCSV);
					}
					dRCSVCorrelation_doubleB->Fill(minCSV, subjetDR);
					if(jtSubJetVtxType->at(ijet).at(0)==0 && jtSubJetVtxType->at(ijet).at(1)==0) dRCSVCorrelation_doubleB_doubleType0->Fill(minCSV, subjetDR);
					if(jtSubJetVtxType->at(ijet).at(0)==0 || jtSubJetVtxType->at(ijet).at(1)==0) dRCSVCorrelation_doubleB_singleType0->Fill(minCSV, subjetDR);
					if(jtSubJetVtxType->at(ijet).at(0)!=0 && jtSubJetVtxType->at(ijet).at(1)!=0) dRCSVCorrelation_doubleB_noType0->Fill(minCSV, subjetDR);
					jetSVMass_dB->Fill(svtxm->at(ijet).at(0));
					singleDoubleNVtx->Fill(nsvtx[ijet], 2);
					if(subjetDR>subjetDRCut && minSubjetVtxDR<subjetVtxDRCut && (!rejectPseudovertex || !foundPseudovtx)){
						for(int ibin=0; ibin<csvBin; ibin++){ 
							nDiBSubjets[ibin]++; 
						}
						for(int ibin=vtxDRBin; ibin>=0; ibin--){
							nDiBVtxDRSubjets[ibin]++;
						}
					}
					totalDiBSubjets++;					
						
					//cout << "nsubjets: "<< jtSubJetCSV->at(ijet).size() << endl;
					//cout << "pt subjets:" << jtSubJetPt->at(ijet).size() << endl;
					//cout << "hadron subjets:" << jtSubJetHadronFlavor->at(ijet).size() << endl;
					//cout << "vertex nsubjets: "<< jtSubJetSvtxm->at(ijet).size() << endl;
					for(unsigned int isubjet=0; isubjet<jtSubJetPt->at(ijet).size(); isubjet++){
						//cout << "hadron flavor subjet: "<< isubjet << " flavor: "<< jtSubJetHadronFlavor->at(ijet).at(isubjet) << endl;
						//cout << "nvertices: " << jtSubJetSvtxm->at(ijet).at(isubjet).size() << endl;
						for(unsigned int isubvtx=0; isubvtx<jtSubJetSvtxm->at(ijet).at(isubjet).size(); isubvtx++){
							if(isubjet==1) subjetSVdBMass->Fill(jtSubJetSvtxm->at(ijet).at(isubjet).at(isubvtx));
						}
					}
				}
				else if(flavorCount==1){
					if(cCount) bcFrac->Fill(minCSV);
					else singleBFrac->Fill(minCSV);
					if(subjetDR>subjetDRCut){
						if(cCount){
							bcFrac_dR->Fill(minCSV);
							if(subjetVtxDR1<subjetVtxDRCut || subjetVtxDR2<subjetVtxDRCut) bcFrac_dR_subVtxOne->Fill(minCSV);
							if(subjetVtxDR1<subjetVtxDRCut && subjetVtxDR2<subjetVtxDRCut) bcFrac_dR_subVtxBoth->Fill(minCSV);
						}
						else{
							singleBFrac_dR->Fill(minCSV);	
							if(subjetVtxDR1<subjetVtxDRCut || subjetVtxDR2<subjetVtxDRCut) singleBFrac_dR_subVtxOne->Fill(minCSV);
							if(subjetVtxDR1<subjetVtxDRCut && subjetVtxDR2<subjetVtxDRCut) singleBFrac_dR_subVtxBoth->Fill(minCSV);
						}						
					}
					dRCSVCorrelation_singleB->Fill(minCSV, subjetDR);
					jetSVMass_sB->Fill(svtxm->at(ijet).at(0));
					
					singleDoubleNVtx->Fill(nsvtx[ijet], 1);
					
					//cout << "inside flavor1" << endl;
					//cout << "sj dr: "<< subjetDR << " minCsv: "<< minCSV << endl;
					
					if(subjetDR>subjetDRCut && minCSV > 0.9 && (!rejectPseudovertex || !foundPseudovtx)) {
						if(jtSubJetVtxType->at(ijet).at(0)==1 || jtSubJetVtxType->at(ijet).at(1)==1){
							
							cout << "found single-B-tagged jet system with big DR and CSV " << endl;
							nSingleJets++;
							cout << "in event : "<< ievt << endl;
							cout << "full jet info: " << endl;
							cout << "jtpt: "<< jtpt[ijet] << " eta: "<< jteta[ijet] << " phi: "<< jtphi[ijet] << " flavor: "<< refparton_flavorForB[ijet] << " hadron flavor " << jtHadronFlavor[ijet] << endl;
							cout <<" full jet SV info: " << endl;
							for(unsigned int isv=0; isv<svtxm->at(ijet).size(); isv++){
								cout << "mass: "<< svtxm->at(ijet).at(isv) << " pt: "<< svtxpt->at(ijet).at(isv) << " dl "<< svtxdl->at(ijet).at(isv) << " dls: "<< svtxdls->at(ijet).at(isv) << " ntrk: " << svtxntrk->at(ijet).at(isv) << endl;
							}
							
							cout << "csv1: "<< csvArr1 << " csv2: "<< csvArr2 << endl;
							
							for(unsigned int iparton=0; iparton<50 && iparton<pt->size(); iparton++){
								if(abs(pdg->at(iparton))==5){
									if(pt->at(iparton)<0.01) continue;
									if(getDR(jteta[ijet], jtphi[ijet], eta->at(iparton), phi->at(iparton)) < 0.6){
										cout << "parton id: "<< iparton << " pdg: "<< pdg->at(iparton) << " pt: "<< pt->at(iparton) << " eta: "<< eta->at(iparton) << " phi: "<< phi->at(iparton) << endl;
									}
								}
							}
							
							if(jtSubJetVtxType->at(ijet).at(0)==0){
								if(jtSubJetVtxType->at(ijet).at(1)==1){
									nSingleJetsPSVtx++;
									if(jtSubJetHadronFlavor->at(ijet).at(0)==5) nSingleJetsPSVtxSingleB++;
								}
								if(jtSubJetVtxType->at(ijet).at(1)==2) nSingleJetsNoVtx++;
							}
							if(jtSubJetVtxType->at(ijet).at(0)==1){
								if(jtSubJetVtxType->at(ijet).at(1)==0){
									nSingleJetsPSVtx++;
									if(jtSubJetHadronFlavor->at(ijet).at(1)==5) nSingleJetsPSVtxSingleB++;
								}								
								if(jtSubJetVtxType->at(ijet).at(1)==1) nSingleJetsDoublePSVtx++;
								if(jtSubJetVtxType->at(ijet).at(1)==2) nSingleJetsPSNoVtx++;
							}
							if(jtSubJetVtxType->at(ijet).at(0)==2){
								if(jtSubJetVtxType->at(ijet).at(1)==0) nSingleJetsNoVtx++;
								if(jtSubJetVtxType->at(ijet).at(1)==1) nSingleJetsPSNoVtx++;
								if(jtSubJetVtxType->at(ijet).at(1)==2) nSingleJetsDoubleNoVtx++;
							}
							
							for(unsigned int isubjet=0; isubjet<jtSubJetHadronFlavor->at(ijet).size(); isubjet++){
								cout << "subjet " << isubjet << " pt: " << jtSubJetPt->at(ijet).at(isubjet) << " eta: "<< jtSubJetEta->at(ijet).at(isubjet) << " phi: "<< jtSubJetPhi->at(ijet).at(isubjet) << " flavor: "<< jtSubJetHadronFlavor->at(ijet).at(isubjet) << " nPartons: "<< jtSubJetPartonDR->at(ijet).at(isubjet).size() << " nHadrons: " << jtSubJetHadronPt->at(ijet).at(isubjet).size() << endl;
								cout << "subjetVtx info: "<< endl;
								cout << "vtx type: "<< jtSubJetVtxType->at(ijet).at(isubjet) << endl;
								for(unsigned int isubvtx=0; isubvtx<jtSubJetSvtxm->at(ijet).at(isubjet).size(); isubvtx++){
									cout << "mass: "<< jtSubJetSvtxm->at(ijet).at(isubjet).at(isubvtx) << " pt: "<< jtSubJetSvtxpt->at(ijet).at(isubjet).at(isubvtx) << " eta: "<< jtSubJetSvtxeta->at(ijet).at(isubjet).at(isubvtx) << " phi: "<< jtSubJetSvtxphi->at(ijet).at(isubjet).at(isubvtx) << " ntrk: "<< jtSubJetSvtxNtrk->at(ijet).at(isubjet).at(isubvtx) << " dl: "<< jtSubJetSvtxdl->at(ijet).at(isubjet).at(isubvtx) << " dls: "<< jtSubJetSvtxdls->at(ijet).at(isubjet).at(isubvtx) << endl;
									if(isubjet==1) subjetSVsBMass->Fill(jtSubJetSvtxm->at(ijet).at(isubjet).at(isubvtx));
								}
								int totalBHadrons=0;
								int totalCHadrons=0;
								int totalBQuarks=0;
								for(unsigned int iparton=0; iparton<jtSubJetPartonPt->at(ijet).at(isubjet).size(); iparton++){
									if(abs(jtSubJetPartonPdg->at(ijet).at(isubjet).at(iparton))==5) totalBQuarks++;
								}

								if(jtSubJetHadronPt->at(ijet).at(isubjet).size() > 1){
									for(unsigned int ihadron=0; ihadron<jtSubJetHadronPt->at(ijet).at(isubjet).size(); ihadron++){
										cout << "associated hadron pt: " << jtSubJetHadronPt->at(ijet).at(isubjet).at(ihadron) << " pdg: "<< jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron) << " eta: "<< jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron) << " phi: "<< jtSubJetHadronPhi->at(ijet).at(isubjet).at(ihadron) << endl;
										if(abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron))>500 && abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron))<600) totalBHadrons++;
										if(abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron))>400 && abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron))<500) totalCHadrons++;
										
										double drToMatch = getDR(jtSubJetEta->at(ijet).at(isubjet), jtSubJetPhi->at(ijet).at(isubjet), jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron), jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron));
										dRtoMatchedSubjet->Fill(drToMatch);
										if(isubjet){
											drToMatch = getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron), jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron));
										}
										else{
											drToMatch = getDR(jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1), jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron), jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron));
										}
										dRtoUnmatchedSubjet->Fill(drToMatch);
									}
								}
								if(jtSubJetPartonPt->at(ijet).at(isubjet).size() > 1){
									for(unsigned int iparton=0; iparton<jtSubJetPartonPt->at(ijet).at(isubjet).size(); iparton++){
										cout << "associated parton pt: " << jtSubJetPartonPt->at(ijet).at(isubjet).at(iparton) << " pdg: "<< jtSubJetPartonPdg->at(ijet).at(isubjet).at(iparton) << endl;
									}
								}
							//only fill the debug histos for the tagged b-subjets....
								if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))==5){
									if(totalBHadrons==1 && totalCHadrons>0) failureModes->Fill(4);
									else if(totalBHadrons==1 && totalBQuarks==1) failureModes->Fill(1);
									
									else if(totalBHadrons==1 && totalBQuarks>1) failureModes->Fill(2);
									else if(totalBHadrons>1) failureModes->Fill(3);
								}
							}
							
						}
					}

				}
				else{
					unknownFrac->Fill(minCSV);
					if(subjetDR>subjetDRCut) unknownFrac_dR->Fill(minCSV);
				}
			}
		}
	}
	
	double purity[40], efficiency[40];
	double purityErr[40], efficiencyErr[40];
	double purVtx[40], effVtx[40];
	
	TH1D *th1DiB = new TH1D("th1DiB","",40,0,1);
	TH1D *th1TotSub = new TH1D("th1TotSub","",40,0,1);
	TH1D *th1TotDiB = new TH1D("th1TotDiB","",40,0,1);
	for(int i=0; i<40; i++){
		th1DiB->SetBinContent(i+1, nDiBSubjets[i]);
		th1TotSub->SetBinContent(i+1, totalSubjets[i]);
		th1TotDiB->SetBinContent(i+1, totalDiBSubjets);
		
		th1DiB->SetBinError(i+1, sqrt(nDiBSubjets[i]));
		th1TotSub->SetBinError(i+1, sqrt(totalSubjets[i]));
		th1TotDiB->SetBinError(i+1, sqrt(totalDiBSubjets));
	}
	for(int ibin=0; ibin<40; ibin++){
		if(totalSubjets[ibin]){
			purity[ibin] = nDiBSubjets[ibin]/totalSubjets[ibin];
			//calculated with binomial errors
			purityErr[ibin] = getBinomialError(ibin, th1DiB, th1TotSub);
			purVtx[ibin] = nDiBVtxDRSubjets[ibin]/totalVtxDRSubjets[ibin];
		}
		else{
			purity[ibin] = 0;
			purVtx[ibin] = 0;
		}
		efficiency[ibin] = nDiBSubjets[ibin]/totalDiBSubjets;
		efficiencyErr[ibin] = getBinomialError(ibin, th1DiB, th1TotDiB);
		effVtx[ibin] = nDiBVtxDRSubjets[ibin]/totalDiBSubjets;
		cout << " ibin: "<< ibin << " nSubjets: "<< totalSubjets[ibin] << " nDiBSubjets" << nDiBSubjets[ibin] << " totalDiBSubjets" << totalDiBSubjets << endl;
	}
	TGraphErrors *purEffCurve_csv = new TGraphErrors(40, purity, efficiency, purityErr, efficiencyErr);
	purEffCurve_csv->SetNameTitle("purEffCurve_csv","purEffCurve_csv");
	
	TGraph *purEffCurve_subjetVtxDr = new TGraph(40, purVtx, effVtx);
	purEffCurve_subjetVtxDr->SetNameTitle("purEffCurve_subjetVtxDr","purEffCurve_subjetVtxDr");
	
	cout << " total singleB large CSV subjet systems: "<< nSingleJets << endl;
	cout << "1 PseudoVtx systems: "<< nSingleJetsPSVtx << " of which have " << nSingleJetsPSVtxSingleB << " systems with full reco b-jets" << endl;
	cout << "1 No-Vtx systems: " << nSingleJetsNoVtx << endl;
	cout << "Double PSVtx systems: "<< nSingleJetsDoublePSVtx << endl;
	cout << "Double No-Vtx systems: "<< nSingleJetsDoubleNoVtx << endl;
	cout << "PS+No Vtx systems: "<< nSingleJetsPSNoVtx << endl;
	
	TFile *fout = new TFile("out.root","recreate");
	fout->cd();
	
	purEffCurve_csv->Write();
	purEffCurve_subjetVtxDr->Write();
	
	lightFrac->Divide(allJets);
	lightFrac->SetFillColor(4);
	charmFrac->Divide(allJets);
	charmFrac->SetFillColor(kGreen+2);
	singleBFrac->Divide(allJets);
	singleBFrac->SetFillColor(2);
	bcFrac->Divide(allJets);
	bcFrac->SetFillColor(kCyan+1);
	doubleBFrac->Divide(allJets);
	doubleBFrac->SetFillColor(kOrange+2);
	unknownFrac->Divide(allJets);
	unknownFrac->SetFillColor(kMagenta+2);
	
	lightFrac_dR->Divide(allJets_dR);
	lightFrac_dR->SetFillColor(4);
	charmFrac_dR->Divide(allJets_dR);
	charmFrac_dR->SetFillColor(kGreen+2);
	singleBFrac_dR->Divide(allJets_dR);
	singleBFrac_dR->SetFillColor(2);
	bcFrac_dR->Divide(allJets_dR);
	bcFrac_dR->SetFillColor(kCyan+1);
	doubleBFrac_dR->Divide(allJets_dR);
	doubleBFrac_dR->SetFillColor(kOrange+2);
	unknownFrac_dR->Divide(allJets_dR);
	unknownFrac_dR->SetFillColor(kMagenta+2);
	
	lightFrac_dR_subVtxOne->Divide(allJets_dR_subVtxOne);
	lightFrac_dR_subVtxOne->SetFillColor(4);
	charmFrac_dR_subVtxOne->Divide(allJets_dR_subVtxOne);
	charmFrac_dR_subVtxOne->SetFillColor(kGreen+2);
	singleBFrac_dR_subVtxOne->Divide(allJets_dR_subVtxOne);
	singleBFrac_dR_subVtxOne->SetFillColor(2);
	bcFrac_dR_subVtxOne->Divide(allJets_dR_subVtxOne);
	bcFrac_dR_subVtxOne->SetFillColor(kCyan+1);
	doubleBFrac_dR_subVtxOne->Divide(allJets_dR_subVtxOne);
	doubleBFrac_dR_subVtxOne->SetFillColor(kOrange+2);
	unknownFrac_dR_subVtxOne->Divide(allJets_dR_subVtxOne);
	unknownFrac_dR_subVtxOne->SetFillColor(kMagenta+2);
	
	lightFrac_dR_subVtxBoth->Divide(allJets_dR_subVtxBoth);
	lightFrac_dR_subVtxBoth->SetFillColor(4);
	charmFrac_dR_subVtxBoth->Divide(allJets_dR_subVtxBoth);
	charmFrac_dR_subVtxBoth->SetFillColor(kGreen+2);
	singleBFrac_dR_subVtxBoth->Divide(allJets_dR_subVtxBoth);
	singleBFrac_dR_subVtxBoth->SetFillColor(2);
	bcFrac_dR_subVtxBoth->Divide(allJets_dR_subVtxBoth);
	bcFrac_dR_subVtxBoth->SetFillColor(kCyan+1);
	doubleBFrac_dR_subVtxBoth->Divide(allJets_dR_subVtxBoth);
	doubleBFrac_dR_subVtxBoth->SetFillColor(kOrange+2);
	unknownFrac_dR_subVtxBoth->Divide(allJets_dR_subVtxBoth);
	unknownFrac_dR_subVtxBoth->SetFillColor(kMagenta+2);
	
	THStack *ths = new THStack("flavorFrac","flavorFrac");
	ths->Add(lightFrac);
	ths->Add(charmFrac);
	ths->Add(singleBFrac);
	ths->Add(bcFrac);
	ths->Add(doubleBFrac);
	ths->Add(unknownFrac);
	ths->Write();
	
	THStack *ths_dr = new THStack("flavorFrac_drCut","flavorFrac_drCut");
	ths_dr->Add(lightFrac_dR);
	ths_dr->Add(charmFrac_dR);
	ths_dr->Add(singleBFrac_dR);
	ths_dr->Add(bcFrac_dR);
	ths_dr->Add(doubleBFrac_dR);
	ths_dr->Add(unknownFrac_dR);
	ths_dr->Write();
	
	THStack *ths_dr_vtxOne = new THStack("flavorFrac_drCut_vtxDRSingle","flavorFrac_drCut_vtxDRSingle");
	ths_dr_vtxOne->Add(lightFrac_dR_subVtxOne);
	ths_dr_vtxOne->Add(charmFrac_dR_subVtxOne);
	ths_dr_vtxOne->Add(singleBFrac_dR_subVtxOne);
	ths_dr_vtxOne->Add(bcFrac_dR_subVtxOne);
	ths_dr_vtxOne->Add(doubleBFrac_dR_subVtxOne);
	ths_dr_vtxOne->Add(unknownFrac_dR_subVtxOne);
	ths_dr_vtxOne->Write();
	
	THStack *ths_dr_vtxBoth = new THStack("flavorFrac_drCut_vtxDRBoth","flavorFrac_drCut_vtxDRBoth");
	ths_dr_vtxBoth->Add(lightFrac_dR_subVtxBoth);
	ths_dr_vtxBoth->Add(charmFrac_dR_subVtxBoth);
	ths_dr_vtxBoth->Add(singleBFrac_dR_subVtxBoth);
	ths_dr_vtxBoth->Add(bcFrac_dR_subVtxBoth);
	ths_dr_vtxBoth->Add(doubleBFrac_dR_subVtxBoth);
	ths_dr_vtxBoth->Add(unknownFrac_dR_subVtxBoth);
	ths_dr_vtxBoth->Write();
	
	dRCSVCorrelation_doubleB->Write();
	dRCSVCorrelation_singleB->Write();
	dRCSVCorrelation_doubleB_doubleType0->Write();
	dRCSVCorrelation_doubleB_singleType0->Write();
	dRCSVCorrelation_doubleB_noType0->Write();
	
	csv12Correlation->Write();
	csv12Correlation_doubleB->Write();
	csv12Correlation_singleB->Write();
	
	subjetVtxToDR->Write();
	subjetVtxToDR_csv09->Write();
	subjetVtxToDR_bFlav->Write();
	
	csvDistr->Write();
	csvDistr_bsubjet->Write();
	
	jetSVMass_dB->Write();
	jetSVMass_sB->Write();
	
	subjetSVsBMass->Write();
	subjetSVdBMass->Write();
	
	dRtoMatchedSubjet->Write();
	dRtoUnmatchedSubjet->Write();
	
	singleDoubleNVtx->Write();
	
	failureModes->Write();
	
	statusType->Write();
	statusJet->Write();
	doublebStatusCorr->Write();
	fout->Close();
}
