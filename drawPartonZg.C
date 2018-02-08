

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "TChain.h"

using namespace std;

double xPtBins[] = {120,160,250,9999};

double getdR(double eta1, double eta2, double phi1, double phi2){
	
	double dr = sqrt(pow(eta1-eta2,2)+pow(TMath::ACos(TMath::Cos(phi1-phi2)),2));
	return dr;
}

double getzg(double j1, double j2){
	return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}

void drawPartonZg(){
	
	TH1D *gspZg = new TH1D("gspZg","",10,0,0.5); gspZg->Sumw2();
	TH1D *gspdR = new TH1D("gspdR","",80,0,3.14); gspdR->Sumw2();
	TH1D *gspdPhi = new TH1D("gspdPhi","",80,0,3.14); gspdPhi->Sumw2();
	TH1D *subjetZg[4];
	TH1D *matchedPartonZg[4];
	for(int i=0; i<4; i++){
		subjetZg[i] = new TH1D(Form("subjetZg_pt%d",i),"",10,0,0.5); subjetZg[i]->Sumw2();
		matchedPartonZg[i] = new TH1D(Form("matchedPartonZg_%d",i),"",10,0,0.5); matchedPartonZg[i]->Sumw2();
	}
	TH1D *subjetdR = new TH1D("subjetdR","",30,0,0.5); subjetdR->Sumw2();
	TH1D *matchedPartondR = new TH1D("matchedPartondR","",80,0,3.14); matchedPartondR->Sumw2();
	TH1D *matchedTruthPartonZg = new TH1D("matchedTruthPartonZg","",10,0,0.5); matchedTruthPartonZg->Sumw2();
	TH1D *matchedTruthPartondR = new TH1D("matchedTruthPartondR","",80,0,3.14); matchedTruthPartondR->Sumw2();
	TH1D *matchedTagPartonZg = new TH1D("matchedTagPartonZg","",10,0,0.5); matchedTagPartonZg->Sumw2();
	TH1D *matchedTagPartondR = new TH1D("matchedTagPartondR","",80,0,3.14); matchedTagPartondR->Sumw2();
	
	TH1D *partonSubjetdR = new TH1D("partonSubjetdR","",50,0,0.5); partonSubjetdR->Sumw2();
	TH1D *recoHadronSubjetdR = new TH1D("recoHadronSubjetdR","",50,0,0.5); recoHadronSubjetdR->Sumw2();
	
	TChain *tin = new TChain("HiGenParticleAna/hi");
	TChain *jetch = new TChain("akSoftDrop4PFJetAnalyzer/t");
	
	vector<float> *pt=0, *eta=0, *phi=0;
	vector<int> *pdg=0, *nMothers=0, *nDaughters=0;
	vector<vector<int> > *daughterIdx=0;
	tin->SetBranchAddress("pt",&pt);
	tin->SetBranchAddress("eta",&eta);
	tin->SetBranchAddress("phi",&phi);
	tin->SetBranchAddress("pdg",&pdg);
	tin->SetBranchAddress("nMothers",&nMothers);
	tin->SetBranchAddress("nDaughters",&nDaughters);
	tin->SetBranchAddress("daughterIdx",&daughterIdx);
	
	vector<vector<float> > *jtSubJetPt=0, *jtSubJetEta=0, *jtSubJetPhi=0, *jtSubJetcsvV2=0;
	vector<vector<float> > *jtSubJetHadronFlavor=0;
	vector<vector<vector<float> > > *jtSubJetHadronEta=0, *jtSubJetHadronPhi=0, *jtSubJetHadronPdg=0, *jtSubJetSvtxm=0;
	int nref;
	float jtpt[1000];
	jetch->SetBranchAddress("nref",&nref);
	jetch->SetBranchAddress("jtpt",jtpt);
	jetch->SetBranchAddress("jtSubJetPt",&jtSubJetPt);
	jetch->SetBranchAddress("jtSubJetEta",&jtSubJetEta);
	jetch->SetBranchAddress("jtSubJetPhi",&jtSubJetPhi);
	jetch->SetBranchAddress("jtSubJetHadronFlavor",&jtSubJetHadronFlavor);
	jetch->SetBranchAddress("jtSubJetHadronPdg",&jtSubJetHadronPdg);
	jetch->SetBranchAddress("jtSubJetHadronEta",&jtSubJetHadronEta);
	jetch->SetBranchAddress("jtSubJetHadronPhi",&jtSubJetHadronPhi);
	jetch->SetBranchAddress("jtSubJetcsvV2",&jtSubJetcsvV2);
	jetch->SetBranchAddress("jtSubJetSvtxm",&jtSubJetSvtxm);

	ifstream filelist("ppMC_genParticleList.txt");
	std::string filename;
	while(filelist>>filename){
		tin->Add(filename.c_str());
		jetch->Add(filename.c_str());
	}		
	int nentries = tin->GetEntries();
	cout << "nentries: "<< nentries << endl;
	for(int i=0; i<nentries-1; i++){
		if(i>0 && i%10000==0) cout << "at entry " << i << " of " << nentries << endl;
		tin->GetEntry(i);
		jetch->GetEntry(i);
		
		for(unsigned int iparton=0; iparton<pt->size(); iparton++){
			
			if(abs(pdg->at(iparton))==21 && nDaughters->at(iparton)>=2 && pt->at(iparton)>50){
				
				//cout << "found candidate, idx: "<< iparton << " pdg: " << pdg->at(iparton) << " pt: "<< pt->at(iparton) << " nDaughters: "<< nDaughters->at(iparton) << endl;
				//if(daughterIdx->at(iparton).at(0)<0 || daughterIdx->at(iparton).at(1)<0) continue;
				//cout << "found candidate, dindx1: "<< daughterIdx->at(iparton).at(0) << endl;
				//cout << "found candidate, dindx2: "<< daughterIdx->at(iparton).at(1) << endl;
				
				int daughtPdgArr[100];
				int daughtCount=0;
				int idummy=0;
				for(int idaught=0; idaught<nDaughters->at(iparton); idaught++){
					if(abs(pdg->at(daughterIdx->at(iparton).at(idaught)))==5){
						daughtCount++;
						daughtPdgArr[idummy++]=idaught;
					}
				}
				
				if(daughtCount>=2){
					
					double pt1 = pt->at(daughterIdx->at(iparton).at(daughtPdgArr[0]));
					double pt2 = pt->at(daughterIdx->at(iparton).at(daughtPdgArr[1]));
					//if((pt1+pt2)/pt->at(iparton) < 0.85) continue;
					
					//cout << "nDaughters: "<< nDaughters->at(iparton) << endl;
					//for(int idaught=0; idaught<daughtCount; idaught++){
					//cout << "parton " << idaught << " pt, eta, phi, pdg: "<< pt->at(daughterIdx->at(iparton).at(daughtPdgArr[idaught])) << " " << eta->at(daughterIdx->at(iparton).at(daughtPdgArr[idaught])) << " " << phi->at(daughterIdx->at(iparton).at(daughtPdgArr[idaught])) << " " << pdg->at(daughterIdx->at(iparton).at(daughtPdgArr[idaught])) << endl;
					//}
					
					double eta1 = eta->at(daughterIdx->at(iparton).at(daughtPdgArr[0]));
					double eta2 = eta->at(daughterIdx->at(iparton).at(daughtPdgArr[1]));
					double phi1 = phi->at(daughterIdx->at(iparton).at(daughtPdgArr[0]));
					double phi2 = phi->at(daughterIdx->at(iparton).at(daughtPdgArr[1]));
					
					//do subjet matching now...
					int jetid=-99;
					int subjet1=-99;
					int subjet2=-99;
					double drDist1=-99;
					double drDist2=-99;
					for(int ijet=0; ijet<nref; ijet++){
						if(jtpt[ijet]<100) continue;
						double minDist1=3.5;
						double minDist2=3.5;
						for(unsigned int isubjet=0; isubjet<jtSubJetPt->at(ijet).size(); isubjet++){
							
							double dist1 = getdR(eta1,jtSubJetEta->at(ijet).at(isubjet),phi1,jtSubJetPhi->at(ijet).at(isubjet));
							double dist2 = getdR(eta2,jtSubJetEta->at(ijet).at(isubjet),phi2,jtSubJetPhi->at(ijet).at(isubjet));
							if(dist1<dist2){
								if(dist1<minDist1){
									subjet1=isubjet;
									minDist1=dist1;
								}
							}
							else if(dist2<minDist2){
								subjet2=isubjet;
								minDist2=dist2;
								//cout << "parton eta1,phi1, eta2,phi2: "<< eta1 << "," << phi1 << "  , " << eta2 << "," << phi2 << endl;
							}
						}
						if(subjet1>=0 && subjet2>=0){
							jetid=ijet;
							drDist1 = minDist1;
							drDist2 = minDist2;
							break;
						}
						else{
							subjet1=-99;
							subjet2=-99;
						}
					
					}
					if(subjet1>=0 && subjet2>=0){
						
						int ptBin = 0;
						cout << "jtpt: "<< jtpt[jetid] << endl;
						while(jtpt[jetid]>xPtBins[ptBin]){
							ptBin++;
							cout << "jet pt larger than " << xPtBins[ptBin] << "..." << endl;
						}
						cout << "ptbin: "<< ptBin << endl;
						//cout << "subjet zg: " << getzg(jtSubJetPt->at(jetid).at(subjet1), jtSubJetPt->at(jetid).at(subjet2)) << endl;
						subjetZg[ptBin]->Fill(getzg(jtSubJetPt->at(jetid).at(subjet1), jtSubJetPt->at(jetid).at(subjet2)));
						double subEta1 = jtSubJetEta->at(jetid).at(subjet1);
						double subEta2 = jtSubJetEta->at(jetid).at(subjet2);
						double subPhi1 = jtSubJetPhi->at(jetid).at(subjet1);
						double subPhi2 = jtSubJetPhi->at(jetid).at(subjet2);
						subjetdR->Fill(getdR(subEta1,subEta2,subPhi1,subPhi2));
						matchedPartonZg[ptBin]->Fill(getzg(pt1,pt2));
						matchedPartondR->Fill(getdR(eta1,eta2,phi1,phi2));
						partonSubjetdR->Fill(drDist1);
						partonSubjetdR->Fill(drDist2);
						if(jtSubJetcsvV2->at(jetid).at(subjet1)>0.8 && jtSubJetcsvV2->at(jetid).at(subjet2)>0.8){
							if(jtSubJetSvtxm->at(jetid).at(subjet1).size() && jtSubJetSvtxm->at(jetid).at(subjet2).size()){
								if(jtSubJetSvtxm->at(jetid).at(subjet1).at(0)>1.5 && jtSubJetSvtxm->at(jetid).at(subjet2).at(0)>1.5){
									//cout << "subjets are tagged!" << endl;
									matchedTagPartonZg->Fill(getzg(pt1,pt2));
									matchedTagPartondR->Fill(getdR(eta1,eta2,phi1,phi2));
								}
							}
						}
						if(jtSubJetHadronFlavor->at(jetid).at(subjet1)==5 && jtSubJetHadronFlavor->at(jetid).at(subjet2)==5){
							//cout << "subjets both have associated b-hadrons!" << endl;
							matchedTruthPartonZg->Fill(getzg(pt1,pt2));
							matchedTruthPartondR->Fill(getdR(eta1,eta2,phi1,phi2));
						}
						
					}
					
					double zg = getzg(pt1,pt2);
					if(zg<0.1) continue;
					gspZg->Fill(zg);
					gspdR->Fill(sqrt(pow(eta1-eta2,2)+pow(phi1-phi2,2)));
					gspdPhi->Fill(TMath::ACos(TMath::Cos(phi1-phi2)));
				}
			}
		}
		
		for(int ijet=0; ijet<nref; ijet++){
			if(jtpt[ijet]<60) continue;
			for(unsigned int isubjet=0; isubjet<jtSubJetHadronFlavor->at(ijet).size(); isubjet++){
				if(jtSubJetHadronFlavor->at(ijet).at(isubjet)==5){
					for(unsigned int ihadron=0; ihadron<jtSubJetHadronPdg->at(ijet).at(isubjet).size(); ihadron++){
						if(abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron))>500 && abs(jtSubJetHadronPdg->at(ijet).at(isubjet).at(ihadron))<600){
							double hadronSubjetdR = getdR(jtSubJetHadronEta->at(ijet).at(isubjet).at(ihadron),jtSubJetEta->at(ijet).at(isubjet),jtSubJetHadronPhi->at(ijet).at(isubjet).at(ihadron),jtSubJetPhi->at(ijet).at(isubjet));
							recoHadronSubjetdR->Fill(hadronSubjetdR);
						}
					}
				}
			}
		}
		
	}
	
	//gspZg->Scale(1./gspZg->Integral());
	//gspZg->Scale(1./gspZg->GetBinWidth(2));
	gspZg->SetXTitle("parton-level z_{G}");
	gspZg->SetYTitle("1/N dN/dz_{G}");
	//gspZg->Draw();
	
	TFile *fout = new TFile("partonCheck_mixedEvt.root","recreate");
	fout->cd();
	gspZg->Write();
	//gspdR->Scale(1./gspdR->Integral());
	gspdR->Write();
	//gspdPhi->Scale(1./gspdPhi->Integral());
	gspdPhi->Write();
	
	for(int i=0; i<4; i++){
		subjetZg[i]->Write();
		matchedPartonZg[i]->Write();
	}
	subjetdR->Write();
	matchedPartondR->Write();
	matchedTruthPartonZg->Write();
	matchedTruthPartondR->Write();
	matchedTagPartonZg->Write();
	matchedTagPartondR->Write();
	partonSubjetdR->Write();
	recoHadronSubjetdR->Write();
	
	fout->Close();
	
	
}