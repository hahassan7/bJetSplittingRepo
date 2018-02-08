
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"

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
//#else
//template class std::vector<std::vector<float> >;
//template class std::vector<std::vector<std::vector<float> > >;
#endif

const int MAXJETS=70;
const double PI=3.1415926;

double getDR(double eta1, double phi1, double eta2, double phi2){
	return(sqrt(pow(eta1-eta2,2)+pow(phi1-phi2,2)));
}

double getzg(float j1, float j2){
	return (j1>j2 ? (j2/(j1+j2)):(j1/(j1+j2)));
}

void testSubjetPartonAssn(){
	
	double jetDR = 0.4;

	TH1D *bbInsideJet = new TH1D("bbInsideJet","",3,0,3);
	TH1D *bbInsideJetFlavor5 = new TH1D("bbInsideJetFlavor5","",3,0,3);
	TH1D *subjetMatch = new TH1D("subjetMatch","",4,0,4);
	TH1D *subjetMatchFlavor5 = new TH1D("subjetMatchFlavor5","",4,0,4);
	
	TH1D *gspDR = new TH1D("gspDR","",50,0,2*PI); gspDR->Sumw2();
	TH1D *gspxj = new TH1D("gspxj","",50,0,1); gspxj->Sumw2();
	
	TH1D *partonSubjetdR = new TH1D("partonSubjetdR","",100,0,2*PI); partonSubjetdR->Sumw2();
	TH1D *dRforbbbar = new TH1D("dRforbbbar","",20,0,0.5); dRforbbbar->Sumw2();
	TH1D *dRforSingleb = new TH1D("dRforSingleb","",20,0,0.5); dRforSingleb->Sumw2();
	TH1D *dRforIncl = new TH1D("dRforIncl","",20,0,0.5); dRforIncl->Sumw2();
	TH1D *dRforNonbbbar = new TH1D("dRforNonbbbar","",20,0,0.5); dRforNonbbbar->Sumw2();
	
	TH1D *zgforbbbar = new TH1D("zgforbbbar","",20,0,0.5); zgforbbbar->Sumw2();
	TH1D *zgforSingleb = new TH1D("zgforSingleb","",20,0,0.5); zgforSingleb->Sumw2();
	TH1D *zgforIncl = new TH1D("zgforIncl","",20,0,0.5); zgforIncl->Sumw2();
	TH1D *zgforNonbbbar = new TH1D("zgforNonbbbar","",20,0,0.5); zgforNonbbbar->Sumw2();
	
	//***************** START THE NEW GHOSTING HISTOS ***************//
	
	TH1D *subjetHadronDR = new TH1D("subjetHadronDR","",30,0.,0.5); subjetHadronDR->Sumw2();
	TH1D *bsubjetHadronDR = new TH1D("bsubjetHadronDR","",30,0.,0.5); bsubjetHadronDR->Sumw2();
	TH1D *subjetPartonDR = new TH1D("subjetPartonDR","",30,0.,0.5); subjetPartonDR->Sumw2();
	TH1D *bsubjetPartonDR = new TH1D("bsubjetPartonDR","",30,0.,0.5); bsubjetPartonDR->Sumw2();
	
	TH1D *subjetHadronFlavor = new TH1D("subjetHadronFlavor","",7,0,7); subjetHadronFlavor->Sumw2();
	TH1D *subjetPartonFlavor = new TH1D("subjetPartonFlavor","",7,0,7); subjetPartonFlavor->Sumw2();
	
	TH2D *subjetZg2DR = new TH2D("subjetZg2DR","",100,0,0.5,100,0,0.5); subjetZg2DR->Sumw2();
	TH2D *singlebSubjetZg2DR = new TH2D("singlebSubjetZg2DR","",100,0,0.5,100,0,0.5); singlebSubjetZg2DR->Sumw2();
	TH2D *doublebSubjetZg2DR = new TH2D("doublebSubjetZg2DR","",100,0,0.5,100,0,0.5); doublebSubjetZg2DR->Sumw2();
	
	TH2D *subjetZg2DR_gsp = new TH2D("subjetZg2DR_gsp","",100,0,0.5,100,0,0.5); subjetZg2DR_gsp->Sumw2();
	TH2D *singlebSubjetZg2DR_gsp = new TH2D("singlebSubjetZg2DR_gsp","",100,0,0.5,100,0,0.5); singlebSubjetZg2DR_gsp->Sumw2();
	TH2D *doublebSubjetZg2DR_gsp = new TH2D("doublebSubjetZg2DR_gsp","",100,0,0.5,100,0,0.5); doublebSubjetZg2DR_gsp->Sumw2();
	
	double nSingleB[40], nSingleB_gsp[40];
	double nDoubleB[40], nDoubleB_gsp[40];
	double nNonB[40], nNonB_gsp[40];
	double testDoubleBGSP_csv05[40];
	double testDBGSPden[40];
	double testnonBGSP_csv05[40];
	double testnonBGSPden[40];
	for(int i=0; i<40; i++){
		nSingleB[i] = 0;;
		nDoubleB[i] = 0;
		nNonB[i] = 0;
		
		nSingleB_gsp[i] = 0;;
		nDoubleB_gsp[i] = 0;
		nNonB_gsp[i] = 0;
	}
	for(int i=0; i<40; i++){
		testDoubleBGSP_csv05[i]=0;
		testDBGSPden[i]=0;
		testnonBGSP_csv05[i]=0;
		testnonBGSPden[i]=0;
	}
	
	//****************************************************************//
	
	int nref;
	float jtpt[MAXJETS], jteta[MAXJETS], jtphi[MAXJETS];
	int refparton_flavorForB[MAXJETS];
	vector<vector<float> > *jtSubJetPt=0, *jtSubJetEta=0, *jtSubJetPhi=0;
	vector<vector<float> > *jtSubJetPartonFlavor=0, *jtSubJetHadronFlavor=0;
	vector<vector<float> > *jtSubJetCSV=0;
	int refparton_flavorProcess[MAXJETS];
	
	vector<vector<vector<float> > > *jtSubJetHadronDR=0, *jtSubJetPartonDR=0;
	vector<vector<vector<float> > > *jtSubJetHadronPt=0, *jtSubJetHadronEta=0, *jtSubJetHadronPhi=0;
	
	vector<float> *pt=0, *eta=0, *phi=0;
	vector<int> *pdg=0, *status=0;
	
	///*************************************************************//

	/*string filelist = "filelist.txt";
	std::ifstream instr("filelist.txt", std::ifstream::in);*/
	std::string filename = "bjet120Forest.root";
	TFile *fin;
	for(int ifile=0; ifile<1; ifile++){
		
		cout << "top of file loop" << endl;
		//instr>>filename;
		cout << "infile loads correctly" << endl;
		fin = TFile::Open(filename.c_str());
		cout << "starting file " << filename << "... " << endl;
		
		TTree *subjetTree = (TTree*)fin->Get("akSoftDrop4PFJetAnalyzer/t");
		TTree *partonTree = (TTree*)fin->Get("HiGenParticleAna/hi");
		
		subjetTree->SetBranchAddress("nref", &nref);
		subjetTree->SetBranchAddress("jtpt", jtpt);
		subjetTree->SetBranchAddress("jteta", jteta);
		subjetTree->SetBranchAddress("jtphi", jtphi);
		subjetTree->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
		subjetTree->SetBranchAddress("refparton_flavorProcess",refparton_flavorProcess);
		
		subjetTree->SetBranchAddress("jtSubJetPt", &jtSubJetPt);
		subjetTree->SetBranchAddress("jtSubJetEta", &jtSubJetEta);
		subjetTree->SetBranchAddress("jtSubJetPhi", &jtSubJetPhi);
		subjetTree->SetBranchAddress("jtSubJetcsvV2",&jtSubJetCSV);
		subjetTree->SetBranchAddress("jtSubJetPartonFlavor", &jtSubJetPartonFlavor);
		subjetTree->SetBranchAddress("jtSubJetHadronFlavor", &jtSubJetHadronFlavor);
		subjetTree->SetBranchAddress("jtSubJetHadronDR",&jtSubJetHadronDR);
		subjetTree->SetBranchAddress("jtSubJetPartonDR",&jtSubJetPartonDR);
		subjetTree->SetBranchAddress("jtSubJetHadronPt",&jtSubJetHadronPt);
		subjetTree->SetBranchAddress("jtSubJetHadronEta",&jtSubJetHadronEta);
		subjetTree->SetBranchAddress("jtSubJetHadronPhi",&jtSubJetHadronPhi);
		
		partonTree->SetBranchAddress("pt",&pt);
		partonTree->SetBranchAddress("eta",&eta);
		partonTree->SetBranchAddress("phi",&phi);
		partonTree->SetBranchAddress("pdg",&pdg);
		//partonTree->SetBranchAddress("sta",&status);
		
		for(int ievt=0; ievt<subjetTree->GetEntries(); ievt++){
			
			if(ievt%10000==0) cout << "at entry " << ievt << " of " << subjetTree->GetEntries() << endl;
			subjetTree->GetEntry(ievt);
			partonTree->GetEntry(ievt);
			
		//cout << "at event " << ievt << endl;
			
			int bIdx1=-1;
			int bIdx2=-1;
		//find first bbbar splitting
			for(unsigned int ipart=1; ipart<pt->size(); ipart++){
			//if found two adjacent b-quarks (from same parent)
				if(abs(pdg->at(ipart))==5 && abs(pdg->at(ipart-1))==5){
					bIdx1=ipart-1; bIdx2=ipart;
					break;
				}
			}
			if(bIdx1!=-1 && bIdx2!=-1){
				
				gspDR->Fill(getDR(eta->at(bIdx1), phi->at(bIdx1), eta->at(bIdx2), phi->at(bIdx2)));
				double pt1 = pt->at(bIdx1);
				double pt2 = pt->at(bIdx2);
				gspxj->Fill( pt1>pt2 ? pt2/pt1: pt1/pt2);
			}
			
			if(nref>MAXJETS) cout << "OH NO!!!" << endl;
			for(int ijet=0; ijet<nref; ijet++){
				
				bool doublebJet=false;
				
				if(jtpt[ijet]<50 || abs(jteta[ijet])>2.0) continue;
				
				double drp1=0, drp2=0;
				if(bIdx1!=-1 && bIdx2!=-1){
					double drp1 = getDR(jteta[ijet], jtphi[ijet], eta->at(bIdx1), phi->at(bIdx1));
					double drp2 = getDR(jteta[ijet], jtphi[ijet], eta->at(bIdx2), phi->at(bIdx2));
				}
				if(drp1<jetDR && drp2<jetDR){
					bbInsideJet->Fill(2);
					if(abs(refparton_flavorForB[ijet])==5) bbInsideJetFlavor5->Fill(2);
					doublebJet = true;
				}
				else if(drp1<jetDR || drp2<jetDR){
					bbInsideJet->Fill(1);
					if(abs(refparton_flavorForB[ijet])==5) bbInsideJetFlavor5->Fill(1);
				}
				else{
					bbInsideJet->Fill(0);
					if(abs(refparton_flavorForB[ijet])==5) bbInsideJetFlavor5->Fill(0);
				}
				
			//check out the subjets...
				if(jtSubJetPt->at(ijet).size()<2) continue;
				
			//******** Start the ghosting info... ******//
				double csvArr1 = (40.*jtSubJetCSV->at(ijet).at(0));
				double csvArr2 = (40.*jtSubJetCSV->at(ijet).at(1));
				int avgCSV = (csvArr2+csvArr1)/2.;
				if(avgCSV<0) avgCSV=0;
				
				double subjetzg = getzg(jtSubJetPt->at(ijet).at(0), jtSubJetPt->at(ijet).at(1));
				double subjetDR = getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1));
				
				//if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 && abs(jtSubJetHadronFlavor->at(ijet).at(1))==5){
				for(int icsv=0; icsv<40; icsv++){
					if((refparton_flavorProcess[ijet]==5 || refparton_flavorProcess[ijet]==6) && avgCSV > icsv) testnonBGSP_csv05[icsv]++;
					if(avgCSV > icsv) testnonBGSPden[icsv]+=1;
					
					if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 && abs(jtSubJetHadronFlavor->at(ijet).at(1))==5){
						if((refparton_flavorProcess[ijet]==5 || refparton_flavorProcess[ijet]==6) && avgCSV > icsv) testDoubleBGSP_csv05[icsv]++;
						if(avgCSV > icsv) testDBGSPden[icsv]+=1;
					}
				}
				
				//}

					if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 && abs(jtSubJetHadronFlavor->at(ijet).at(1))==5){
						doublebSubjetZg2DR->Fill(subjetzg, subjetDR);
						for(int icsv=avgCSV; icsv<40; icsv++){ nDoubleB[icsv]++; }
					}
				      else if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 || abs(jtSubJetHadronFlavor->at(ijet).at(1))==5){
				      	singlebSubjetZg2DR->Fill(subjetzg, subjetDR);
				      	for(int icsv=avgCSV; icsv<40; icsv++){ nSingleB[icsv]++; }
				      }
				      else{
				      	subjetZg2DR->Fill(subjetzg, subjetDR);
				      	for(int icsv=avgCSV; icsv<40; icsv++){ nNonB[icsv]++; }
				      }
			
			      if((refparton_flavorProcess[ijet]==5 || refparton_flavorProcess[ijet]==6)){
			      	if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 && abs(jtSubJetHadronFlavor->at(ijet).at(1))==5){
			      		doublebSubjetZg2DR_gsp->Fill(subjetzg, subjetDR);
			      		for(int icsv=avgCSV; icsv<40; icsv++){ nDoubleB_gsp[icsv]++; }
			      	}
			            else if(abs(jtSubJetHadronFlavor->at(ijet).at(0))==5 || abs(jtSubJetHadronFlavor->at(ijet).at(1))==5){
			            	singlebSubjetZg2DR_gsp->Fill(subjetzg, subjetDR);
			            	for(int icsv=avgCSV; icsv<40; icsv++){ nSingleB_gsp[icsv]++; }
			            }
     			            else{
     			            	subjetZg2DR_gsp->Fill(subjetzg, subjetDR);
     			            	for(int icsv=avgCSV; icsv<40; icsv++){ nNonB_gsp[icsv]++; }
     			            }
     			      }

				
				for(int isubjet=0; isubjet<2; isubjet++){
					for(unsigned int ihadron=0; ihadron<jtSubJetHadronDR->at(ijet).at(isubjet).size(); ihadron++){
						subjetHadronDR->Fill(jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron));
						if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))==5) bsubjetHadronDR->Fill(jtSubJetHadronDR->at(ijet).at(isubjet).at(ihadron));
					}
					for(unsigned int iparton=0; iparton<jtSubJetPartonDR->at(ijet).at(isubjet).size(); iparton++){
						subjetPartonDR->Fill(jtSubJetPartonDR->at(ijet).at(isubjet).at(iparton));
						if(abs(jtSubJetPartonFlavor->at(ijet).at(isubjet))==5) bsubjetPartonDR->Fill(jtSubJetPartonDR->at(ijet).at(isubjet).at(iparton));
					}
					if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))<6) subjetHadronFlavor->Fill(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet)));
					else if(abs(jtSubJetHadronFlavor->at(ijet).at(isubjet))==21) subjetHadronFlavor->Fill(6);
					
					if(abs(jtSubJetPartonFlavor->at(ijet).at(isubjet))<6) subjetPartonFlavor->Fill(abs(jtSubJetPartonFlavor->at(ijet).at(isubjet)));
					else if(abs(jtSubJetPartonFlavor->at(ijet).at(isubjet))==21) subjetPartonFlavor->Fill(6);
				}
				
				
			//********** finish ghosting info ***********//
				
				if(bIdx1!=-1 && bIdx2!=-1){
					
					dRforIncl->Fill(getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1)));
					zgforIncl->Fill(getzg(jtSubJetPt->at(ijet).at(0), jtSubJetPt->at(ijet).at(1)) );
					
					double closestSubjet_p1=999;
					double closestSubjet_p2=999;
					int clSubjetIdx1=-1, clSubjetIdx2=-1;
					for(unsigned int isubjet=0; isubjet<jtSubJetPt->at(ijet).size(); isubjet++){
						
				//cout << "at subjet " << isubjet << endl;
						
						double drsubp1 = getDR(jtSubJetEta->at(ijet).at(isubjet), jtSubJetPhi->at(ijet).at(isubjet), eta->at(bIdx1), phi->at(bIdx1));
						double drsubp2 = getDR(jtSubJetEta->at(ijet).at(isubjet), jtSubJetPhi->at(ijet).at(isubjet), eta->at(bIdx2), phi->at(bIdx2));
						
						if(drsubp1<closestSubjet_p1){
							closestSubjet_p1 = drsubp1;
							clSubjetIdx1=isubjet;
						}
						
						if(drsubp2<closestSubjet_p2){
							closestSubjet_p2 = drsubp2;
							clSubjetIdx2=isubjet;
						}
					}
					
					bool subjetMatchFlag=false;
					if(clSubjetIdx1!=-1){
						if(clSubjetIdx1!=clSubjetIdx2){
							partonSubjetdR->Fill(closestSubjet_p1);
							partonSubjetdR->Fill(closestSubjet_p2);
						}
						
						if(clSubjetIdx1!=clSubjetIdx2 && closestSubjet_p1<0.4 && closestSubjet_p2<0.4){
							subjetMatch->Fill(3);
							if(abs(refparton_flavorForB[ijet])==5) subjetMatchFlavor5->Fill(3);
							subjetMatchFlag=true;
							dRforbbbar->Fill( getDR(jtSubJetEta->at(ijet).at(clSubjetIdx1), jtSubJetPhi->at(ijet).at(clSubjetIdx1), jtSubJetEta->at(ijet).at(clSubjetIdx2), jtSubJetPhi->at(ijet).at(clSubjetIdx2)));
							zgforbbbar->Fill( getzg(jtSubJetPt->at(ijet).at(clSubjetIdx1),jtSubJetPt->at(ijet).at(clSubjetIdx2)));
						}
						else if(clSubjetIdx1!=clSubjetIdx2 && (closestSubjet_p1<0.4 || closestSubjet_p2<0.4) ) {
							subjetMatch->Fill(2);
							if(abs(refparton_flavorForB[ijet])==5) subjetMatchFlavor5->Fill(2);
							dRforSingleb->Fill( getDR(jtSubJetEta->at(ijet).at(clSubjetIdx1), jtSubJetPhi->at(ijet).at(clSubjetIdx1), jtSubJetEta->at(ijet).at(clSubjetIdx2), jtSubJetPhi->at(ijet).at(clSubjetIdx2)));
							zgforSingleb->Fill( getzg(jtSubJetPt->at(ijet).at(clSubjetIdx1),jtSubJetPt->at(ijet).at(clSubjetIdx2)));
						}
						else if(clSubjetIdx1!=clSubjetIdx2){
							subjetMatch->Fill(1);
							if(abs(refparton_flavorForB[ijet])==5) subjetMatchFlavor5->Fill(1);
						}
						else{
							cout << "subjet dr1: "<< closestSubjet_p1 << " dr2: "<< closestSubjet_p2 << endl;
							cout << " subjet idx1: "<< clSubjetIdx1 << " idx2: "<< clSubjetIdx2 << endl;
							subjetMatch->Fill(0);
							if(abs(refparton_flavorForB[ijet])==5) subjetMatchFlavor5->Fill(0);
							dRforNonbbbar->Fill( getDR(jtSubJetEta->at(ijet).at(0), jtSubJetPhi->at(ijet).at(0), jtSubJetEta->at(ijet).at(1), jtSubJetPhi->at(ijet).at(1)));
							zgforNonbbbar->Fill( getzg(jtSubJetPt->at(ijet).at(0), jtSubJetPt->at(ijet).at(1)) );
						}
					}
				}
			}
		}
		cout << "closing file " << filename << endl;
		fin->Close();
		cout << "closed " << endl;
	}
	
	cout << " numbers: "<< endl;
	cout << testDoubleBGSP_csv05[0] << " " << testDBGSPden[0] << endl;
	cout << testDoubleBGSP_csv05[1] << " " << testDBGSPden[1] << endl;
	cout << testDoubleBGSP_csv05[2] << " " << testDBGSPden[2] << endl;
	
	for(int icsv=0; icsv<40; icsv++){ 
		if(testDBGSPden[icsv]) testDoubleBGSP_csv05[icsv]/=testDBGSPden[icsv]; 
		else testDoubleBGSP_csv05[icsv]=0; 
		
		if(testnonBGSPden[icsv]) testnonBGSP_csv05[icsv]/=testnonBGSPden[icsv]; 
		else testnonBGSP_csv05[icsv]=0;
	}
	cout << " test CSV > 0.2 double-b-tagged GSP efficiency: "<< testDoubleBGSP_csv05[8] << " CSV > 0.5: "<< testDoubleBGSP_csv05[20] << " CSV > 0.9: "<< testDoubleBGSP_csv05[36] << endl;
	
	double doubleTagEff[40];
	double singleTagEff[40];
	double singleTagRej[40];
	double untagEff[40];
	
	double doubleTagEff_gsp[40];
	double singleTagEff_gsp[40];
	double AllRej_gsp[40];
	double untagEff_gsp[40];
	double SingleRej_gsp[40];
	double DoubleRej_gsp[40];
	
	double bareCSV[40];
	
	cout << "starting final arrays" << endl;
	for(int i=0; i<40; i++){
		bareCSV[i] = (double)i/40.;
		
		doubleTagEff[i] = 1.-(nDoubleB[i]/nDoubleB[0]);
		singleTagEff[i] = 1.-(nSingleB[i]/nSingleB[0]);
		untagEff[i] = 1.-(nNonB[i]/(nNonB[0]));
		singleTagRej[i] = nSingleB[i]/nSingleB[0];
		
		doubleTagEff_gsp[i] = nDoubleB_gsp[i]/nDoubleB[i];
		singleTagEff_gsp[i] = nSingleB_gsp[i]/nSingleB[i];
		untagEff_gsp[i] = (nDoubleB_gsp[i]+nSingleB_gsp[i]+nNonB_gsp[i])/(nDoubleB[i]+nSingleB[i]+nNonB[i]);
		double untagJets = (nDoubleB[i]+nSingleB[i]+nNonB[i]) - (nDoubleB_gsp[i]+nSingleB_gsp[i]+nNonB_gsp[i]);
		AllRej_gsp[i] = (untagJets/(nDoubleB[i]+nSingleB[i]+nNonB[i]));

		SingleRej_gsp[i] = ((nSingleB[i]-nSingleB_gsp[i])/nSingleB[i]);
		DoubleRej_gsp[i] = ((nDoubleB[i]-nDoubleB_gsp[i])/nDoubleB[i]);
		
		//protect against statistical fluctuations that make the untag eff go nuts
		if(untagEff[i]<0) untagEff[i]=0.;
		if(untagEff[i]>1) untagEff[i]=1.;
		if(AllRej_gsp[i]<0) AllRej_gsp[i]=1.;
		
		cout << "csv value > " << (double)i/40. << " has single-tag eff of " << SingleRej_gsp[i] << " double tag of " << DoubleRej_gsp[i] << " light jet eff of " << AllRej_gsp[i] << endl; 
		cout << " untag jets: "<< untagJets << " total GSP jets: "<< (nDoubleB_gsp[i]+nSingleB_gsp[i]+nNonB_gsp[i]) << " total non-GSP jets: "<< (nDoubleB[i]+nSingleB[i]+nNonB[i]) << endl;
		//cout << "csv value > " << (double)i/40. << " has single-tag eff of " << singleTagEff[i] << " double tag of " << doubleTagEff[i] << " light jet eff of " << untagEff[i] << endl; 
		
	}
	TGraph *doubleTagROC = new TGraph(40, doubleTagEff, untagEff);
	doubleTagROC->SetNameTitle("doubleTagROC","");
	TGraph *singleTagROC = new TGraph(40, singleTagEff, untagEff);
	singleTagROC->SetNameTitle("singleTagROC","");
	TGraph *doubleVSingleROC = new TGraph(40, doubleTagEff, singleTagRej);
	doubleVSingleROC->SetNameTitle("doubleVSingleROC","");
	
	TGraph *doubleTagROC_gsp = new TGraph(40, doubleTagEff_gsp, AllRej_gsp);
	doubleTagROC_gsp->SetNameTitle("doubleTagROC_gsp","");
	TGraph *singleTagROC_gsp = new TGraph(40, singleTagEff_gsp, AllRej_gsp);
	singleTagROC_gsp->SetNameTitle("singleTagROC_gsp","");
	TGraph *noTagROC_gsp = new TGraph(40, untagEff_gsp, AllRej_gsp);
	noTagROC_gsp->SetNameTitle("noTagROC_gsp","");
	
	TGraph *csvGSPeff = new TGraph(40, bareCSV, testDoubleBGSP_csv05);
	csvGSPeff->SetNameTitle("csvGSPeff","");
	TGraph *csvGSPeff_db = new TGraph(40, bareCSV, testnonBGSP_csv05);
	csvGSPeff_db->SetNameTitle("csvGSPeff_db","");
	
	TFile *fout = new TFile("gspDebugHistos_singleFile.root","recreate");
	fout->cd();
	bbInsideJet->Write();
	bbInsideJetFlavor5->Write();
	subjetMatch->Write();
	subjetMatchFlavor5->Write();
	
	gspDR->Write();
	gspxj->Write();
	
	partonSubjetdR->Write();
	dRforbbbar->Scale(1./(double)dRforbbbar->GetEntries());
	dRforbbbar->Write();
	dRforSingleb->Scale(1./(double)dRforSingleb->GetEntries());
	dRforSingleb->Write();
	dRforIncl->Scale(1./(double)dRforIncl->GetEntries());
	dRforIncl->Write();
	dRforNonbbbar->Scale(1./(double)dRforNonbbbar->GetEntries());
	dRforNonbbbar->Write();
	
	zgforbbbar->Scale(1./(double)zgforbbbar->GetEntries());
	zgforbbbar->Write();
	zgforSingleb->Scale(1./(double)zgforSingleb->GetEntries());
	zgforSingleb->Write();
	zgforIncl->Scale(1./(double)zgforIncl->GetEntries());
	zgforIncl->Write();
	zgforNonbbbar->Scale(1./(double)zgforNonbbbar->GetEntries());
	zgforNonbbbar->Write();
	
	//new ghosting histograms...
	subjetHadronDR->Scale(1./(double)subjetHadronDR->GetEntries());
	subjetHadronDR->Write();
	bsubjetHadronDR->Scale(1./(double)bsubjetHadronDR->GetEntries());
	bsubjetHadronDR->Write();
	subjetPartonDR->Scale(1./(double)subjetPartonDR->GetEntries());
	subjetPartonDR->Write();
	bsubjetPartonDR->Scale(1./(double)bsubjetPartonDR->GetEntries());
	bsubjetPartonDR->Write();
	
	subjetHadronFlavor->Scale(1./(double)subjetHadronFlavor->GetEntries());
	subjetHadronFlavor->Write();
	subjetPartonFlavor->Scale(1./(double)subjetPartonFlavor->GetEntries());
	subjetPartonFlavor->Write();
	
	subjetZg2DR->Scale(1./(double)subjetZg2DR->GetEntries());
	subjetZg2DR->Write();
	singlebSubjetZg2DR->Scale(1./(double)singlebSubjetZg2DR->GetEntries());
	singlebSubjetZg2DR->Write();
	doublebSubjetZg2DR->Scale(1./(double)doublebSubjetZg2DR->GetEntries());
	doublebSubjetZg2DR->Write();
	
	subjetZg2DR_gsp->Scale(1./(double)subjetZg2DR_gsp->GetEntries());
	subjetZg2DR_gsp->Write();
	singlebSubjetZg2DR_gsp->Scale(1./(double)singlebSubjetZg2DR_gsp->GetEntries());
	singlebSubjetZg2DR_gsp->Write();
	doublebSubjetZg2DR_gsp->Scale(1./(double)doublebSubjetZg2DR_gsp->GetEntries());
	doublebSubjetZg2DR_gsp->Write();
	
	singleTagROC->Write();
	doubleTagROC->Write();
	doubleVSingleROC->Write();
	
	singleTagROC_gsp->Write();
	doubleTagROC_gsp->Write();
	noTagROC_gsp->Write();
	
	csvGSPeff->Write();
	csvGSPeff_db->Write();
	
	fout->Close();
	
}