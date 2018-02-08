
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TLatex.h"

void makeMultiPanelCanvasWithGap(TCanvas*& canv,
								 const Int_t columns,
								 const Int_t rows,
								 const Float_t leftOffset,
								 const Float_t bottomOffset,
								 const Float_t leftMargin,
								 const Float_t bottomMargin,
								 const Float_t edge, const Float_t asyoffset) {
	if (canv==0) {
		Error("makeMultiPanelCanvas","Got null canvas.");
		return;
	}
	canv->Clear();
	
	TPad* pad[columns][rows];
	
	Float_t Xlow[columns];
	Float_t Xup[columns];
	Float_t Ylow[rows];
	Float_t Yup[rows];
	
	Float_t PadWidth =
	(1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
					  (1.0/(1.0-edge))+(Float_t)columns-2.0);
	Float_t PadHeight =
	(1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
						(1.0/(1.0-edge))+(Float_t)rows-2.0);
	
	//PadHeight = 0.5*PadHeight;
	
	Xlow[0] = leftOffset;
	Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
	Xup[columns-1] = 1;
	Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
	
	Yup[0] = 1;
	Ylow[0] = 1.0-PadHeight/(1.0-edge);
	Ylow[rows-1] = bottomOffset;
	Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
	
	for(Int_t i=1;i<columns-1;i++) {
		Xlow[i] = Xup[0] + (i-1)*PadWidth;
		Xup[i] = Xup[0] + (i)*PadWidth;
	}
	Int_t ct = 0;
	for(Int_t i=rows-2;i>0;i--) {
		if(i==rows-2){
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
		}else{
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
			//Yup[i] = 0.2*Yup[i];
		}
		ct++;
	}
	
	TString padName;
	for(Int_t i=0;i<columns;i++) {
		for(Int_t j=0;j<rows;j++) {
			canv->cd();
			padName = Form("p_%d_%d",i,j);
			//pad[i][j] = new TPad(padName.Data(),padName.Data(),
			//Xlow[i],Ylow[j],Xup[i],Yup[j]);
			// this is hacked version to create aysmmetric pads around low 
			if(j==0){
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
			}else{
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
			}
			
			
			if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
			else pad[i][j]->SetLeftMargin(0);
			
			if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
			else pad[i][j]->SetRightMargin(0);
			
			if(j==0) pad[i][j]->SetTopMargin(edge);
			//else pad[i][j]->SetTopMargin(0.01);
			else pad[i][j]->SetTopMargin(0.02);
			
			//if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
			//else pad[i][j]->SetTopMargin(0.0);
			
			//put the gap here
			//if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
			/*else*/ pad[i][j]->SetBottomMargin(0.15);
			
			pad[i][j]->Draw();
			pad[i][j]->cd();
			pad[i][j]->SetNumber(columns*j+i+1);
			
		}
	}
}

//newton's convergent method... x_k+1 = 1/n[(n-1)x_k + A/x_k^(n-1)]
double nthRoot(double p, int n){
	
	double initGuess = 0.9;
	
	double deltaX=999;
	int iter=1;
	while(abs(deltaX)>0.01){
		deltaX = (1./n)*((p/TMath::Power(initGuess,(n-1))) - initGuess);
		iter++;
		initGuess+=deltaX;
	}
	
	return initGuess;
	
}

void toyZg(){
	
	const int Njets = 1000000;
	
	TFile *fin = new TFile("~/zooKeeper/JetRAA_datapoints.root");
	TGraphAsymmErrors *jetRAA = (TGraphAsymmErrors*)fin->Get("RAA_R3_cent0");
	TF1 *raaFit = new TF1("raaFit","[0]+x*[1]+x*x*[2]+x*x*x*[3]",40,300);
	
	TF1 *splitFun = new TF1("splitFun","[0]+[1]/x+[2]/x/x",0.01,0.5);
	
	double xbins[8]={7.4, 3.8, 2.2, 1.6, 1.3, 1.1, 1.05, 1.0};
	TH1D *zgFromData = new TH1D("zgFromData","",8,0.1,0.5);
	for(int i=1; i<=8; i++){ zgFromData->SetBinContent(i, xbins[i-1]); zgFromData->SetBinError(i, 0.01);}
	
	TH1D *energyLost = new TH1D("energyLost","",100,0,50); energyLost->Sumw2();
	TH2D *energyLossMap = new TH2D("elossMap","",100,0,500,100,0,40);
	
	TH2D *zgCorr = new TH2D("zgCorr","",50,0,0.5,50,0,0.5);
	
	//140-160, 160-180, 180-200, 200-250, 250-300, 300-500
	double binning[6] = {160,180,200,250,300,500};
	TH1D *zgppDist[6], *zgPbPbDist[6];
	for(int i=0; i<6; i++){
		zgppDist[i] = new TH1D(Form("zgDist_pp_%d",i),"",10,0.,0.5); zgppDist[i]->Sumw2();
		zgPbPbDist[i] = new TH1D(Form("zgDist_PbPb_%d",i),"",10,0.,0.5); zgPbPbDist[i]->Sumw2();
	}
	
	TF1 *jetSpectrum = new TF1("jetSpectrum","TMath::Power(x,-5)",140,500);
	TF1 *gausSmear = new TF1("gausSmear","gaus",0.5,1.5);
	gausSmear->SetParameter(0,1); gausSmear->SetParameter(1,1); gausSmear->SetParameter(2,0.05);
	
	//raaFit->SetParameter(0,0.4); raaFit->SetParameter(1,0.1); raaFit->SetParameter(2,0); raaFit->SetParameter(3,0);
	jetRAA->Fit(raaFit,"","",40,250);
	TCanvas *c1 = new TCanvas("c1","",600,600);
	c1->cd();
	jetRAA->Draw();
	jetRAA->GetXaxis()->SetTitle("Jet p_{T}");
	jetRAA->GetYaxis()->SetTitle("RAA");
	
	//splitFun->SetParameter(0,1); splitFun->SetParameter(1,0); splitFun->SetParameter(2,0);
	zgFromData->Fit(splitFun,"qN0","",0.01,0.5);
	//zgFromData->Draw();
	zgFromData->GetYaxis()->SetTitle("1/N_{jet} dN/dz_{g}");
	zgFromData->GetXaxis()->SetTitle("z_{g}");
	
	TRandom3 *rand = new TRandom3(0);
	
	for(int ijet=0; ijet<Njets; ijet++){
		
		//double zgsplit = rand->Rndm()/2.;
		double zgsplit = splitFun->GetRandom();
		//if(zgsplit<0.1) continue;
		double fullJetPt = jetSpectrum->GetRandom();
		
		int bin=0;
		while(fullJetPt > binning[bin]) bin++;
		
		double subleadSub = fullJetPt*zgsplit;
		double leadSub = fullJetPt*(1-zgsplit);
		
		//if(subleadSub<40) continue;
		
		//cout << "full jet: "<< fullJetPt << " split: "<< zgsplit << endl;
		//cout << "leading: "<< leadSub << " sublead: "<< subleadSub << endl;
		//cout << "bin: "<< bin << endl;
		
		double ppZg = subleadSub/(leadSub+subleadSub);
		//cout << "pp zg: "<< ppZg << endl;
		zgppDist[bin]->Fill(ppZg);
		
		//do energy loss...
		double elossToEval = leadSub > 250 ? raaFit->Eval(250):raaFit->Eval(leadSub);
		double leadSub_post=leadSub*nthRoot(elossToEval, 5)*gausSmear->GetRandom();
		elossToEval = subleadSub > 250 ? raaFit->Eval(250):raaFit->Eval(subleadSub);
		double subleadSub_post=subleadSub*nthRoot(elossToEval, 5)*gausSmear->GetRandom();
		
		//energyLost->Fill(leadSub*(1-nthRoot(raaFit->Eval(leadSub), 5)));
		energyLost->Fill(subleadSub-subleadSub_post);
		energyLossMap->Fill(subleadSub, subleadSub-subleadSub_post);
		
		double PbPbZg = subleadSub_post/(leadSub_post+subleadSub_post);
		
		if(PbPbZg<0.1) continue;
		//cout << "PbPb zg: "<< PbPbZg << endl;
		zgPbPbDist[bin]->Fill(PbPbZg);
				
		zgCorr->Fill(ppZg, PbPbZg);
	}
	
	TCanvas *c2 = new TCanvas("c2","",1200,800);
	c2->Divide(3,2);
	for(int i=0; i<6; i++){
		c2->cd(i+1);
		zgPbPbDist[i]->Scale(1./zgPbPbDist[i]->Integral(3,10));
		zgppDist[i]->SetMarkerColor(2); zgppDist[i]->SetLineColor(2);
		zgppDist[i]->Scale(1./zgppDist[i]->Integral(3,10));
		zgPbPbDist[i]->Draw();
		zgppDist[i]->Draw("same");
	}
	
	TLatex *l1[6];
	TH1D *zgRatio[6];
	TCanvas *cc = new TCanvas("cc","",1200,800);
	makeMultiPanelCanvasWithGap(cc, 3, 2, 0., 0., 0.18, 0.18, 0.1, 0);
	//cc->Divide(3,2);
	for(int i=0; i<6; i++){
		cc->cd(i+1);
		zgRatio[i] = (TH1D*)zgPbPbDist[i]->Clone(Form("zgRatio_%d",i));
		zgRatio[i]->Divide(zgppDist[i]);
		zgRatio[i]->SetYTitle("PbPb / pp");
		zgRatio[i]->SetXTitle("z_{g}");
		zgRatio[i]->SetMaximum(1.75);
		zgRatio[i]->SetMinimum(0.25);
		zgRatio[i]->Draw();
		if(i==0) l1[i] = new TLatex(0.2,0.4,Form("140 < p_{T,jet} < %.0f",binning[i]));
		else l1[i] = new TLatex(0.2,0.4,Form("%.0f < p_{T,jet} < %.0f",binning[i-1],binning[i]));
		l1[i]->SetTextSize(0.08);
		l1[i]->Draw("Same");
	}
	
	TCanvas *c4 = new TCanvas("c4","",600,600);
	c4->cd();
	energyLost->Draw();
	
	TCanvas *c5 = new TCanvas("c5","",600,600);
	c5->cd();
	energyLossMap->Draw("colz");
	
	TCanvas *c6 = new TCanvas("c6","",600,600);
	c6->cd();
	zgCorr->Draw("colz");
}