#include "CFIT/cfit.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "headers.h"

double getHistErr(TH1D *histo){
    double err=0;
    histo->IntegralAndError(1,histo->GetNbinsX(), err);
    return err;
}

void getResults(CFIT::cfit *cf,float *par,float *err);

//decide between calculating the b or light scale factor
void cfit_gsp(bool calcSFb=1, double cut=0.7, bool calcVsZg=false)
{
    gROOT->SetBatch();
   
   gSystem->Load("CFIT/libCFIT.so");

   float par[100];
   float err[100];
   float par_tag[100];
   float err_tag[100];
   string sfTypeString = calcSFb ? "btagSF" : "btagLight";
  
   int bins=ptBins;
   double *binRange = xbins;
   if(calcVsZg){
        bins=zgBins;
        binRange=xzgBins;
   }

   TH1D *scaleFactors = new TH1D(Form("scaleFactors_%s",sfTypeString.c_str()),"",bins,binRange);  scaleFactors->Sumw2(); 
   TH1D *sfLight = new TH1D("sfLight","",bins,binRange); sfLight->Sumw2();
   TH1D *mcEff = new TH1D("mcEff","",bins,binRange); mcEff->Sumw2();
   TH1D *dataEff = new TH1D("dataEff","",bins,binRange); dataEff->Sumw2();
   TH1D *mcPur = new TH1D("mcPur","",bins,binRange); mcPur->Sumw2();
   TH1D *dataPur = new TH1D("dataPur","",bins,binRange); dataPur->Sumw2();
   
   TH1D *mcLightFrac = new TH1D("mcLightFrac","",bins,binRange); mcLightFrac->Sumw2();
   TH1D *dataLightFrac = new TH1D("dataLightFrac","",bins,binRange); dataLightFrac->Sumw2();
   TH1D *mcBFrac = new TH1D("mcBFrac","",bins,binRange); mcBFrac->Sumw2();
   TH1D *dataBFrac = new TH1D("dataBFrac","",bins,binRange); dataBFrac->Sumw2();

   TCanvas *cMaster = new TCanvas("cMaster","cMaster");
    
   string filename;
   if(doSubjets) filename = "cfitInput_subjets.root";
   else filename = "cfitInput_fullJets.root";
    
   string appendix="";
   if(calcVsZg) appendix="zg_";
   int nBins = calcVsZg ? zgBins : ptBins; 
   for(int ipt=0; ipt<nBins; ipt++){
  
       if(calcVsZg && ipt<2) continue; //skipping over bins without content (zg<0.1)

   	CFIT::cfit *cf = new CFIT::cfit("Jet Probability");
   	cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
   	cf->SetMorphing(OPTMORPH_GEOMETRIC);
   	cf->ProducePlots(1, Form("_%sbin%d",appendix.c_str(),ipt));
        cf->SetVerbose(1);

        cf->SetInputFile(filename.c_str());

   //cf->AddSys("SYS1","_sys1_down","_sys1_up");
   //cf->AddSys("SYS2","_sys2_down","_sys2_up");

   	cf->SetMatrixOption("WRITE");
   	
   	cf->SetData(Form("h_data_%s%d",appendix.c_str(),ipt));
   	cf->SetDataTag(Form("h_data_tag_%s%d",appendix.c_str(),ipt));
   	cf->SetDataUntag(Form("h_data_untag_%s%d",appendix.c_str(),ipt));
   	
   	cf->AddTemplate("b",Form("h_b_mc_%s%d",appendix.c_str(),ipt),2);
   	cf->AddTemplate("c",Form("h_c_mc_%s%d",appendix.c_str(),ipt),3);
   	cf->AddTemplate("usdg",Form("h_usdg_mc_%s%d",appendix.c_str(),ipt),4);

   	cf->AddTemplateTag("b",Form("h_b_mc_tag_%s%d",appendix.c_str(),ipt),2);
   	cf->AddTemplateTag("c",Form("h_c_mc_tag_%s%d",appendix.c_str(),ipt),3);
   	cf->AddTemplateTag("usdg",Form("h_usdg_mc_tag_%s%d",appendix.c_str(),ipt),4);

   	cf->AddTemplateUntag("b",Form("h_b_mc_untag_%s%d",appendix.c_str(),ipt),2);
   	cf->AddTemplateUntag("c",Form("h_c_mc_untag_%s%d",appendix.c_str(),ipt),3);
   	cf->AddTemplateUntag("usdg",Form("h_usdg_mc_untag_%s%d",appendix.c_str(),ipt),4);
   	
   	std::cout << "bin : "<< ipt << endl;
   	
   	cf->Run("");      	
   	getResults(cf,par,err);
   	float chi2 = cf->GetChisq();
   	float ndata = cf->GetNData();
   	float nmc1;
        if(calcSFb) nmc1 = cf->GetNTemplate("b");
        else nmc1 = cf->GetNTemplate("usdg");
        float nmcb = cf->GetNTemplate("b");
        std::cout << "nmc1 = " << nmc1 << std::endl;
   	float nmc = cf->GetNTemplate("usdg")+cf->GetNTemplate("b")+cf->GetNTemplate("c");
   	float nmc_err = TMath::Sqrt(pow(cf->GetErrTemplate("usdg"),2)+pow(cf->GetErrTemplate("b"),2)+pow(cf->GetErrTemplate("c"),2));

   	cout << "untag chi2: "<< chi2 << endl;
   	cout << "untag ndof: "<< cf->GetNDOF() << endl;
   	cout << "p (bjet): " << par[0] << " +/- " << err[0] << endl;
   	cout << "p (cjet): " << par[1] << " +/- " << err[1] << endl;
   	cout << "p (ljet): " << par[2] << " +/- " << err[2] << endl;

   	cf->Run("tag");
   	getResults(cf,par_tag,err_tag);
   	float chi2_tag = cf->GetChisq();
   	float ndata_tag = cf->GetNData();
   	float nmc1_tag, nmc1_tag_err;
        if(calcSFb){
            nmc1_tag = cf->GetNTemplate("b");
            nmc1_tag_err = cf->GetErrTemplate("b");
        }
        else{
            nmc1_tag = cf->GetNTemplate("usdg");
            nmc1_tag_err = cf->GetErrTemplate("usdg");
        }
        float nmcb_tag = cf->GetNTemplate("b");
        std::cout << "nmc1_tag = " << nmc1_tag << std::endl;
   	float nmc_tag = cf->GetNTemplate("usdg")+cf->GetNTemplate("b")+cf->GetNTemplate("c");
        float nmc_tag_err = TMath::Sqrt(pow(cf->GetErrTemplate("usdg"),2)+pow(cf->GetErrTemplate("b"),2)+pow(cf->GetErrTemplate("c"),2));
   	
   	cout << "tag chi2: "<< chi2_tag << endl;
   	cout << "tag ndof: "<< cf->GetNDOF() << endl;
   	cout << "p (bjet): " << par_tag[0] << " +/- " << err_tag[0] << endl;
   	cout << "p (cjet): " << par_tag[1] << " +/- " << err_tag[1] << endl;
   	cout << "p (ljet): " << par_tag[2] << " +/- " << err_tag[2] << endl;

   	float fr = nmc1/nmc;
   	float fr_tag = nmc1_tag/nmc_tag;
   
        mcLightFrac->SetBinContent(ipt+1,(1-fr));
        mcLightFrac->SetBinError(ipt+1, (1-fr)*(1./nmc1 + 1./nmc));
        mcBFrac->SetBinContent(ipt+1,fr);
        mcBFrac->SetBinError(ipt+1, fr*(1./nmc1 + 1./nmc));

   	float effMC = nmc1_tag/nmc1;
   	float effErrMC = sqrt(pow(effMC,2)*(1./nmc1_tag + 1./nmc1));
   	cout << "nmc_tag: "<< nmc_tag << " nmc: "<< nmc << " par_tag[0]: "<< par_tag[0] << " par[0] " << par[0] << " fr_tag: "<< fr_tag << " fr: " << fr << endl;
        float effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr;

        float purMC = nmcb_tag/nmcb;

        //do the mistag stuff too...
        TFile *fin = new TFile(filename.c_str());
        TH1D *h_usdg_mc_negTag = (TH1D*)fin->Get(Form("h_usdg_mc_negTag_%s%d",appendix.c_str(),ipt))->Clone("h_usdg_mc_negTag");
        TH1D *h_c_mc_negTag = (TH1D*)fin->Get(Form("h_c_mc_negTag_%s%d",appendix.c_str(),ipt))->Clone("h_c_mc_negTag");
        TH1D *h_b_mc_negTag = (TH1D*)fin->Get(Form("h_b_mc_negTag_%s%d",appendix.c_str(),ipt))->Clone("h_b_mc_negTag");
        
        TH1D *h_usdg_mc_negTagCut = (TH1D*)fin->Get(Form("h_usdg_mc_negTagCut_%s%d",appendix.c_str(),ipt))->Clone("h_usdg_mc_negTagCut");
        TH1D *h_c_mc_negTagCut = (TH1D*)fin->Get(Form("h_c_mc_negTagCut_%s%d",appendix.c_str(),ipt))->Clone("h_c_mc_negTagCut");
        TH1D *h_b_mc_negTagCut = (TH1D*)fin->Get(Form("h_b_mc_negTagCut_%s%d",appendix.c_str(),ipt))->Clone("h_b_mc_negTagCut");

        TH1D *h_data_negTag = (TH1D*)fin->Get(Form("h_data_negTag_%s%d",appendix.c_str(),ipt))->Clone("h_data_negTag");
        TH1D *h_data_negTagCut = (TH1D*)fin->Get(Form("h_data_negTagCut_%s%d",appendix.c_str(),ipt))->Clone("h_data_negTagCut");

        float mistagMC = (nmc_tag - nmc1_tag)/nmc_tag;
        float mistagMCErr = TMath::Sqrt(pow(nmc_tag_err,2)+pow(nmc1_tag_err,2));
        mistagMCErr = mistagMC*TMath::Sqrt(pow(mistagMCErr/(nmc_tag - nmc1_tag),2)+pow(nmc_tag_err/nmc_tag,2));
        float negEffMC_top = (h_usdg_mc_negTagCut->Integral()+h_b_mc_negTagCut->Integral()+h_c_mc_negTagCut->Integral());
        float negEffMC_bot = (h_usdg_mc_negTag->Integral()+h_b_mc_negTag->Integral()+h_c_mc_negTag->Integral());
        float negEffMCErr_top = TMath::Sqrt(pow(getHistErr(h_usdg_mc_negTagCut),2)+pow(getHistErr(h_b_mc_negTagCut),2)+pow(getHistErr(h_c_mc_negTagCut),2)); 
        float negEffMCErr_bot = TMath::Sqrt(pow(getHistErr(h_usdg_mc_negTag),2)+pow(getHistErr(h_b_mc_negTag),2)+pow(getHistErr(h_c_mc_negTag),2));
        float negEffMCErr = (negEffMC_top/negEffMC_bot)*TMath::Sqrt(pow(negEffMCErr_top/negEffMC_top,2)+pow(negEffMCErr_bot/negEffMC_bot,2));
        float negEffData = h_data_negTagCut->Integral()/h_data_negTag->Integral();
        float negEffDataErr = negEffData*TMath::Sqrt(pow(getHistErr(h_data_negTagCut)/h_data_negTagCut->Integral(),2)+pow(getHistErr(h_data_negTag)/h_data_negTag->Integral(),2));
        float mistagData = negEffData*mistagMC/(negEffMC_top/negEffMC_bot);
        float mistagDataErr = mistagData*TMath::Sqrt(pow(negEffDataErr/negEffData,2)+pow(mistagMCErr/mistagMC,2)+pow(negEffMCErr_top/negEffMC_top,2)+pow(negEffMCErr_bot/negEffMC_bot,2));

   	sfLight->SetBinContent(ipt+1, mistagData/mistagMC);
        sfLight->SetBinError(ipt+1, mistagData/mistagMC*TMath::Sqrt(pow(mistagDataErr/mistagData,2)+pow(mistagMCErr/mistagMC,2)));
        
        //std::cout << "effMC = " << effMC << std::endl;
   	//std::cout << "effDATA = " << effDATA << std::endl;
   	//std::cout << "sf = " << effDATA/effMC << std::endl;

   // perform statistical variation
   	cf->SetMatrixOption("READ");
   	cf->SetStatVariation(667);
   	cf->ProducePlots(0);   
   	cf->Run("");
   	getResults(cf,par,err);
   	float nmc1_stat;
        if(calcSFb) nmc1_stat = cf->GetNTemplate("b");
        else nmc1_stat = cf->GetNTemplate("usdg");
        std::cout << "nmc1_stat = " << nmc1_stat << std::endl;
   	float nmc_stat = cf->GetNTemplate("usdg")+cf->GetNTemplate("b")+cf->GetNTemplate("c");
   	
   	cf->SetStatVariation(667);
   	cf->Run("tag");
   	getResults(cf,par_tag,err_tag);
   	float nmc1_tag_stat;
        if(calcSFb) nmc1_tag_stat = cf->GetNTemplate("b");
        else nmc1_tag_stat = cf->GetNTemplate("usdg");
        std::cout << "nmc1_tag_stat = " << nmc1_tag_stat << std::endl;
        float nmc_tag_stat = cf->GetNTemplate("usdg")+cf->GetNTemplate("b")+cf->GetNTemplate("c");
   	
   // do the calculation of SF here again ....
   	float fr_stat = nmc1_stat/nmc_stat;
   	float fr_tag_stat = nmc1_tag_stat/nmc_tag_stat;
   	float effMC_stat = nmc1_tag_stat/nmc1_stat;
   	float effErrMC_stat = sqrt(pow(effMC_stat,2)*(1./nmc1_tag_stat + 1./nmc1_stat));
   	float effDATA_stat = nmc_tag_stat/nmc_stat*par_tag[0]/par[0]*fr_tag_stat/fr_stat;
   	std::cout << "effMC = " << effMC << " +/- " << abs(effErrMC - effErrMC_stat) << std::endl;
   	std::cout << "effDATA = " << effDATA << " +/- " << abs(effDATA-effDATA_stat) << std::endl;
   	
   	float sfErr = sqrt(pow(effDATA/effMC,2)*(pow(abs(effErrMC_stat-effErrMC)/effMC,2)+pow(abs(effDATA-effDATA_stat)/effDATA,2)));
   	std::cout << "sf stat = " << effDATA/effMC << " +/- "<< sfErr << std::endl;
   
        mcEff->SetBinContent(ipt+1,effMC);
        mcEff->SetBinError(ipt+1,0.001+abs(effMC-effMC_stat));

        dataEff->SetBinContent(ipt+1,effDATA);
        dataEff->SetBinError(ipt+1,abs(effDATA-effDATA_stat));

        mcPur->SetBinContent(ipt+1, 1-mistagMC);
        mcPur->SetBinError(ipt+1, mistagMCErr);

        dataPur->SetBinContent(ipt+1, 1-mistagData);
        dataPur->SetBinError(ipt+1, mistagDataErr);

        dataBFrac->SetBinContent(ipt+1, (nmc_tag/nmc)*(1-mistagData)/(effDATA));
        dataBFrac->SetBinError(ipt+1, dataBFrac->GetBinContent(ipt+1)*TMath::Sqrt(pow(nmc_tag_err/nmc_tag,2)+pow(nmc_err/nmc,2)+pow(mistagDataErr/mistagData,2)+pow(abs(effDATA-effDATA_stat)/effDATA,2)));
        dataLightFrac->SetBinContent(ipt+1, 1-dataBFrac->GetBinContent(ipt+1));
        dataLightFrac->SetBinError(ipt+1, dataBFrac->GetBinError(ipt+1));

   	scaleFactors->SetBinContent(ipt+1,effDATA/effMC);
        scaleFactors->SetBinError(ipt+1,sfErr);

        std::cout << endl << endl;
   	
   // perform systematic variation
   	cf->SetMatrixOption("READ");
 //  cf->SetSysVariation("_sys1_up");
   //cf->Run();
   //getResults(cf,par,err);
 //  cf->SetSysVariation("_sys1_up");
   //cf->Run("tag");
   //getResults(cf,par,err);
   	
   // do the calculation of SF here again ....
   	
   	delete cf;
   }
   cMaster->cd();
   scaleFactors->Draw(""); 
   scaleFactors->GetYaxis()->SetRangeUser(0,1.1);
   cMaster->Print("pics/scaleFactors_vsPt.eps"); 
   string jetTypeString = doSubjets ? "subjets" : "fullJets";
   TFile *fout = new TFile(Form("sf_%s_%s_muFiltered_csv%g_sv1p5.root",jetTypeString.c_str(),appendix.c_str(),cut),"recreate");
   fout->cd();
   mcEff->Write();
   dataEff->Write();
   mcPur->Write();
   dataPur->Write();
   mcLightFrac->Write();
   dataLightFrac->Write();
   mcBFrac->Write();
   dataBFrac->Write();
   scaleFactors->Write();
   sfLight->Write();
   fout->Close();
   
   gApplication->Terminate();
}

void getResults(CFIT::cfit *cf,float *par,float *err)
{   
   float chi2 = cf->GetChisq();
   int nPar = cf->GetNPar();
   for(int i=0;i<nPar;i++)
     {
	par[i] = cf->GetPar(i);
	err[i] = cf->GetParErr(i);	
     }
}
