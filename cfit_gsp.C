#include "CFIT/cfit.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

void getResults(CFIT::cfit *cf,float *par,float *err);

void cfit_gsp()
{
   gROOT->SetBatch();
   
   gSystem->Load("CFIT/libCFIT.so");

   float par[100];
   float err[100];
   float par_tag[100];
   float err_tag[100];
   
   const int ptBins =3;
   
   for(int ipt=0; ipt<ptBins; ipt++){
   
   
   	CFIT::cfit *cf = new CFIT::cfit("Jet Probability");
   	cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
   	cf->SetMorphing(OPTMORPH_GEOMETRIC);
   	cf->ProducePlots(1);
   	
   	cf->SetInputFile("cfitInput.root");

   //cf->AddSys("SYS1","_sys1_down","_sys1_up");
   //cf->AddSys("SYS2","_sys2_down","_sys2_up");

   	cf->SetMatrixOption("WRITE");
   	
   	cf->SetData(Form("h_data_%d",ipt));
   	cf->SetDataTag(Form("h_data_tag_%d",ipt));
   	cf->SetDataUntag(Form("h_data_untag_%d",ipt));
   	
   	cf->AddTemplate("usdg",Form("h_usdg_mc_%d",ipt),2);
   	cf->AddTemplate("c",Form("h_c_mc_%d",ipt),3);
   	cf->AddTemplate("b",Form("h_b_mc_%d",ipt),4);

   	cf->AddTemplateTag("usdg",Form("h_usdg_mc_tag_%d",ipt),2);
   	cf->AddTemplateTag("c",Form("h_c_mc_tag_%d",ipt),3);
   	cf->AddTemplateTag("b",Form("h_b_mc_tag_%d",ipt),4);

   	cf->AddTemplateUntag("usdg",Form("h_usdg_mc_untag_%d",ipt),2);
   	cf->AddTemplateUntag("c",Form("h_c_mc_untag_%d",ipt),3);
   	cf->AddTemplateUntag("b",Form("h_b_mc_untag_%d",ipt),4);
   	
   	std::cout << "ptBin : "<< ipt << endl;
   	
   	cf->Run();      	
   	getResults(cf,par,err);
   	float chi2 = cf->GetChisq();
   	float ndata = cf->GetNData();
   	float nmc1 = cf->GetNTemplate("b");
   	float nmc = cf->GetNTemplate("usdg")+cf->GetNTemplate("c")+cf->GetNTemplate("b");
   	
   	cout << "untag chi2: "<< chi2 << endl;
   	cout << "untag ndof: "<< cf->GetNDOF() << endl;
   	cout << "p (bjet): " << par[0] << " +/- " << err[0] << endl;
   	cout << "p (cjet): " << par[1] << " +/- " << err[1] << endl;
   	cout << "p (ljet): " << par[2] << " +/- " << err[2] << endl;

   	cf->Run("tag");   
   	getResults(cf,par_tag,err_tag);
   	float chi2_tag = cf->GetChisq();
   	float ndata_tag = cf->GetNData();
   	float nmc1_tag = cf->GetNTemplate("b");
   	float nmc_tag = cf->GetNTemplate("usdg")+cf->GetNTemplate("c")+cf->GetNTemplate("b");
   	
   	cout << "tag chi2: "<< chi2_tag << endl;
   	cout << "tag ndof: "<< cf->GetNDOF() << endl;
   	cout << "p (bjet): " << par_tag[0] << " +/- " << err_tag[0] << endl;
   	cout << "p (cjet): " << par_tag[1] << " +/- " << err_tag[1] << endl;
   	cout << "p (ljet): " << par_tag[2] << " +/- " << err_tag[2] << endl;

   	float fr = nmc1/nmc;
   	float fr_tag = nmc1_tag/nmc_tag;
   	
   	float effMC = nmc1_tag/nmc1;
   	float effErrMC = sqrt(pow(effMC,2)*(1./nmc1_tag + 1./nmc1));
   	float effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr;
   	
   	//std::cout << "effMC = " << effMC << std::endl;
   	//std::cout << "effDATA = " << effDATA << std::endl;
   	//std::cout << "sf = " << effDATA/effMC << std::endl;

   // perform statistical variation
   	cf->SetMatrixOption("READ");
   	cf->SetStatVariation(667);
   	cf->ProducePlots(0);   
   	cf->Run();
   	getResults(cf,par,err);
   	float nmc1_stat = cf->GetNTemplate("b");
   	float nmc_stat = cf->GetNTemplate("usdg")+cf->GetNTemplate("c")+cf->GetNTemplate("b");
   	
   	cf->SetStatVariation(667);
   	cf->Run("tag");
   	getResults(cf,par_tag,err_tag);
   	float nmc1_tag_stat = cf->GetNTemplate("b");
   	float nmc_tag_stat = cf->GetNTemplate("usdg")+cf->GetNTemplate("c")+cf->GetNTemplate("b");
   	
   // do the calculation of SF here again ....
   	float fr_stat = nmc1_stat/nmc_stat;
   	float fr_tag_stat = nmc1_tag_stat/nmc_tag_stat;
   	float effMC_stat = nmc1_tag_stat/nmc1_stat;
   	float effErrMC_stat = sqrt(pow(effMC_stat,2)*(1./nmc1_tag_stat + 1./nmc1_stat));
   	float effDATA_stat = nmc_tag_stat/nmc_stat*par_tag[0]/par[0]*fr_tag_stat/fr_stat;
   	std::cout << "effMC = " << effMC << " +/- " << abs(effMC-effMC_stat) << std::endl;
   	std::cout << "effDATA = " << effDATA << " +/- " << abs(effDATA-effDATA_stat) << std::endl;
   	
   	float sfErr = sqrt(pow(effDATA/effMC,2)*(pow(abs(effMC-effMC_stat)/effMC,2)+pow(abs(effDATA-effDATA_stat)/effDATA,2)));
   	std::cout << "sf stat = " << effDATA/effMC << " +/- "<< sfErr << std::endl;
   	
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
