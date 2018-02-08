
{
TFile *f0 = new TFile("out_PbPbMC_QCDjet_noCuts.root");
TFile *f1 = new TFile("out_PbPbMC_QCDjet_subjetDR0p2.root");
TFile *f2 = new TFile("gspTreeOut_PbPbMC_explicitJTA.root");
TFile *f3 = new TFile("out_ppMC_QCDjet_noCuts_fixROC.root");
TFile *f4 = new TFile("gspTreeOut_ppMC_Pythia6_OfficialQCDJetMC_explicitJTA_recalibJP_withNegCSV.root");
TFile *f5 = new TFile("out_ppMC_QCDjet_rejectPSVtx_fixROC.root");
TFile *f6 = new TFile("out_ppMC_QCDjet_rejectPSVtx_subjetDR0p2_fixROC.root");

TGraph *PbPb_cent0 = (TGraph*)f2->Get("purEffCurve_drCut_cent0");
TGraph *PbPb_cent1 = (TGraph*)f2->Get("purEffCurve_drCut_cent1");
TGraph *PbPb_cent2 = (TGraph*)f2->Get("purEffCurve_drCut_cent2");
TGraph *ppNo = (TGraph*)f4->Get("purEffCurve_noCuts_cent0");
TGraph *ppDef = (TGraph*)f4->Get("purEffCurve_drCut_cent0");
TGraph *ppFullVtx = (TGraph*)f4->Get("purEffCurve_vtxTypeCut_cent0");
TGraph *ppFullVtxDR = (TGraph*)f4->Get("purEffCurve_vtxType_drCut_cent0");

PbPb_cent0->SetPoint(39,1,0);
PbPb_cent1->SetPoint(39,1,0);
PbPb_cent2->SetPoint(39,1,0);
ppNo->SetPoint(39,1,0);
ppDef->SetPoint(39,1,0);
ppFullVtx->SetPoint(39,1,0);
ppFullVtxDR->SetPoint(39,1,0);

TH1D *hh = new TH1D("hh","",1,0,1);
hh->SetMaximum(1);
hh->SetMinimum(0);
hh->SetXTitle("CSV GSP-tagging purity");
hh->SetYTitle("CSV GSP-tagging efficiency");

hh->Draw();
PbPb_cent0->SetLineColor(kMagenta);
PbPb_cent0->SetMarkerColor(kMagenta);
//PbPb_cent0->Draw("pl*,same");
PbPb_cent1->SetLineColor(kMagenta+1);
PbPb_cent1->SetMarkerColor(kMagenta+1);
//PbPb_cent1->Draw("pl*,same");
PbPb_cent2->SetLineColor(kMagenta+3);
PbPb_cent2->SetMarkerColor(kMagenta+3);
//PbPb_cent2->Draw("pl*,same");
ppNo->SetLineColor(kGreen);
ppNo->SetMarkerColor(kGreen);
//ppNo->Draw("pl*,same");
ppFullVtxDR->SetLineColor(kBlue+3);
ppFullVtxDR->SetMarkerColor(kBlue+3);
//ppFullVtxDR->Draw("pl*,same");
ppDef->SetLineColor(kGreen+3);
ppDef->SetMarkerColor(kGreen+3);
ppDef->Draw("pl*,same");
ppFullVtx->SetLineColor(kBlue);
ppFullVtx->SetMarkerColor(kBlue);
//ppFullVtx->Draw("pl*,same");

TLegend *l1 = new TLegend(0.1,0.1,0.4,0.3);
/*l1->AddEntry(PbPb_cent0,"PbPb MC, 0-10% subjet #DeltaR > 0.2","lp");
l1->AddEntry(PbPb_cent1,"PbPb MC, 10-30% subjet #DeltaR > 0.2","lp");
l1->AddEntry(PbPb_cent2,"PbPb MC, 30-100% subjet #DeltaR > 0.2","lp");*/
//l1->AddEntry(ppNo,"pp MC","lp");
l1->AddEntry(ppDef,"pp MC, subjet #DeltaR > 0.1","lp");
//l1->AddEntry(ppFullVtx,"pp MC, Full-Vtx Only","lp");
//l1->AddEntry(ppFullVtxDR,"pp MC, Full-Vtx + #DeltaR > 0.1","lp");
l1->Draw("same");

}