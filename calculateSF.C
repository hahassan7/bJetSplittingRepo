void calculateSF(double cut=0.8){

    bool calcVsZg = true; //switch between calculating the scaleFactors vs Zg instead of pT

    //DONT change this - the light SF's are calculated through the negativeCSV tagger now!
    bool doBSF = true; //switch between b or light SFs

    gROOT->ProcessLine(Form(".x setupCFIT.C++(%g)",cut));
    gSystem->Load("CFIT/libCFIT.so");
    gROOT->ProcessLine(Form(".x cfit_gsp.C(%d,%g,%d)",doBSF,cut,calcVsZg));

}
