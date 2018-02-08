#ifndef HEADERS__H_
#define HEADERS__H_

const int ptBins = 6;
const bool doSubjets = true;
//#if doSubjets 
double xbins[ptBins+1] = {10,30,50,70,100,150,200};
//#else
//double xbins[ptBins+1] = {80, 100, 120, 150, 400};
//#endif

int zgBins = 8;
//double xzgBins[8] = {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5}; //zgbins
double xzgBins[9] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4}; //dr-bins (easier to keep same name)

#endif
