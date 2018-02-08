# bJetSplittingRepo
HIN-18-002 analysis macros

Typical Workflow:
1) scanGSPEvents.C (creates working ntuples)
2) calculateSF.C (wrapper for setupCFIT.C and cfit_gsp.C to create the data efficiency templates)
3) Rename output of calculateSF.C to something useful
4) Load the output as part of the efficiency correction tables in L70-91 of drawZg.C
5) drawZg.C (creates final plots)
