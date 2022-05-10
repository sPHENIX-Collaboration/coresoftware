
//double fitf(Double_t * f, Double_t *p);
//TGraph * grff = 0;
#include <GlobalVariables.C>

#include "LiteCaloEval.h"
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

void do_eta_fit2(const char * reffile, const char * infile, const char * modfile)
{
  gSystem->Load("libLiteCaloEvalTowSlope.so");
  LiteCaloEval reflce, modlce;
  reflce.CaloType(LiteCaloEval::CEMC);
  modlce.CaloType(LiteCaloEval::CEMC);
  reflce.Get_Histos(reffile);
  modlce.Get_Histos(infile,modfile);
  modlce.fitmin = 0.12;
  modlce.fitmax = 0.7;
  modlce.FitRelativeShifts(&reflce,1);
  


}


