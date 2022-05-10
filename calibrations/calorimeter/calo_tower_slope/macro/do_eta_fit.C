
//double fitf(Double_t * f, Double_t *p);
//TGraph * grff = 0;
#include <GlobalVariables.C>

#include "LiteCaloEval.h"
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

void do_eta_fit(const char * reffile, const char * modfile)
{
  gSystem->Load("libLiteCaloEvalTowSlope.so");
  LiteCaloEval reflce, modlce;
  reflce.CaloType(LiteCaloEval::HCALOUT);
  modlce.CaloType(LiteCaloEval::HCALOUT);
  reflce.Get_Histos(reffile);
  modlce.Get_Histos(modfile);
  modlce.FitRelativeShifts(&reflce,21);

}


