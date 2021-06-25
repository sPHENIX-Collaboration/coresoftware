// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_FUN4ALLG4SLOPECAL_C
#define MACRO_FUN4ALLG4SLOPECAL_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <litecaloeval/LiteCaloEval.h>

R__LOAD_LIBRARY(libfun4all.so)

void Fun4All_G4_SlopeCal(
    const int nEvents = 1,
    const string &inputClustersFileList = "dst_calo_cluster.list",
    const string &outputFile = "test1")
{
  // this convenience library knows all our i/o objects so you don't
  // have to figure out what is in each dst type
  gSystem->Load("libg4dst.so");
  gSystem->Load("libLiteCaloEvalTowSlope.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity();  // set it to 1 if you want event printouts

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTClusters");
  in->AddListFile(inputClustersFileList);
  se->registerInputManager(in);

  string outputfile = outputFile + "_g4cemc_eval.root";
  string outputfile2 = outputFile + "_g4hcalin_eval.root";
  string outputfile3 = outputFile + "_g4hcalout_eval.root";

  LiteCaloEval *eval = new LiteCaloEval("CEMCEVALUATOR", "CEMC", outputfile);
  //  eval->Verbosity(verbosity);
  eval->CaloType(LiteCaloEval::CEMC);
  se->registerSubsystem(eval);

  LiteCaloEval *eval2 = new LiteCaloEval("HINEVALUATOR", "HCALIN", outputfile2);
  eval2->CaloType(LiteCaloEval::HCALIN);
  //  eval->Verbosity(verbosity);
  se->registerSubsystem(eval2);

  LiteCaloEval *eval3 = new LiteCaloEval("HOUTEVALUATOR", "HCALOUT", outputfile3);
  //  eval->Verbosity(verbosity);
  eval3->CaloType(LiteCaloEval::HCALOUT);
  se->registerSubsystem(eval3);

  se->run(nEvents);
  se->End();
  delete se;
  gSystem->Exit(0);
}

#endif  //MACRO_FUN4ALLG4SLOPECAL_C
