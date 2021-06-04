// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_FUN4ALLANALYSIS_C
#define MACRO_FUN4ALLANALYSIS_C

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <LiteCaloEval.h>


R__LOAD_LIBRARY(libfun4all.so)

void Fun4All_G4_Analysis(
    const int nEvents = 1,
    //   const string &inputTracksFileList = "dst_tracks.list",
    const string &inputClustersFileList = "dst_calo_cluster.list",
    const char * outputFile =  "test1"
)
{
// this convenience library knows all our i/o objects so you don't
// have to figure out what is in each dst type 
  gSystem->Load("libg4dst.so");
  gSystem->Load("libLiteCaloEvalTowSlope.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(); // set it to 1 if you want event printouts

// here you create and register your analysis module like:
// MyModule *mod = new MyModule();
//  se->registerSubsystem(mod);

/*
   Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTTrks");
   in->AddListFile(inputTracksFileList);
   se->registerInputManager(in);
*/
   Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTClusters");
  //in = new Fun4AllDstInputManager("DSTClusters");
   in->AddListFile(inputClustersFileList);
   se->registerInputManager(in);
   
   string outputfile = string(outputFile) + "_g4cemc_eval.root";
   string outputfile2 = string(outputFile) + "_g4hcalin_eval.root";
   string outputfile3 = string(outputFile) + "_g4hcalout_eval.root";
 
//  if (do_hcalin_eval) HCALInner_Eval(string(outputFile) + "_g4hcalin_eval.root")
  
//f  if (do_hcalout_eval) HCALOuter_Eval(string(outputFile) + "_g4hcalout_eval.root");
//string(outputFile) + "_g4cemc_eval.root");

   LiteCaloEval *eval = new LiteCaloEval("CEMCEVALUATOR", "CEMC", outputfile);
//  eval->Verbosity(verbosity);
   eval->user_choice = 'e';
   se->registerSubsystem(eval);

   LiteCaloEval *eval2 = new LiteCaloEval("HINEVALUATOR", "HCALIN", outputfile2);   
   eval2->user_choice = 'i';
//  eval->Verbosity(verbosity);
   se->registerSubsystem(eval2);

   LiteCaloEval *eval3 = new LiteCaloEval("HOUTEVALUATOR", "HCALOUT", outputfile3);
//  eval->Verbosity(verbosity);
   eval3->user_choice = 'o';
   se->registerSubsystem(eval3);

   se->run(nEvents);
   se->End();
   delete se;
   gSystem->Exit(0);

}

#endif //MACRO_FUN4ALLANALYSIS_C
