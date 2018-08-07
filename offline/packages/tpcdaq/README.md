# Welcome

This is an example analys module which iterate through reconstructed jet, match it with cloest truth jet and tracks within a dR cone.
The output is a ROOT file with histograms of reconstructed jet sepectrum and a TTree with per-jet kinematics and matched truth jet kinematics. 
The matched track pT and Delta-R is also saved in an array in the TTree. 

## Compile this module 

To compile this module, please follow instruction of [building a package](https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes#Building_a_package), 

After compiling, it should produce this file: `$MYINSTALL/lib/libmyjetanalysis.so`

## Run this module

To run this module, please insert this block of code in the Fun4All macro for your jet analysis:
```c++

Fun4All_Macro(...)
{
  //.... ( Run this after jet and tracking recontstruction  )
  
  // this loads the library for my jet analysis you have just compiled. 
  gSystem->Load("libmyjetanalysis");

  // this loads the library for my jet analysis you have just compiled. 
  MyJetAnalysis *myJetAnalysis = new MyJetAnalysis("AntiKt_Tower_r04","AntiKt_Truth_r04","myjetanalysis.root");
  se->registerSubsystem(myJetAnalysis);
  
  //.... (other IO operations)
  
  // this run number of events
  se->run(nEvents);
 
  // histogram and TTree are saved at this step
  se->End();
}
```
Furthermore, here is a full set of example macro to run this analysis module 
by analyzing existing single particle (32 GeV pi-) simulation production 
(therefore always makes single particle jet with 1 matched track).
* Repository: `blackcathj/macros.git`, branch: `CD-1-reference-production-readback-my-jet-analysis`
* i.e. https://github.com/blackcathj/macros/tree/CD-1-reference-production-readback-my-jet-analysis/macros/g4simulations
* Main macro: [Fun4All_G4_sPHENIX.C](https://github.com/blackcathj/macros/blob/CD-1-reference-production-readback-my-jet-analysis/macros/g4simulations/Fun4All_G4_sPHENIX.C)

## Check the output

The main output file is `myjetanalysis.root`, which contains
```
root [1] .ls
TFile**		myjetanalysis.root	
 TFile*		myjetanalysis.root	
  KEY: TH1F	hInclusive_E;1	AntiKt_Tower_r04 inclusive jet E
  KEY: TH1F	hInclusive_eta;1	AntiKt_Tower_r04 inclusive jet #eta
  KEY: TH1F	hInclusive_phi;1	AntiKt_Tower_r04 inclusive jet #phi
  KEY: TTree	T;1	MyJetAnalysis Tree
```
The T tree contains one entry per reconstructed jet. Here is its first entry:
```
root [2] T->Show(0)
======> EVENT:0
 event           = 0          <- event number
 id              = 49         <- jet ID in this event
 nComponent      = 51         <- number of component in reco jet, i.e. # of towers in jet
 event           = 0
 phi             = -1.43342
 e               = 19.4207
 pt              = 18.2251
 truthID         = 0
 truthNComponent = 1
 truthEta        = 0.348408
 truthPhi        = -1.42448
 truthE          = 32.0003
 truthPt         = 30.1514
 nMatchedTrack   = 1
 trackdR         = 0.00869318
 trackpT         = 30.4937
```

## Make this your module

Please copy this folder to a new folder under your local [analysis repository](https://github.com/sPHENIX-Collaboration/analysis) and change package and class name. 

# Read more

Many next-step topics are listed in the [software](https://wiki.bnl.gov/sPHENIX/index.php/Software) page. And specifically, to use the simulation for your study, a few thing you might want to try:

* [More on write your analysis module for more dedicated analysis](https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes), 
* Examples to use evaluators of track, caloiemter, and jets in your analysis module ([CaloEvaluator](https://www.phenix.bnl.gov/WWW/sPHENIX/doxygen/html/dd/d59/classCaloEvaluator.html), [JetEvaluator](https://www.phenix.bnl.gov/WWW/sPHENIX/doxygen/html/d1/df4/classJetEvaluator.html), [SvtxEvaluator](https://www.phenix.bnl.gov/WWW/sPHENIX/doxygen/html/d6/d11/classSvtxEvaluator.html)). Welcome to copy and resuse blocks of the code from them.
