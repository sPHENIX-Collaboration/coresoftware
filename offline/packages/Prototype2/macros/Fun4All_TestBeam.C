void Fun4All_TestBeam(
          const char *input_file= "/sphenix/data/data01/t1044-2016a/cosmics/cosmics_00001554-0000.prdf",
          const char *output_file = "TB_DST.root",
	  int nEvents = 4000)
{
 gSystem->Load("libfun4all");
 gSystem->Load("/direct/phenix+u/abhisek/HCAL_DST/build/lib/libHCal.so");

 Fun4AllServer *se = Fun4AllServer::instance();
 se->Verbosity(0);

 recoConsts *rc = recoConsts::instance();
 //rc->set_IntFlag("RUNNUMBER",0);

 SubsysReco *unpack = new CaloUnpackPRDF();
// unpack->Verbosity(1);
 se->registerSubsystem( unpack );

 Fun4AllDstOutputManager *out_Manager  = new Fun4AllDstOutputManager("DSTOUT",output_file);
 se->registerOutputManager( out_Manager );

 Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
 in->fileopen(input_file);
 se->registerInputManager(in);

 se->run(nEvents);

 se->End(); 

}
