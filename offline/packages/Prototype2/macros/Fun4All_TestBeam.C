#include <string>

using namespace std;

void Fun4All_TestBeam(
          const char *input_file= "/gpfs/mnt/gpfs02/sphenix/data/data01/t1044-2016a/fnal/beam/beam_00002070-0000.prdf",
          const char *output_file = "data/TB_DST.root",
	  int nEvents = 1000)
{
 gSystem->Load("libfun4all");
 gSystem->Load("libPrototype2.so");

 Fun4AllServer *se = Fun4AllServer::instance();
 se->Verbosity(Fun4AllServer::VERBOSITY_SOME);

 recoConsts *rc = recoConsts::instance();
 //rc->set_IntFlag("RUNNUMBER",0);

 SubsysReco *unpack = new CaloUnpackPRDF();
// unpack->Verbosity(1);
 se->registerSubsystem( unpack );

 //main DST output
 Fun4AllDstOutputManager *out_Manager  = new Fun4AllDstOutputManager("DSTOUT",output_file);
// se->registerOutputManager( out_Manager );

 //alternatively, fast check on DST using DST Reader:
 Prototype2DSTReader *reader = new Prototype2DSTReader(string(output_file) + string("_DSTReader.root"));
 reader->AddTower("RAW_LG_HCALIN");
 reader->AddTower("RAW_HG_HCALIN");
 reader->AddTower("RAW_LG_HCALOUT");
 reader->AddTower("RAW_HG_HCALOUT");
 reader->AddTower("RAW_CEMC");
 se->registerSubsystem( reader );


 Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
 in->fileopen(input_file);
 se->registerInputManager(in);

 se->run(nEvents);

 se->End(); 

}
