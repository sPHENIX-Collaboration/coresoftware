#include <string>

using namespace std;

void
Fun4All_TestBeam(int nEvents = 100,
    const char *input_file =
        "/gpfs/mnt/gpfs02/sphenix/data/data01/t1044-2016a/fnal/beam/beam_00002078-0000.prdf",
    const char *output_file = "beam_00002078.root")
{
  gSystem->Load("libfun4all");
  gSystem->Load("/gpfs/mnt/gpfs02/phenix/scratch/abhisek/coresoftware/offline/packages/Prototype2/build/lib/libPrototype2.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(Fun4AllServer::VERBOSITY_SOME);

  recoConsts *rc = recoConsts::instance();
  //rc->set_IntFlag("RUNNUMBER",0);

  SubsysReco *unpack = new CaloUnpackPRDF();
// unpack->Verbosity(1);
  se->registerSubsystem(unpack);

  CaloCalibration * calib = NULL;

  calib = new CaloCalibration("CEMC");
  //calib->GetCalibrationParameters().ReadFromFile("xml",
  //    string(getenv("CALIBRATIONROOT")) + string("/Prototype2/Calibration/")); // calibration database
  se->registerSubsystem(calib);

  calib = new CaloCalibration("HCALIN");
  calib->set_calib_tower_node_prefix("CALIB_LG");
  calib->set_raw_tower_node_prefix("RAW_LG");
  //calib->Verbosity(true);
  se->registerSubsystem(calib);

  calib = new CaloCalibration("HCALIN");
  calib->set_calib_tower_node_prefix("CALIB_HG");
  calib->set_raw_tower_node_prefix("RAW_HG");
  se->registerSubsystem(calib);

  calib = new CaloCalibration("HCALOUT");
  calib->set_calib_tower_node_prefix("CALIB_LG");
  calib->set_raw_tower_node_prefix("RAW_LG");
  se->registerSubsystem(calib);

  calib = new CaloCalibration("HCALOUT");
  calib->set_calib_tower_node_prefix("CALIB_HG");
  calib->set_raw_tower_node_prefix("RAW_HG");
  se->registerSubsystem(calib);

  //main DST output
  Fun4AllDstOutputManager *out_Manager = new Fun4AllDstOutputManager("DSTOUT",
      output_file);
  se->registerOutputManager(out_Manager);

  //alternatively, fast check on DST using DST Reader:
  Prototype2DSTReader *reader = new Prototype2DSTReader(
      string(output_file) + string("_DSTReader.root"));
  reader->AddTower("RAW_LG_HCALIN");
  reader->AddTower("RAW_HG_HCALIN");
  reader->AddTower("RAW_LG_HCALOUT");
  reader->AddTower("RAW_HG_HCALOUT");
  reader->AddTower("RAW_CEMC");
  reader->AddTower("CALIB_LG_HCALIN");
  reader->AddTower("CALIB_HG_HCALIN");
  reader->AddTower("CALIB_LG_HCALOUT");
  reader->AddTower("CALIB_HG_HCALOUT");
  reader->AddTower("CALIB_CEMC");
  se->registerSubsystem(reader);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
  in->fileopen(input_file);
  se->registerInputManager(in);

  se->run(nEvents);

  se->End();

}
