#include <string>

using namespace std;

void
Fun4All_TestBeam(int nEvents = 10,
    const char *input_file =
        "/gpfs/mnt/gpfs02/sphenix/data/data01/t1044-2016a/fnal/beam/beam_00002078-0000.prdf",
    const char *output_file = "beam_00002078.root")
{
  gSystem->Load("libfun4all");
  gSystem->Load("libPrototype2.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(Fun4AllServer::VERBOSITY_SOME);

  recoConsts *rc = recoConsts::instance();
  //rc->set_IntFlag("RUNNUMBER",0);

  // ------------------- HCal and EMcal -------------------
  SubsysReco *unpack = new CaloUnpackPRDF();
// unpack->Verbosity(1);
  se->registerSubsystem(unpack);

  CaloCalibration * calib = NULL;

  calib = new CaloCalibration("CEMC");
  calib->GetCalibrationParameters().ReadFromFile("xml",
      string(getenv("CALIBRATIONROOT")) + string("/Prototype2/Calibration/")); // calibration database
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

  // ------------------- Hodoscpes -------------------

  GenericUnpackPRDF *gunpack = NULL;

  const int N_hodo = 8;

  gunpack = new GenericUnpackPRDF("HODO_VERTICAL");
  for (int i = 0; i < N_hodo; ++i)
    gunpack->add_channel(21101, 96 + i, i); // 24 Cerenkov 1
  se->registerSubsystem(gunpack);

  gunpack = new GenericUnpackPRDF("HODO_HORIZONTAL");
  for (int i = 0; i < N_hodo; ++i)
    gunpack->add_channel(21101, 104 + i, i); // 24 Cerenkov 1
  se->registerSubsystem(gunpack);

  calib = new CaloCalibration("HODO_VERTICAL");
  calib->GetCalibrationParameters().set_int_param("calib_const_scale", 1);
  calib->GetCalibrationParameters().set_int_param("use_chan_calibration", 1);
  // Martin find that even channel has negative polarity and odd channel has positive polarity
  for (int i = 0; i < N_hodo; ++i)
    calib->GetCalibrationParameters().set_double_param(
        Form("calib_const_column0_row%d", i), ((i % 2 > 0) ? +1 : -1));
  se->registerSubsystem(calib);

  calib = new CaloCalibration("HODO_HORIZONTAL");
  calib->GetCalibrationParameters().set_int_param("calib_const_scale", 1);
  calib->GetCalibrationParameters().set_int_param("use_chan_calibration", 1);
  // Martin find that even channel has negative polarity and odd channel has positive polarity
  for (int i = 0; i < N_hodo; ++i)
    calib->GetCalibrationParameters().set_double_param(
        Form("calib_const_column0_row%d", i), ((i % 2 > 0) ? +1 : -1));
  se->registerSubsystem(calib);

  // ------------------- Other detectors -------------------

  gunpack = new GenericUnpackPRDF("C1");
// unpack->Verbosity(1);
  gunpack->add_channel(21102, 24, 0); // 24 Cerenkov 1
  se->registerSubsystem(gunpack);

  calib = new CaloCalibration("C1");
  se->registerSubsystem(calib);

  gunpack = new GenericUnpackPRDF("C2");
// unpack->Verbosity(1);
  gunpack->add_channel(21102, 25, 0); //25 Cerenkov 2 Inner
  gunpack->add_channel(21102, 26, 1); //26  Cerenkov 2 Outer
  se->registerSubsystem(gunpack);

  calib = new CaloCalibration("C2");
  se->registerSubsystem(calib);

//  16  HCAL S1 In front of Inner HCAL
//  17  HCAL S2 In "cryostat"
//  18  HCAL S3 On wall behind Outer HCAL
  gunpack = new GenericUnpackPRDF("HCAL_SCINT");
// unpack->Verbosity(1);
  gunpack->add_channel(21102, 16, 1);
  gunpack->add_channel(21102, 17, 2);
  gunpack->add_channel(21102, 18, 3);
  se->registerSubsystem(gunpack);

  calib = new CaloCalibration("HCAL_SCINT");
  se->registerSubsystem(calib);

  gunpack = new GenericUnpackPRDF("PbGL");
// unpack->Verbosity(1);
  gunpack->add_channel(21102, 0, 0); // 0 PbGL  Only inserted in beam for testing
  se->registerSubsystem(gunpack);

  calib = new CaloCalibration("PbGL");
  se->registerSubsystem(calib);

//  28  Bottom trigger veto
//  29  Top trigger veto
//  30  Left trigger veto
//  31  Right trigger veto
  gunpack = new GenericUnpackPRDF("TRIGGER_VETO");
// unpack->Verbosity(1);
  gunpack->add_channel(21102, 28, 0);
  gunpack->add_channel(21102, 29, 1);
  gunpack->add_channel(21102, 30, 2);
  gunpack->add_channel(21102, 31, 3);
  se->registerSubsystem(gunpack);

  calib = new CaloCalibration("TRIGGER_VETO");
  se->registerSubsystem(calib);

  const int N_TileMapper = 16;

  gunpack = new GenericUnpackPRDF("TILE_MAPPER");
  for (int i = 0; i < N_hodo; ++i)
    gunpack->add_channel(21102, 32 + i, i); // 24 Cerenkov 1
  se->registerSubsystem(gunpack);

  calib = new CaloCalibration("TILE_MAPPER");
  se->registerSubsystem(calib);

  // -------------------  Output -------------------
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

  reader->AddTower("RAW_HODO_VERTICAL");
  reader->AddTower("RAW_HODO_HORIZONTAL");
  reader->AddTower("CALIB_HODO_VERTICAL");
  reader->AddTower("CALIB_HODO_HORIZONTAL");

  reader->AddTower("RAW_C1");
  reader->AddTower("CALIB_C1");

  reader->AddTower("RAW_C2");
  reader->AddTower("CALIB_C2");

  reader->AddTower("RAW_HCAL_SCINT");
  reader->AddTower("CALIB_HCAL_SCINT");

  reader->AddTower("RAW_PbGL");
  reader->AddTower("CALIB_PbGL");

  reader->AddTower("RAW_TRIGGER_VETO");
  reader->AddTower("CALIB_TRIGGER_VETO");

  reader->AddTower("RAW_TILE_MAPPER");
  reader->AddTower("CALIB_TILE_MAPPER");

  se->registerSubsystem(reader);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
  in->fileopen(input_file);
  se->registerInputManager(in);

  se->run(nEvents);

  se->End();

}
