R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libPHGarfield.so)

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phgarfield/PHGarfieldRossegger.h>

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

int RunPHGarfieldRossegger(
    const std::string& garfield_output = "simple_2d_rossegger.root",
    const bool use_real_tpc_source_geometry = true,
    const bool use_frame_charge_model = false,
    const double frame_reference_phi = 0.0,
    const std::string& pad_placement_file = "include/TPC_pad_placement.txt",
    const std::string& gain_map_file = "include/layer_gain_79513_Mariia_side01.root",
    const bool write_diagnostic_field_3d = true,
    const std::string& diagnostic_field_3d_output = "sphenix_rossegger_fields_3d.root",
    const bool write_phgarfield_field_3d = true,
    const std::string& phgarfield_field_3d_side0 = "sphenix_3d_ibf_field_side0_South_v3.root",
    const std::string& phgarfield_field_3d_side1 = "sphenix_3d_ibf_field_side1_North_v3.root",
    const unsigned int tpc_side_for_2d_output = 0,
    const unsigned int radial_job_index = 0,
    const unsigned int n_radial_jobs = 1)
{
  if (fs::exists(garfield_output) &&
      (!write_diagnostic_field_3d || fs::exists(diagnostic_field_3d_output)) &&
      (!write_phgarfield_field_3d ||
       (fs::exists(phgarfield_field_3d_side0) && fs::exists(phgarfield_field_3d_side1))))
  {
    std::cout << "Requested Rossegger map outputs already exist. Skipping calculation.\n"
              << "  " << garfield_output << "\n";
    if (write_diagnostic_field_3d) { std::cout << "  " << diagnostic_field_3d_output << "\n"; }
    if (write_phgarfield_field_3d)
    {
      std::cout << "  " << phgarfield_field_3d_side0 << "\n"
                << "  " << phgarfield_field_3d_side1 << "\n";
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto* se = Fun4AllServer::instance();
  auto* rossegger = new PHGarfieldRossegger("PHGarfieldRossegger");

  rossegger->setGeometryCm(21.78, 76.28, 102.325);
  rossegger->setSourceRadiusCm(use_frame_charge_model ? 21.78 : 22.8, use_frame_charge_model ? 76.28 : 75.43);


  rossegger->setDensity(20.0,1.0,1.8);

  rossegger->setSourceGrid(48, 216, 48);
  rossegger->setObservationGrid(96, 108, 48);///takes like 5 min but it is too precise
  rossegger->setModeTruncation(24, 12, 12);
  rossegger->setAutoAxisymmetric(false);

  rossegger->setTpcSide(tpc_side_for_2d_output);
  rossegger->setRadialJob(radial_job_index, n_radial_jobs);

  rossegger->setUseRealTpcSourceGeometry(use_real_tpc_source_geometry && !use_frame_charge_model);
  rossegger->setUseFrameChargeModel(use_frame_charge_model);
  rossegger->setFrameReferencePhi(frame_reference_phi);
  rossegger->setPadPlacementFile(pad_placement_file);
  rossegger->setGainMapFile(gain_map_file);
  rossegger->setGainHistograms("hGainMap_side0_South", "hGainMap_side1_North");
  rossegger->setDivideChargeByGain(true);
  rossegger->setNormalizeGainWeightedTotal(true);

  rossegger->setGarfieldOutputFile(garfield_output);
  rossegger->setWriteField3D(write_diagnostic_field_3d);
  rossegger->setField3DOutputFile(diagnostic_field_3d_output);
  rossegger->setWritePHGarfieldField3D(write_phgarfield_field_3d);
  rossegger->setPHGarfieldField3DOutputFiles(phgarfield_field_3d_side0, phgarfield_field_3d_side1);
  rossegger->setVerifyOutput(true);

  se->registerSubsystem(rossegger);

  const int status = se->run(1);
  se->End();
  return status;
}
