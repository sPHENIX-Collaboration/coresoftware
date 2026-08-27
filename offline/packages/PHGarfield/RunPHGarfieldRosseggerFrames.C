R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libPHGarfield.so)

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phgarfield/PHGarfieldRossegger.h>

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

int RunPHGarfieldRosseggerFrames(
    const std::string& garfield_output = "frames_side0_2d_rossegger_v6_4_m1Ez_geom.root",
    const std::string& diagnostic_field_3d_output = "frames_side0_3d_v6_4_m1Ez_geom.root",
    const double frame_reference_phi = 0.0,
    const bool write_diagnostic_field_3d = true,
    const unsigned int tpc_side_for_2d_output = 0,
    const unsigned int mode_phi_max = 24,
    const unsigned int n_radial_modes = 64,
    const unsigned int n_longitudinal_modes = 12,
    const unsigned int radial_job_index = 0,
    const unsigned int n_radial_jobs = 1)
{
  if (fs::exists(garfield_output) &&
      (!write_diagnostic_field_3d || fs::exists(diagnostic_field_3d_output)))
  {
    std::cout << "Requested frame-only Rossegger outputs already exist. Skipping calculation.\n"
              << "  " << garfield_output << "\n";
    if (write_diagnostic_field_3d) { std::cout << "  " << diagnostic_field_3d_output << "\n"; }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto* se = Fun4AllServer::instance();
  auto* rossegger = new PHGarfieldRossegger("PHGarfieldRosseggerFrames");
 
  rossegger->setGeometryCm(21.78, 76.28, 102.325);
  rossegger->setSourceRadiusCm(21.78, 76.28);

  rossegger->setPhiModulation(0.0, 0.0, 0.0, 0.0);


  rossegger->setSourceGrid(18, 144, 24);
  rossegger->setObservationGrid(192, 240, 24);
  rossegger->setModeTruncation(mode_phi_max, n_radial_modes, n_longitudinal_modes);


  rossegger->setAutoAxisymmetric(false);

  rossegger->setTpcSide(tpc_side_for_2d_output);
  rossegger->setRadialJob(radial_job_index, n_radial_jobs);

  rossegger->setUseRealTpcSourceGeometry(false);
  rossegger->setUseFrameChargeModel(true);
  rossegger->setFrameReferencePhi(frame_reference_phi);
  rossegger->setFrameBoundaryPotential(1.0);
  rossegger->setFrameChargeWeighting(PHGarfieldRossegger::FrameChargeWeighting::ProportionalToArea);
  rossegger->setDivideChargeByGain(false);
  rossegger->setNormalizeGainWeightedTotal(false);

  rossegger->setGarfieldOutputFile(garfield_output);
  rossegger->setWriteField3D(write_diagnostic_field_3d);
  rossegger->setField3DOutputFile(diagnostic_field_3d_output);
  rossegger->setWritePHGarfieldField3D(false);
  rossegger->setVerifyOutput(true);

  se->registerSubsystem(rossegger);

  const int status = se->run(1);
  se->End();
  return status;
}