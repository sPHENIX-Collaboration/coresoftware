// PHGarfield Rossegger using the r-phi density maps produced by the notebook.

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libPHGarfield.so)

#include <fun4all/Fun4AllServer.h>
#include <phgarfield/PHGarfieldRossegger.h>

#include <string>

int RunPHGarfieldRosseggerDensityMap(
    const std::string& density_file = "include/TPC_IBF_primary_maps.root",
    const bool use_ibf = true,
    const bool normalize_to_reference = true,
    const double reference_density_nC_m3 = 20.0,
    const std::string& garfield_output = "notebook_ibf_rossegger.root",
    const std::string& field3d_output = "notebook_ibf_fields_3d.root",
    const std::string& phgarfield_side0 = "notebook_ibf_field_side0_South.root",
    const std::string& phgarfield_side1 = "notebook_ibf_field_side1_North.root")
{
  auto* se = Fun4AllServer::instance();
  auto* r = new PHGarfieldRossegger("PHGarfieldRossegger");

  const std::string stem = use_ibf ? "h_ibf_final_rphi_side" : "h_primary_final_rphi_side";
  r->setGeometryCm(21.78,76.28,102.325);
  r->setSourceRadiusCm(21.78,75.43);
  r->setDensity(reference_density_nC_m3,1.0,2.0);
  r->setUseDensityMap(true);
  r->setDensityMapFile(density_file,stem+"0",stem+"1");
  r->setNormalizeDensityMap(normalize_to_reference);

  r->setSourceGrid(48,216,48);
  r->setObservationGrid(96,108,48);
  r->setModeTruncation(24,12,12);
  r->setAutoAxisymmetric(false);

  r->setGarfieldOutputFile(garfield_output);
  r->setWriteField3D(true);
  r->setField3DOutputFile(field3d_output);
  r->setWritePHGarfieldField3D(true);
  r->setPHGarfieldField3DOutputFiles(phgarfield_side0,phgarfield_side1);
  r->setVerifyOutput(true);

  se->registerSubsystem(r);
  const int status=se->run(1);
  se->End();
  return status;
}
