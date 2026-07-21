#include "PHGarfield.h"
#include <phool/phool.h>
#include <cdbobjects/CDBTTree.h>

#include <phfield/PHField3DCartesian.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TPolyLine3D.h>
#include <TRotation.h>
#include <TVector3.h>
#include <TFile.h>
#include <TH2.h>
#include <TAxis.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <Garfield/ComponentUser.hh>
#include <Garfield/MediumMagboltz.hh>

#include <cmath>
#include <limits>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <memory>  // for basic_ostream, operat...
#include <vector>
#include <regex>
#include <string>
#include <utility>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace fs = std::filesystem;

PHGarfield::PHGarfield(const std::string& name,
                       const std::string& electricFieldMap,
                       double spaceChargeScale_side0,
                       double spaceChargeScale_side1)
  : SubsysReco(name),
    m_defaultGasfile("/sphenix/user/hemmick/gasfiles_20260624"),
    m_electricFieldMap(electricFieldMap),
    m_spaceChargeScale_side0(spaceChargeScale_side0),
    m_spaceChargeScale_side1(spaceChargeScale_side1)
{
}

PHGarfield::~PHGarfield()
{
  //  Housekeeping.
  delete m_field;
  delete m_cdbTPCMAPttree;
  delete m_component;
  delete m_gas;
  delete m_erCorrection;
  delete m_ezCorrection;
}

int PHGarfield::InitRun(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHGarfield::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;
  }
  CDBInterface* m_cdb = CDBInterface::instance();

  //  Here we use the CDBInterface to set up the magnetic field map:
  std::string url = m_cdb->getUrl("FIELDMAP_TRACKING");
  m_field = new PHField3DCartesian(url, 1.0);

  //  Here we use the CDBInterface to set up the channel making of the TPC:
  std::string text = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");
  m_cdbTPCMAPttree = new CDBTTree(text);
  m_cdbTPCMAPttree->LoadCalibrations();

  // Load the optional axisymmetric space-charge correction map.
  // Failure is non-fatal: Garfield then uses only the nominal 400 V/cm field.
  if (!m_electricFieldMap.empty())
  {
    if (!LoadElectricFieldCorrections(m_electricFieldMap))
    {
      std::cerr << PHWHERE << " Failed to load electric-field correction map: "
                << m_electricFieldMap << std::endl;
    }
  }

  //  Make the Garfield Component and register the methods that will interface to our fields...
  m_component = new Garfield::ComponentUser();
  m_component->SetMagneticField([this](double x, double y, double z, double& bx, double& by, double& bz)
                                { GetMagneticFieldTesla(x, y, z, bx, by, bz); });
  m_component->SetElectricField([this](double x, double y, double z, double& ex, double& ey, double& ez)
                                { GetElectricFieldVcm(x, y, z, ex, ey, ez); });

  // Here we fetch the gas from the CDB
  std::string gasfile = m_cdb->getUrl("PHGARFIELD_GAS");
  if (gasfile.empty() || !fs::exists(gasfile))
    {
      std::cerr << PHWHERE << " Missing CDB gasfile: " << gasfile << std::endl;
      std::cerr << PHWHERE << " Using default gasfile: " << m_defaultGasfile << std::endl;
      gasfile = m_defaultGasfile;
    }
  InitializeGas(gasfile);

  //  Diagnostic during code development...
  FillRadii();
  if (Verbosity() > 1)
  {
    PrintMaps();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGarfield::FillRadii()
{
  //  Unload the pad map to get the radii in a handy location:
  for (unsigned int side = 0; side < 2; side++)
  {
    for (unsigned int sector = 0; sector < 12; sector++)
    {
      for (unsigned int fee = 0; fee < 26; fee++)
      {
        for (unsigned int channel = 0; channel < 256; channel++)
        {
          unsigned int key = (256 * (fee)) + channel;
          int layer = m_cdbTPCMAPttree->GetIntValue(key, "layer");
          double r = m_cdbTPCMAPttree->GetDoubleValue(key, "R") / CLHEP::cm;
          if (layer > 6)
          {
            radii[layer - 7] = r;
          }
        }
      }
    }
  }
}

void PHGarfield::PrintGarfield(double x, double y, double z) const
{
  double ex;
  double ey;
  double ez;
  double bx;
  double by;
  double bz;
  double vx;
  double vy;
  double vz;
  GetElectricFieldVcm(x, y, z, ex, ey, ez);
  GetMagneticFieldTesla(x, y, z, bx, by, bz);
  m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  std::cout << " x:" << x
            << " y:" << y
            << " z:" << z
            << " ex:" << ex
            << " ey:" << ey
            << " ez:" << ez
            << " bx:" << bx
            << " by:" << by
            << " bz:" << bz
            << " vx:" << vx
            << " vy:" << vy
            << " vz:" << vz
            << std::endl;
}

void PHGarfield::PrintGasSummary() const
{
  if (!m_GasFilesLoaded)
    {
      std::cerr << PHWHERE << "No Gas File(s) have been successfully loaded." << std::endl;
      return;
    }

  std::vector<double> nE;
  std::vector<double> nB;
  std::vector<double> nA;
  m_gas->GetFieldGrid(nE, nB, nA);
  
  std::cout << "Gas File Grid Dimensions: " << std::endl;
  std::cout << nE.size() << " E-fields ranging from " << nE.front() << " to " << nE.back() << std::endl; 
  std::cout << nB.size() << " B-fields ranging from " << nB.front() << " to " << nB.back() << std::endl; 
  std::cout << nA.size() << " Angles   ranging from " << nA.front() << " to " << nA.back() << std::endl; 
}

void PHGarfield::PrintMaps() const
{
  //  Print out a few test points of the Garfield information
  PrintGarfield(0.0, 0.0, 0.1);
  PrintGarfield(0.0, 0.0, 100.0);
  PrintGarfield(0.0, 40.0, 100.1);
  PrintGarfield(0.0, 78.0, 010.1);

  //  Print out the pad coordinate map:
  int MAX = 10;
  int prints = 0;
  for (unsigned int side = 0; side < 2; side++)
  {
    for (unsigned int sector = 0; sector < 12; sector++)
    {
      for (unsigned int fee = 0; fee < 26; fee++)
      {
        for (unsigned int channel = 0; channel < 256; channel++)
        {
          unsigned int key = (256 * (fee)) + channel;
          int layer = m_cdbTPCMAPttree->GetIntValue(key, "layer");
          double phi = ((side == 1 ? 1 : -1) * (m_cdbTPCMAPttree->GetDoubleValue(key, "phi") - std::numbers::pi / 2.)) + ((sector % 12) * std::numbers::pi / 6);
          double r = m_cdbTPCMAPttree->GetDoubleValue(key, "R") / CLHEP::cm;

          phi = bounder(phi, PHI_MIN);

          if (layer > 6)
          {
            if (prints < MAX)
            {
              prints++;
              std::cout << " side: " << side;
              std::cout << " sector: " << sector;
              std::cout << " fee: " << fee;
              std::cout << " channel: " << channel;
              std::cout << " layer: " << layer;
              std::cout << " phi: " << phi;
              std::cout << " r: " << r;
              std::cout << std::endl;
            }
          }
        }
      }
    }
  }
}


void PHGarfield::MoveMagnet(double x_cm, double y_cm, double z_cm)
{
  m_magpos.SetXYZ(x_cm, y_cm, z_cm);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: magnetic-field map translation = ("
              << x_cm << ", " << y_cm << ", " << z_cm << ") cm" << std::endl;
  }
}

void PHGarfield::RotateMagnet(double theta_x, double theta_y, double theta_z)
{
  m_magrot.RotateX(theta_x);
  m_magrot.RotateY(theta_y);
  m_magrot.RotateZ(theta_z);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: magnetic-field map rotation increment = ("
              << theta_x << ", " << theta_y << ", " << theta_z << ") rad" << std::endl;
  }
}

void PHGarfield::MoveTpc(double x_cm, double y_cm, double z_cm)
{
  m_tpcpos.SetXYZ(x_cm, y_cm, z_cm);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: TPC translation = ("
              << x_cm << ", " << y_cm << ", " << z_cm << ") cm" << std::endl;
  }
}

void PHGarfield::RotateTpc(double theta_x, double theta_y, double theta_z)
{
  m_tpcrot.RotateX(theta_x);
  m_tpcrot.RotateY(theta_y);
  m_tpcrot.RotateZ(theta_z);
  if (Verbosity() > 0)
  {
    std::cout << "PHGarfield: TPC rotation increment = ("
              << theta_x << ", " << theta_y << ", " << theta_z << ") rad" << std::endl;
  }
}

TVector3 PHGarfield::TpcPointToGlobalPoint(double x_cm, double y_cm, double z_cm) const
{
  // Transform a point from the local TPC/Garfield frame to the global detector frame.
  return m_tpcrot * TVector3(x_cm, y_cm, z_cm) + m_tpcpos;
}

TVector3 PHGarfield::GlobalPointToTpcPoint(double x_cm, double y_cm, double z_cm) const
{
  // Transform a point from the global detector frame to the local TPC/Garfield frame.
  return m_tpcrot.Inverse() * (TVector3(x_cm, y_cm, z_cm) - m_tpcpos);
}

TVector3 PHGarfield::TpcPointToMagnetFieldMapPoint(double x_cm, double y_cm, double z_cm) const
{
  // Garfield coordinates are kept in the local TPC frame.  To sample the
  // sPHENIX magnetic-field map, first place that point into the global frame,
  // then express the global point in the magnetic-field-map frame.
  const TVector3 p_global = TpcPointToGlobalPoint(x_cm, y_cm, z_cm);
  return m_magrot.Inverse() * (p_global - m_magpos);
}

TVector3 PHGarfield::MagnetFieldMapVectorToTpcVector(double bx, double by, double bz) const
{
  // GetFieldValue returns the B vector in the magnetic-field-map axes.  Rotate
  // it to global axes and then into the local TPC axes before handing it to
  // Garfield.  No translation is applied to vectors.
  const TVector3 b_map(bx, by, bz);
  const TVector3 b_global = m_magrot * b_map;
  return m_tpcrot.Inverse() * b_global;
}

void PHGarfield::GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double& bx_t, double& by_t, double& bz_t) const
{
  // NOTE: Garfield uses cm, V/cm, and Tesla.
  //       PHField3DCartesian/CLHEP coordinates use CLHEP length units and
  //       magnetic fields in CLHEP field units.
  //
  // Important: x_cm,y_cm,z_cm are local TPC coordinates.  We keep Garfield and
  // the gas tables in this local frame.  Only the point used to query the
  // magnetic-field map is transformed.

  const TVector3 p_map_cm = TpcPointToMagnetFieldMapPoint(x_cm, y_cm, z_cm);

  double point[4] =
      {
          p_map_cm.X() * CLHEP::cm,
          p_map_cm.Y() * CLHEP::cm,
          p_map_cm.Z() * CLHEP::cm,
          0.0};

  double bfield_map[3] = {0.0, 0.0, 0.0};

  // Get the magnetic field via the PHField3DCartesian object constructed using
  // the CDB url reference.
  m_field->GetFieldValue(point, bfield_map);

  const TVector3 b_tpc = MagnetFieldMapVectorToTpcVector(
      bfield_map[0], bfield_map[1], bfield_map[2]);

  bx_t = b_tpc.X() / CLHEP::tesla;
  by_t = b_tpc.Y() / CLHEP::tesla;
  bz_t = b_tpc.Z() / CLHEP::tesla;
}

void PHGarfield::GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm) const
{
  // NOTE: Garfield uses cm, V/cm, and Tesla.
  // The notebook maps use cm on their axes and V/m in their bins.
  // The map is produced for one TPC half using s = |z|, measured from
  // the central membrane toward the readout plane.

  const double r_cm = std::hypot(x_cm, y_cm);
  const double abs_z_cm = std::abs(z_cm);

  // Nominal drift field.  Electrons drift away from the central membrane,
  // therefore E_z points toward the central membrane on both sides.
  ex_vcm = 0.0;
  ey_vcm = 0.0;
  ez_vcm = z_cm > 0.0 ? -m_CMVoltageDefault : m_CMVoltageDefault;

  const double spaceChargeScale = (z_cm < 0.0) ? m_spaceChargeScale_side0
                                               : m_spaceChargeScale_side1;

  if (!m_erCorrection || !m_ezCorrection || spaceChargeScale == 0.0)
  {
    return;
  }

  const double delta_er_vcm = spaceChargeScale *
                              InterpolateCorrectionVcm(m_erCorrection, r_cm, abs_z_cm);
  const double delta_ez_local_vcm = spaceChargeScale *
                                    InterpolateCorrectionVcm(m_ezCorrection, r_cm, abs_z_cm);

  // Convert the cylindrical radial correction to Cartesian components.
  if (r_cm > 0.0)
  {
    ex_vcm += delta_er_vcm * x_cm / r_cm;
    ey_vcm += delta_er_vcm * y_cm / r_cm;
  }

  // hEzDefault is expressed along the local coordinate s = |z|.
  // Convert it to the global Cartesian z direction.
  ez_vcm += z_cm >= 0.0 ? delta_ez_local_vcm : -delta_ez_local_vcm;
}

bool PHGarfield::LoadElectricFieldCorrections(const std::string& filename)
{
  std::unique_ptr<TFile> input(TFile::Open(filename.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cerr << PHWHERE << " Could not open electric-field map: "
              << filename << std::endl;
    return false;
  }

  auto* er = dynamic_cast<TH2*>(input->Get("QA/hErDefault"));
  auto* ez = dynamic_cast<TH2*>(input->Get("QA/hEzDefault"));

  // Also allow maps written at the ROOT-file top level.
  if (!er) er = dynamic_cast<TH2*>(input->Get("hErDefault"));
  if (!ez) ez = dynamic_cast<TH2*>(input->Get("hEzDefault"));

  if (!er || !ez)
  {
    std::cerr << PHWHERE
              << " Missing QA/hErDefault or QA/hEzDefault in "
              << filename << std::endl;
    return false;
  }

  delete m_erCorrection;
  delete m_ezCorrection;
  m_erCorrection = dynamic_cast<TH2*>(er->Clone("PHGarfield_ErCorrection"));
  m_ezCorrection = dynamic_cast<TH2*>(ez->Clone("PHGarfield_EzCorrection"));

  if (!m_erCorrection || !m_ezCorrection)
  {
    delete m_erCorrection;
    delete m_ezCorrection;
    m_erCorrection = nullptr;
    m_ezCorrection = nullptr;
    return false;
  }

  m_erCorrection->SetDirectory(nullptr);
  m_ezCorrection->SetDirectory(nullptr);

  std::cout << "Loaded axisymmetric electric-field corrections from "
            << filename << std::endl;
  std::cout << "  scale k_eff side0/south/z<0 = "
          << m_spaceChargeScale_side0 << std::endl;
  std::cout << "  scale k_eff side1/north/z>0 = "
          << m_spaceChargeScale_side1 << std::endl;
  std::cout << "  r range [cm] = ["
            << m_erCorrection->GetXaxis()->GetXmin() << ", "
            << m_erCorrection->GetXaxis()->GetXmax() << "]" << std::endl;
  std::cout << "  |z| range [cm] = ["
            << m_erCorrection->GetYaxis()->GetXmin() << ", "
            << m_erCorrection->GetYaxis()->GetXmax() << "]" << std::endl;

  return true;
}

double PHGarfield::InterpolateCorrectionVcm(
    const TH2* hist,
    const double r_cm,
    const double abs_z_cm) const
{
  if (!hist)
  {
    return 0.0;
  }

  const TAxis* xaxis = hist->GetXaxis();
  const TAxis* yaxis = hist->GetYaxis();

  if (!xaxis || !yaxis)
  {
    return 0.0;
  }

  // Interpolate safely between bin centers, not at histogram edges.
  const double r_min = xaxis->GetBinCenter(1);
  const double r_max = xaxis->GetBinCenter(xaxis->GetNbins());

  const double z_min = yaxis->GetBinCenter(1);
  const double z_max = yaxis->GetBinCenter(yaxis->GetNbins());

  constexpr double epsilon = 1.0e-6;

  const double r_eval =
      std::clamp(r_cm, r_min + epsilon, r_max - epsilon);

  const double z_eval =
      std::clamp(std::abs(abs_z_cm),
                 z_min + epsilon,
                 z_max - epsilon);

  // Input histogram is in V/m. Garfield expects V/cm.
  return 0.01 * hist->Interpolate(r_eval, z_eval);
}

void PHGarfield::InitializeGas(const std::string &name)
{
  //  Create and fill the gas object so that we can trace particles through the gas...
  m_gas = new Garfield::MediumMagboltz();

  if (!std::filesystem::exists(name))
  {
    std::cerr << "Missing gas file or gas directory: " << name << std::endl;
    return;
  }

  if (fs::is_regular_file(name))
    {
      std::cout << "Loading Garfield gas from file: " << name << std::endl;
      if (!m_gas->LoadGasFile(name))
	{
	  std::cerr << "Failed to load " << name << std::endl;
	  return;
	}
      m_GasFilesLoaded = true;
    }
  else if (fs::is_directory(name))
    {
      std::cout << "Loading Garfield gas from directory: " << name << std::endl;
      std::regex filePattern(R"(^MERGED_E([0-9]{3})\.gas$)");
      std::smatch matchResults;

      // Iterate through all items in the directory
      // NOTE:  Map assures that files are properly ordered when merged...
      std::map<unsigned int, std::string> FilesToMerge;
      for (const auto& entry : fs::directory_iterator(name))
	{
	  // Only process regular files
	  if (entry.is_regular_file())
	    {
	      std::string filepath = entry.path().string();
	      std::string filename = entry.path().filename().string();
	      
	      // Check if the filename matches our target pattern
	      if (std::regex_match(filename, matchResults, filePattern))
		{
		  //std::cout << "matchResults: " << matchResults[1].str() << std::endl;
		  FilesToMerge[std::stoul( matchResults[1].str() )]=filepath;
		}
	    }
	}
      bool firstE = true;
      for (const auto& [key, filepath] : FilesToMerge)
	{
	  if (firstE)
	    {
	      m_gas->LoadGasFile(filepath);
	      firstE = false;
	      m_GasFilesLoaded = true;
	    }
	  else
	    {
	      m_gas->MergeGasFile(filepath, true);
	      m_GasFilesLoaded = true;
	    }
	}
    }

  PrintGasSummary();
}

int PHGarfield::process_event(PHCompositeNode* topNode)
{
  // Initial implementation doesn't do anything event-by-event.
  // Nonetheless, a future user might want do do something here...
  (void) topNode;
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHGarfield::bounder(double phi, double phi_min)
{
  double phi_max = phi_min + 2.0 * std::numbers::pi;
  while (phi < phi_min)
  {
    phi = phi + 2.0 * std::numbers::pi;
  }
  while (phi >= phi_max)
  {
    phi = phi - 2.0 * std::numbers::pi;
  }

  return phi;
}

TPolyLine3D* PHGarfield::ReverseDrift(double x, double y, double z, double step_ns)
{
  // Public/default drift function: input and output are local TPC coordinates.
  // Do not transform these coordinates before asking Garfield for E/gas values;
  // otherwise the pre-generated gas tables and nominal TPC E field can be sampled
  // outside their intended frame.
  std::vector<double> xlist;
  std::vector<double> ylist;
  std::vector<double> zlist;

  xlist.push_back(x);
  ylist.push_back(y);
  zlist.push_back(z);

  double ex = 0.0;
  double ey = 0.0;
  double ez = 0.0;
  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;

  // Use NaN so the first StopHere call does not trigger the z==zPrevious guard.
  double zPrevious = std::numeric_limits<double>::quiet_NaN();

  constexpr int maxSteps = 20000;
  int step = 0;
  while (!StopHere(x, y, z, zPrevious) && step < maxSteps)
  {
    zPrevious = z;

    GetMagneticFieldTesla(x, y, z, bx, by, bz);
    GetElectricFieldVcm(x, y, z, ex, ey, ez);
    m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);

    if (!std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz) ||
        (vx == 0.0 && vy == 0.0 && vz == 0.0))
    {
      std::cerr << PHWHERE
                << " Garfield returned invalid/zero electron velocity; stopping drift. "
                << " p_tpc=(" << x << ", " << y << ", " << z << ") cm"
                << " E=(" << ex << ", " << ey << ", " << ez << ") V/cm"
                << " B=(" << bx << ", " << by << ", " << bz << ") T"
                << std::endl;
      break;
    }

    x = x - vx * step_ns;
    y = y - vy * step_ns;
    z = z - vz * step_ns;

    xlist.push_back(x);
    ylist.push_back(y);
    zlist.push_back(z);
    ++step;
  }

  if (step >= maxSteps)
  {
    std::cerr << PHWHERE << " ReverseDrift reached maxSteps=" << maxSteps
              << "; stopping drift at p_tpc=(" << x << ", " << y << ", " << z
              << ") cm" << std::endl;
  }

  // Keep all stored points, including the initial point.  This avoids returning
  // an empty TPolyLine3D if the starting point is already outside the active volume.
    // Drop the final point. It is the first point that triggered StopHere
  // on the next loop iteration, for example after crossing the central membrane.
  // This restores the old PHGarfield behavior.
  const size_t nGoodPoints = (xlist.size() > 1) ? xlist.size() - 1 : xlist.size();

  TPolyLine3D* poly = new TPolyLine3D(static_cast<int>(nGoodPoints));
  for (size_t i = 0; i < nGoodPoints; ++i)
  {
    poly->SetPoint(static_cast<int>(i), xlist[i], ylist[i], zlist[i]);
  }

  return poly;
}

TPolyLine3D* PHGarfield::ReverseDriftGlobalCoords(double x_cm, double y_cm, double z_cm, double step_ns)
{
  // Caller gives a GLOBAL point and receives a GLOBAL-coordinate polyline.
  // Actual Garfield drift is still done in local TPC coordinates.

  const TVector3 p_tpc = GlobalPointToTpcPoint(x_cm, y_cm, z_cm);

  TPolyLine3D* poly = ReverseDrift(p_tpc.X(), p_tpc.Y(), p_tpc.Z(), step_ns);
  if (!poly)
  {
    return nullptr;
  }

  const int nPoints = poly->GetN();
  Float_t* points = poly->GetP();

  for (int i = 0; i < nPoints; ++i)
  {
    const double px = points[3 * i + 0];
    const double py = points[3 * i + 1];
    const double pz = points[3 * i + 2];

    const TVector3 p_global = TpcPointToGlobalPoint(px, py, pz);
    poly->SetPoint(i, p_global.X(), p_global.Y(), p_global.Z());
  }

  return poly;
}


bool PHGarfield::StopHere(const double x, const double y, const double z,
                          const double zPrevious)
{
  const double r = std::hypot(x, y);

  if (r < 18.0)
  {
    return true;
  }
  if (r > 82.0)
  {
    return true;
  }
  if (z > 120.0)
  {
    return true;
  }
  if (z < -120.0)
  {
    return true;
  }

  // The first step passes zPrevious=NaN.  Do not apply crossing/stuck checks then.
  if (std::isfinite(zPrevious))
  {
    // z crossed the central membrane.
    if (z * zPrevious < 0.0)
    {
      return true;
    }

    // Protect against infinite loops if Garfield returns no longitudinal motion.
    if (z == zPrevious)
    {
      return true;
    }
  }

  return false;
}
