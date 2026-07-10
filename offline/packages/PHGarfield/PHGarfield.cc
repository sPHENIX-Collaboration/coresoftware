#include "PHGarfield.h"
#include <phool/phool.h>
#include <cdbobjects/CDBTTree.h>

#include <phfield/PHField3DCartesian.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TPolyLine3D.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <Garfield/ComponentUser.hh>
#include <Garfield/MediumMagboltz.hh>

#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>  // for basic_ostream, operat...
#include <vector>
#include <regex>
#include <string>
#include <utility>
#include <sstream>
#include <iomanip>

namespace fs = std::filesystem;

PHGarfield::PHGarfield(const std::string& name)
  : SubsysReco(name),
    //m_defaultGasfile("/sphenix/user/hemmick/gasfiles_20260624/Ar75_CF20_iso5.gas")
    m_defaultGasfile("/sphenix/user/hemmick/gasfiles_20260624")
{
}

PHGarfield::~PHGarfield()
{
  //  Housekeeping.
  delete m_field;
  delete m_cdbTPCMAPttree;
  delete m_component;
  delete m_gas;
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

void PHGarfield::GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double& bx_t, double& by_t, double& bz_t) const
{
  // NOTE:  Garfield uses  cm, V/cm, and Tesla.
  //        CLHEP    uses  mm, V/mm, and kiloTesla
  //        PHField3DCartesian follows the CLHEP conventions for magnetic fields.

  double point[4] =
      {
          x_cm * CLHEP::cm,
          y_cm * CLHEP::cm,
          z_cm * CLHEP::cm,
          //(z_cm-20.0) * CLHEP::cm,
          0.0};

  double bfield[3] = {0.0, 0.0, 0.0};

  //  Get the magnetic field via the PHField3DCartesian object constructed usinf the CDB url reference.
  m_field->GetFieldValue(point, bfield);

  bx_t = bfield[0] / CLHEP::tesla;
  by_t = bfield[1] / CLHEP::tesla;
  bz_t = bfield[2] / CLHEP::tesla;
}

void PHGarfield::GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm) const
{
  // NOTE:  Garfield uses  cm, V/cm, and Tesla.
  (void) x_cm;
  (void) y_cm;

  ex_vcm = 0.0;
  ey_vcm = 0.0;
  ez_vcm = z_cm > 0 ? -400.0 : 400.0;
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
  std::vector<double> xlist;
  std::vector<double> ylist;
  std::vector<double> zlist;

  xlist.push_back(x);
  ylist.push_back(y);
  zlist.push_back(z);

  double ex;
  double ey;
  double ez;
  double bx;
  double by;
  double bz;
  double vx;
  double vy;
  double vz;

  double zPrevious = z;
  while (!StopHere(x, y, z, zPrevious))
  {
    zPrevious = z;
    GetMagneticFieldTesla(x, y, z, bx, by, bz);
    GetElectricFieldVcm(x, y, z, ex, ey, ez);
    m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);

    x = x - vx * step_ns;
    y = y - vy * step_ns;
    z = z - vz * step_ns;

    xlist.push_back(x);
    ylist.push_back(y);
    zlist.push_back(z);
  }

  TPolyLine3D* poly = new TPolyLine3D(xlist.size() - 1);
  for (unsigned int i = 0; i < xlist.size() - 1; i++)
  {
    poly->SetPoint(i, xlist[i], ylist[i], zlist[i]);
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

  // z crossed the central membrane.
  if (z * zPrevious < 0.0)
  {
    return true;
  }

  return false;
}
