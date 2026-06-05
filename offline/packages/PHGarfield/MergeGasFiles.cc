#include <filesystem>
#include <iostream>
#include <string>
#include <TRandom3.h>
#include <ctime>
#include <TMath.h>

#include <TNtuple.h>
#include <TFile.h>

#include "Garfield/MediumMagboltz.hh"

using namespace std;

int main()
{

  //  This is a utility to test whether the process of "merging" files is actually different from a single file.
  //  It may be of no further use afterthe development was complete.
  //                                  TKH 6/2/2026
  int nValid=10000;
  TNtuple *Validity = new TNtuple("Validity", "Validity", "Valid:e:b:a:Vxerr:Vyerr:Vzerr");

  // New version chooses to NOT write output to a file (which seems broken),
  // but to instead just tries to merge the files and validate the copy in memory.
  const std::string dir = "gasfiles";
  const std::string out = "Ar75_CF20_iso5.gas";

  auto filename = [&](const int i)
  {
    return dir + "/PART_" + std::to_string(i) + ".gas";
  };

  Garfield::MediumMagboltz gas;
  Garfield::MediumMagboltz gas0;

  const std::string first = filename(0);
  if (!std::filesystem::exists(first))
    {
      std::cerr << "Missing first gas file: " << first << std::endl;
      return 1;
    }

  if (!gas.LoadGasFile(first))
    {
      std::cerr << "Failed to load " << first << std::endl;
      return 1;
    }

  // Gas 0 only loads the FIRST file.  This will test the memory validity of the merge...
  if (!gas0.LoadGasFile(first))
    {
      std::cerr << "Failed to load " << first << std::endl;
      return 1;
    }

  int nFiles = 1;
  for (int i = 1; ; ++i)
    {
      const std::string file = filename(i);
      
      if (!std::filesystem::exists(file))
	{
	  std::cout << "Stopping at first missing file: " << file << std::endl;
	  break;
	}
      
      std::cout << "Merging " << file << std::endl;
      
      if (!gas.MergeGasFile(file, true))
	{
	  std::cerr << "Failed to merge " << file << std::endl;
	  return 1;
	}
      
      ++nFiles;
    }

  //  Don't write out since it crashes?
  //gas.WriteGasFile("test.gas");

  // Now perform the validation test...
  double emin=400;
  double emax=400;
  //double ne=1;

  double bmin=1.15;
  double bmax=1.45;
  double nb=50;

  double amin=0.0;
  double amax=0.2;
  //double na=50;

  // Initialize using the current system time
  TRandom3 Randy(time(0));  //new initialization each run
  cout << endl << endl << "Valid Calls: "<< endl;
  for (int i=0; i<nValid; i++)
    {
      double eMag   = Randy.Uniform(emin, emax);
      double bMag   = Randy.Uniform(bmin, bmin+0.2*(bmax-bmin)/nb);  // Comes from file0...
      double a      = Randy.Uniform(amin, amax);
      double PHI    = Randy.Uniform(0.0,2.0*TMath::Pi());

      double ex     = 0;
      double ey     = 0;
      double ez     = eMag;

      double bx     = bMag * sin(a) * cos(PHI);
      double by     = bMag * sin(a) * sin(PHI);
      double bz     = bMag * cos(a);

      double vx;      
      double vy;      
      double vz;      

      double vx0;      
      double vy0;      
      double vz0;      

      gas.ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      gas0.ElectronVelocity(ex, ey, ez, bx, by, bz, vx0, vy0, vz0);

      /*
      cout << "  i:" <<i;
      cout << "  Ex:"<<ex;
      cout << "  Ey:"<<ey;
      cout << "  Ez:"<<ez;
      cout << "  Bx:"<<bx;
      cout << "  By:"<<by;
      cout << "  Bz:"<<bz;
      cout << "  Vx:"<<vx;
      cout << "  Vy:"<<vy;
      cout << "  Vz:"<<vz;
      cout << endl;
      cout << "  i:" <<i;
      cout << "  Ex:"<<ex;
      cout << "  Ey:"<<ey;
      cout << "  Ez:"<<ez;
      cout << "  Bx:"<<bx;
      cout << "  By:"<<by;
      cout << "  Bz:"<<bz;
      cout << "  Vx0:"<<vx0;
      cout << "  Vy0:"<<vy0;
      cout << "  Vz0:"<<vz0;
      cout << "  Vx%:"<<(vx-vx0)/vx;
      cout << "  Vy%:"<<(vy-vy0)/vy;
      cout << "  Vy%:"<<(vz-vz0)/vz;
      cout << endl;
      cout << endl;
      */
      
      double DelVx = (vx-vx0)/((vx+vx0)/2.);
      double DelVy = (vy-vy0)/((vy+vy0)/2.);
      double DelVz = (vz-vz0)/((vz+vz0)/2.);
      Validity->Fill(1,sqrt(ex*ex + ey*ey + ez*ez),sqrt(bx*bx + by*by + bz*bz),a,DelVx, DelVy, DelVz);
    }
  
  cout << endl << endl << "Invalid Calls: "<< endl;
  for (int i=0; i<nValid; i++)
    {
      double eMag   = Randy.Uniform(emin, emax);
      double bMag   = Randy.Uniform(bmin+5.0*(bmax-bmin)/nb, bmax);  // Comes from beyond file0...
      double a      = Randy.Uniform(amin, amax);
      double PHI    = Randy.Uniform(0.0,2.0*TMath::Pi());

      double ex     = 0;
      double ey     = 0;
      double ez     = eMag;

      double bx     = bMag * sin(a) * cos(PHI);
      double by     = bMag * sin(a) * sin(PHI);
      double bz     = bMag * cos(a);

      double vx;      
      double vy;      
      double vz;      
      double vx0;      
      double vy0;      
      double vz0;      
      gas.ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      gas0.ElectronVelocity(ex, ey, ez, bx, by, bz, vx0, vy0, vz0);
      /*
      cout << "  i:" <<i;
      cout << "  Ex:"<<ex;
      cout << "  Ey:"<<ey;
      cout << "  Ez:"<<ez;
      cout << "  Bx:"<<bx;
      cout << "  By:"<<by;
      cout << "  Bz:"<<bz;
      cout << "  Vx:"<<vx;
      cout << "  Vy:"<<vy;
      cout << "  Vz:"<<vz;
      cout << endl;
      cout << "  i:" <<i;
      cout << "  Ex:"<<ex;
      cout << "  Ey:"<<ey;
      cout << "  Ez:"<<ez;
      cout << "  Bx:"<<bx;
      cout << "  By:"<<by;
      cout << "  Bz:"<<bz;
      cout << "  Vx0:"<<vx0;
      cout << "  Vy0:"<<vy0;
      cout << "  Vz0:"<<vz0;
      cout << "  Vx%:"<<(vx-vx0)/vx;
      cout << "  Vy%:"<<(vy-vy0)/vy;
      cout << "  Vy%:"<<(vz-vz0)/vz;
      cout << endl;
      cout << endl;
      */
      double DelVx = (vx-vx0)/((vx+vx0)/2.);
      double DelVy = (vy-vy0)/((vy+vy0)/2.);
      double DelVz = (vz-vz0)/((vz+vz0)/2.);
      Validity->Fill(0,sqrt(ex*ex + ey*ey + ez*ez),sqrt(bx*bx + by*by + bz*bz),a,DelVx, DelVy, DelVz);
    }
  

  TFile *output= new TFile("GarfieldValidity.root","RECREATE");
  Validity->Write();
  output->Close();
  
  return 0;
}
