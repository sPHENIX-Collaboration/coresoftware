#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/IO_BaseClass.h>  // for IO_BaseClass
#include <HepMC/IO_GenEvent.h>
#include <HepMC/SimpleVector.h>  // for FourVector
#include <HepMC/Units.h>         // for GEV, MM

#include <CLHEP/Vector/LorentzVector.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>  // for TParticlePDG

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

TDatabasePDG* PDGdb;

// Function to parse the STARlight output file and return a vector of HepMC::GenParticle pointers
// std::vector<HepMC::GenParticle*> parseStarlightOutput(const std::string& filename)
int fillEvent(HepMC::GenEvent* evt, std::ifstream& file)
{
  std::string label;

  if (file >> label)
  {
    // first line should be the event
    // EVENT: n ntracks nvertices ,
    if (label != "EVENT:")
    {
      std::cout << "label is " << label << std::endl;
      return -1;
    }

    int nevt, ntrk, nvtx;
    file >> nevt >> ntrk >> nvtx;
    if (nevt % 100 == 0)
    {
      std::cout << nevt << std::endl;
    }

    evt->set_event_number(nevt);

    // Next line should be the vertex
    // Note: starlight currently only has one vertex, but a future version could have more
    // VERTEX: x y z t nv nproc nparent ndaughters
    double x, y, z, t;
    int nv, nproc, npar, ndau;
    file >> label >> x >> y >> z >> t >> nv >> nproc >> npar >> ndau;
    if (label != "VERTEX:")
    {
      return -1;
    }
    CLHEP::HepLorentzVector newVertex;
    newVertex = CLHEP::HepLorentzVector(x, y, z, t);
    HepMC::GenVertex* v0 = new HepMC::GenVertex(newVertex);
    evt->add_vertex(v0);

    // TRACK: GPID px py py nev ntr stopv PDGPID
    double px, py, pz;  // three vector components of the track's momentum
    int gpid, nev, ntr, stopv, pdgpid;
    for (int itrk = 0; itrk < ndau; itrk++)
    {
      file >> label >> gpid >> px >> py >> pz >> nev >> ntr >> stopv >> pdgpid;
      if (label != "TRACK:")
      {
        return -1;
      }

      // Get mass of particle
      TParticlePDG* partinfo = PDGdb->GetParticle(pdgpid);
      double mass = partinfo->Mass();

      // Calculate energy
      double E = sqrt(mass * mass + px * px + py * py + pz * pz);

      CLHEP::HepLorentzVector p4(px, py, pz, E);

      const int particle_status = 1;

      // Create a HepMC::GenParticle using status code 1 for stable particle
      v0->add_particle_out(new HepMC::GenParticle(p4, pdgpid, particle_status));
    }
  }
  else
  {
    // end of file, or failure to read file
    return -2;
  }

  return 1;
}

// Main function to convert STARlight output to HEPMC format
int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <HepMC_output_file>" << std::endl;
    return 1;
  }

  // std::string starlightfname = "slight.out"; // Slight Input file Hardcoded.
  // This way the starlight input file is not hardcoded which is bad.
  std::string starlightfname = argv[1];
  std::string hepmcfname = argv[2];

  // Open the STARlight output file
  std::ifstream starlightfile(starlightfname);
  if (!starlightfile.is_open())
  {
    std::cerr << " STARlight file not found: " << starlightfname << std::endl;
    return -1;
  }

  // Open pdg database
  PDGdb = TDatabasePDG::Instance();

  // Open the HepMC output file
  std::ofstream outputFile(hepmcfname);
  HepMC::IO_GenEvent ascii_io(hepmcfname.c_str(), std::ios::out);

//  unsigned int events = 0;

  while (true)
  {
//    events++;

    HepMC::GenEvent* evt = new HepMC::GenEvent();
    evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

    int status = fillEvent(evt, starlightfile);
    if (status < 0)
    {
      delete evt;
      break;
    }

    ascii_io << evt;

    delete evt;
  }

  std::cout << "Conversion completed successfully. HEPMC file saved to " << hepmcfname << std::endl;

  starlightfile.close();

  return 0;
}
