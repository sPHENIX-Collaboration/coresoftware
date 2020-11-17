// Allows the user to generate hijing events and store the results in
// a HepMC file.
//
// Inspired by code from ATLAS.  Thanks!
//
#include "HiMain1.h"
#include "HiMain2.h"
#include "HiParnt.h"
#include "HiStrng.h"
#include "HijCrdn.h"
#include "HijJet1.h"
#include "HijJet2.h"
#include "HijJet4.h"
#include "RanSeed.h"

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/IO_AsciiParticles.h>
#include <HepMC/IO_GenEvent.h>

#include <CLHEP/Geometry/Point3D.h>
#include <CLHEP/Random/MTwistEngine.h>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Vector/LorentzVector.h>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#define f2cFortran
#define gFortran
// cppcheck-suppress *
#include <cfortran.h>

PROTOCCALLSFSUB3(HIJING, hijing, STRING, FLOAT, FLOAT)
#define HIJING(FRAME, BMIN0, BMAX0) \
  CCALLSFSUB3(HIJING, hijing, STRING, FLOAT, FLOAT, FRAME, BMIN0, BMAX0)

PROTOCCALLSFSUB8(HIJSET, hijset, FLOAT, STRING, STRING, STRING, INT, INT, INT, INT)
#define HIJSET(EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)                      \
  CCALLSFSUB8(HIJSET, hijset, FLOAT, STRING, STRING, STRING, INT, INT, INT, INT, \
              EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)

using namespace std;
using namespace boost;
using namespace boost::property_tree;

void hijfst_control(int, vector<string>, vector<float>, vector<int>, vector<float>, vector<float>, vector<float>);

CLHEP::HepRandomEngine *engine;

float atl_ran(int *)
{
  return (float) CLHEP::RandFlat::shoot(engine);
}
// This prevents cppcheck to flag the next line as error
// cppcheck-suppress *
FCALLSCFUN1(FLOAT, atl_ran, ATL_RAN, atl_ran, PINT)

typedef HepGeom::Point3D<double> HepPoint3D;
typedef HepMC::GenEvent::particle_iterator piter;

int fillEvent(HepMC::GenEvent *evt);

// Accessor to HIJING Options and parameters COMMON
HiParnt m_hiparnt;

// Accessor to HIJING Random number seed COMMON
RanSeed m_ranseed;

// Accessors to HIJING Event information COMMONS
HiMain1 m_himain1;
HiMain2 m_himain2;
HijJet1 m_hijjet1;
HijJet2 m_hijjet2;
HijJet4 m_hijjet4;
HiStrng m_histrng;
HijCrdn m_hijcrdn;

float efrm;
std::string m_frame;
std::string m_proj;
std::string m_targ;
int iap;
int iat;
int izp;
int izt;
int spec;
double m_vertexOffsetCut = 1.0E-7;
bool keepSpectators;

int main(int argc, char **argv)
{
  string config_filename = "sHijing.xml";
  string output;
  long randomSeed = 0;
  bool randomseed_set = false;
  unsigned int N = 1;
  bool NEvents_set = false;
  for (int i = 1; i < argc; ++i)
  {
    std::string optionstring = argv[i];
    if (optionstring == "-h")
    {
      cout << endl
           << "Usage: sHijing <config xmlfile [sHijing.xml]>" << endl;
      cout << endl;
      cout << "Parameters:" << endl;
      cout << "-n <number of events [1]>" << endl;
      cout << "-o <outputfile [sHijing.dat]>" << endl;
      cout << "-s <random seet [std::random_device]>" << endl;
      exit(0);
    }
    else if (optionstring == "-o")
    {
      if (i + 1 < argc)
      {                      // Make sure we aren't at the end of argv!
        output = argv[++i];  // Increment 'i' so we get the argument
      }
      else
      {  // Uh-oh, there was no argument to the destination option.
        std::cerr << "-o option requires one argument." << std::endl;
        exit(1);
      }
      continue;
    }
    else if (optionstring == "-s")
    {
      if (i + 1 < argc)
      {                                     // Make sure we aren't at the end of argv!
        randomSeed = std::stol(argv[++i]);  // Increment 'i' so get the argument.
        randomseed_set = true;
      }
      else
      {  // Uh-oh, there was no argument to the destination option.
        std::cerr << "-s option requires one argument." << std::endl;
        exit(1);
      }
      continue;
    }
    else if (optionstring == "-n")
    {
      if (i + 1 < argc)
      {                             // Make sure we aren't at the end of argv!
        N = std::stoul(argv[++i]);  // Increment 'i' so get the argument.
        NEvents_set = true;
      }
      else
      {  // Uh-oh, there was no argument to the destination option.
        std::cerr << "-s option requires one argument." << std::endl;
        exit(1);
      }
      continue;
    }
    else
    {
      config_filename = argv[i];
    }
  }
  char frame[] = "        ";
  char proj[] = "        ";
  char targ[] = "        ";

  using boost::property_tree::ptree;
  iptree pt, null;
  std::ifstream config_file(config_filename);

  if (config_file)
  {
    // Read XML configuration file.
    read_xml(config_file, pt);
    cout << "using config file: " << config_filename << endl;
  }
  else
  {
    cout << "no xml config file - using internal values" << endl;
  }
  efrm = pt.get("HIJING.EFRM", 200.0);
  m_frame = pt.get("HIJING.FRAME", "CMS");
  m_proj = pt.get("HIJING.PROJ", "A");
  m_targ = pt.get("HIJING.TARG", "A");
  iap = pt.get("HIJING.IAP", 197);
  izp = pt.get("HIJING.IZP", 79);
  iat = pt.get("HIJING.IAT", 197);
  izt = pt.get("HIJING.IZT", 79);
  float bmin = pt.get("HIJING.BMIN", 0.0);
  float bmax = pt.get("HIJING.BMAX", 0.0);
  if (!NEvents_set)
  {
    N = pt.get("HIJING.N", 1);
  }
  keepSpectators = pt.get("HIJING.KEEP_SPECTATORS", 1);
  if (output.empty())
  {
    output = pt.get("HIJING.OUTPUT", "sHijing.dat");
  }
  if (!randomseed_set)
  {
    std::random_device rdev;
    randomSeed = pt.get("HIJING.RANDOM.SEED", rdev());
  }
  engine = new CLHEP::MTwistEngine(randomSeed);

  // See if there are any sections for HIPR1, IHPR2
  iptree &pr1 = pt.get_child("HIJING.HIPR1", null);
  BOOST_FOREACH (iptree::value_type &v, pr1)
  {
    int key = boost::lexical_cast<int>(v.first.data());
    float value = boost::lexical_cast<float>(v.second.data());
    m_hiparnt.hipr1(key) = value;
  }
  iptree &pr2 = pt.get_child("HIJING.IHPR2", null);
  BOOST_FOREACH (iptree::value_type &v, pr2)
  {
    int key = boost::lexical_cast<int>(v.first.data());
    int value = boost::lexical_cast<int>(v.second.data());
    m_hiparnt.ihpr2(key) = value;
  }

  int fastjet_enable_p = 0;
  std::vector<string> algorithm_v;
  std::vector<float> R_v;
  std::vector<int> PID_v;
  std::vector<float> EtaMin_v;
  std::vector<float> EtaMax_v;
  std::vector<float> EtMin_v;

  iptree &it = pt.get_child("HIJING.FASTJET", null);
  BOOST_FOREACH (iptree::value_type &v, it)
  {
    if (to_upper_copy(v.first) != "ALGORITHM") continue;
    algorithm_v.push_back(to_upper_copy(v.second.get("NAME", "ANTIKT")));
    R_v.push_back(v.second.get("R", 0.2));
    PID_v.push_back(v.second.get("PID", 2000000));

    EtaMin_v.push_back(v.second.get("EtaMin", -2));
    EtaMax_v.push_back(v.second.get("EtaMax", 2));
    EtMin_v.push_back(v.second.get("EtMin", 5));

    fastjet_enable_p = 1;
  }
  cout << "seed: " << randomSeed << endl;
  cout << "output: " << output << endl;
  cout << "Number of Events: " << N << endl;

  hijfst_control(fastjet_enable_p, algorithm_v, R_v, PID_v, EtaMin_v, EtaMax_v, EtMin_v);

  // The call to Hijing needs simple C-style strings.
  m_frame.copy(frame, m_frame.size());
  m_proj.copy(proj, m_proj.size());
  m_targ.copy(targ, m_targ.size());

  HIJSET(efrm, frame, proj, targ, iap, izp, iat, izt);

  //  int status;
  HepMC::IO_GenEvent ascii_io(output.c_str(), std::ios::out);
  unsigned int events = 0;
  do
  {
    HIJING(frame, bmin, bmax);

    if (m_himain1.ierrstat() != 0)
    {
      continue;
    }
    events++;

    HepMC::GenEvent *evt = new HepMC::GenEvent();
    evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);
    evt->set_event_number(events);

    fillEvent(evt);

    ascii_io << evt;

    delete evt;
  } while (events < N);

  return 0;
}

int fillEvent(HepMC::GenEvent *evt)
{
  int np = m_himain1.np();
  int nt = m_himain1.nt();
  int n0 = m_himain1.n0();
  int n01 = m_himain1.n01();
  int n10 = m_himain1.n10();
  int n11 = m_himain1.n11();

  // The Hijing documentation is a little unclear about this, but the
  // four values, n0, n01, n10 and n11, count collisions between
  // nucleons that have already (or have not already had) collisions
  // with other nucleons.  N_coll, the way we use it, is the sum of
  // all these values.
  int ncoll = n0 + n01 + n10 + n11;

  int natt = m_himain1.natt();
  float b = m_hiparnt.hint1(19);
  float bphi = m_hiparnt.hint1(20);

  // Determine the participant eccentricity.
  std::vector<float> x;
  std::vector<float> y;
  float tx, ty;
  float bbx = b * cos(bphi);
  float bby = b * sin(bphi);
  // std::cout << "x,y,c" << std::endl;
  for (int i = 0; i < m_hiparnt.ihnt2(1); i++)
  {
    if (m_histrng.nfp(i, 5) > 2)
    {
      tx = m_hijcrdn.yp(1, i) + bbx;
      ty = m_hijcrdn.yp(2, i) + bby;
      x.push_back(tx);
      y.push_back(ty);

      // std::cout << tx << "," << ty << "," << "0" << std::endl;
    }
  }
  for (int i = 0; i < m_hiparnt.ihnt2(3); i++)
  {
    if (m_histrng.nft(i, 5) > 2)
    {
      tx = m_hijcrdn.yt(1, i);
      ty = m_hijcrdn.yt(2, i);
      x.push_back(tx);
      y.push_back(ty);

      // std::cout << tx << "," << ty << "," << "1" << std::endl;
    }
  }

  float e_part = 0.0;
  if (x.size() != 0)
  {
    float N = x.size();
    float xa = std::accumulate(x.begin(), x.end(), 0) / N;
    float ya = std::accumulate(y.begin(), y.end(), 0) / N;

    float sxx = 0;
    float syy = 0;
    float sxy = 0;
    for (int i = 0; i < N; i++)
    {
      sxx += (x[i] - xa) * (x[i] - xa);
      syy += (y[i] - ya) * (y[i] - ya);
      sxy += (x[i] - xa) * (y[i] - ya);
    }
    sxx /= N;
    syy /= N;
    sxy /= N;
    e_part = std::sqrt((sxx - syy) * (sxx - syy) + 4 * sxy * sxy) / (sxx + syy);
  }

  // Need to calculate a few things missing from this list
  HepMC::HeavyIon hi(1, np, nt, ncoll, 0, 0, n01, n10, n11, b, bphi, e_part, 42.0);

  evt->set_heavy_ion(hi);

  // Set the random generator seeds.  How do people handle the fact
  // that CLHEP produced unsigned longs and HepMC wants signed ones?
  // evt->set_random_states(engine->put());

  //  Did we keep decay history?
  bool keptHistory = (m_hiparnt.ihpr2(21) == 1);

  //  The number of particles in the hijing output
  int numHijingPart = m_himain1.natt();

  // Vectors that will keep track of where particles originate from and die
  std::vector<int> partOriginVertex_vec;
  std::vector<int> partDecayVertex_vec;
  std::vector<HepMC::GenParticle *> particleHepPartPtr_vec;

  partOriginVertex_vec.assign(numHijingPart, 0);
  partDecayVertex_vec.assign(numHijingPart, -1);
  particleHepPartPtr_vec.assign(numHijingPart, (HepMC::GenParticle *) 0);

  // Vector that will keep pointers to generated vertices
  std::vector<HepMC::GenVertex *> vertexPtrVec;

  // Create the event vertex
  CLHEP::HepLorentzVector newVertex;
  newVertex = CLHEP::HepLorentzVector(0., 0., 0., 0.);
  HepMC::GenVertex *v1 = new HepMC::GenVertex(newVertex);

  evt->add_vertex(v1);
  vertexPtrVec.push_back(v1);

  double eproj = (m_frame == "CMS") ? efrm / 2.0 : efrm;
  int proj_id = (m_proj == "A") ? 3000000 + iap : 2212;
  v1->add_particle_in(new HepMC::GenParticle(CLHEP::HepLorentzVector(0., 0., eproj, eproj),
                                             proj_id, 101));

  double etarg = (m_frame == "CMS") ? efrm / 2.0 : efrm;
  int targ_id = (m_targ == "A") ? 3000000 + iap : 2212;
  v1->add_particle_in(new HepMC::GenParticle(CLHEP::HepLorentzVector(0., 0., -etarg, etarg),
                                             targ_id, 102));

  // Loop on all Hijing particles and put them all as outgoing from the event vertex
  for (int i = 1; i <= natt; ++i)
  {
    if (m_himain2.katt(i, 4) == 103)
    {
      // We handle jets a little differently.
      CLHEP::HepLorentzVector jetP4(m_himain2.patt(i, 1),
                                    m_himain2.patt(i, 2),
                                    m_himain2.patt(i, 3),
                                    m_himain2.patt(i, 4));
      v1->add_particle_out(new HepMC::GenParticle(jetP4,
                                                  m_himain2.katt(i, 1),
                                                  m_himain2.katt(i, 4)));

      continue;
    }

    //  Skip non-interacting projectile and target nucleons
    if (!keepSpectators &&
        ((m_himain2.katt(i, 2) == 0) || (m_himain2.katt(i, 2)) == 10))
      continue;

    //  Find the vertex of the parent particle
    int parentIndex = m_himain2.katt(i, 3) - 1;
    int parentOriginIndex = 0;
    int parentDecayIndex = -1;

    //  If the particle has a true parent, record its vertex info
    if (parentIndex >= 0)
    {
      parentOriginIndex = partOriginVertex_vec[parentIndex];
      parentDecayIndex = partDecayVertex_vec[parentIndex];
    }

    //  A CLHEP vector containing the particle start point
    CLHEP::HepLorentzVector particleStart(m_himain2.vatt(i, 1),
                                          m_himain2.vatt(i, 2),
                                          m_himain2.vatt(i, 3),
                                          m_himain2.vatt(i, 4));

    //  By default, the particle originates from the primary vertex
    int particleVertexIndex = 0;

    //  Have we kept the history?
    if (keptHistory)
    {
      //  Check to see if we've already generated a decay vertex
      //  for this parent
      if (parentDecayIndex != -1)
      {
        //  Make sure it is consistent with this track origin
        HepPoint3D vertex_pos(vertexPtrVec[parentDecayIndex]->point3d().x(),
                              vertexPtrVec[parentDecayIndex]->point3d().y(),
                              vertexPtrVec[parentDecayIndex]->point3d().z());
        double distance = vertex_pos.distance(particleStart.vect());

        if (distance > m_vertexOffsetCut)
        {
          HepMC::GenVertex::particles_out_const_iterator iter;
          for (iter = vertexPtrVec[parentDecayIndex]->particles_out_const_begin();
               iter != vertexPtrVec[parentDecayIndex]->particles_out_const_end();
               iter++)
          {
            std::cout << (*iter)->barcode() << ", ";
          }

          std::cout << std::endl;
        }

        //  Nonetheless set the parent decay vertex to be this particle's vertex
        particleVertexIndex = parentDecayIndex;
      }
      else
      {
        //  Now compare the distance between the vertex FROM
        //  which the parent originates and the start of this
        //  particle
        HepPoint3D vertex_pos(vertexPtrVec[parentOriginIndex]->point3d().x(),
                              vertexPtrVec[parentOriginIndex]->point3d().y(),
                              vertexPtrVec[parentOriginIndex]->point3d().z());
        double distance = vertex_pos.distance(particleStart.vect());

        if (distance > m_vertexOffsetCut && parentIndex == -1)
        {
          //  *** Explicitly handle Hijing bug which generates
          //  *** particle with displaced vertex but no parent
          //  *** ????

          std::cout
              << "HIJING BUG:: Particle found with displaced vertex but no parent "
              << std::endl;
        }
        if (distance > m_vertexOffsetCut || parentOriginIndex != 0)
        {
          // We need to create a new vertex
          HepMC::GenVertex *newVertex_p = new HepMC::GenVertex(particleStart);

          evt->add_vertex(newVertex_p);
          vertexPtrVec.push_back(newVertex_p);
          particleVertexIndex = vertexPtrVec.size() - 1;

          //  Now we indicate that the parent has a decay vertex
          partDecayVertex_vec[parentIndex] = particleVertexIndex;

          //  Now tell the vertex about the particle that created it
          newVertex_p->add_particle_in(particleHepPartPtr_vec[parentIndex]);
        }
        else
        {
          //  Assign the particle to its parent's vertex
          particleVertexIndex = parentOriginIndex;
        }
      }
    }
    else
    {
      //  We have to brute-force search for a vertex that might match this particle
      int foundVert = -1;
      for (unsigned int ivert = 0; ivert < vertexPtrVec.size(); ivert++)
      {
        HepPoint3D vertex_pos(vertexPtrVec[ivert]->point3d().x(),
                              vertexPtrVec[ivert]->point3d().y(),
                              vertexPtrVec[ivert]->point3d().z());
        double distance = vertex_pos.distance(particleStart.vect());
        if (distance < m_vertexOffsetCut)
        {
          foundVert = ivert;
          break;
        }
      }

      if (foundVert == -1)
      {
        //  We need to create a new vertex
        HepMC::GenVertex *newVertex_p = new HepMC::GenVertex(particleStart);

        evt->add_vertex(newVertex_p);
        vertexPtrVec.push_back(newVertex_p);
        particleVertexIndex = vertexPtrVec.size() - 1;
      }
      else
      {
        particleVertexIndex = foundVert;
      }
    }

    //  If the Hijing particle has decayed, set the status appropriately
    int particleId = m_himain2.katt(i, 1);
    int particleStatus = 1;
    if (m_himain2.katt(i, 4) == 11)
    {
      particleStatus = 2;
    }

    // Create the new particle
    CLHEP::HepLorentzVector particleP4(m_himain2.patt(i, 1),
                                       m_himain2.patt(i, 2),
                                       m_himain2.patt(i, 3),
                                       m_himain2.patt(i, 4));

    HepMC::GenParticle *newParticle_p = new HepMC::GenParticle(particleP4,
                                                               particleId,
                                                               particleStatus);

    //  Record the particle in the vector of pointers (ostensibly we
    //  only need this when we have the history but for simplicity
    //  always do it)
    particleHepPartPtr_vec[i - 1] = newParticle_p;

    // Now add the particle to its vertex
    vertexPtrVec[particleVertexIndex]->add_particle_out(newParticle_p);
    partOriginVertex_vec[i - 1] = particleVertexIndex;
  }

  return 0;
}
