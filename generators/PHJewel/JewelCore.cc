#include "JewelCore.h"

#include <dlfcn.h>

#include <initializer_list>
#include <stdexcept>
#include <string>

template <class F>
void JewelCore::sym(F &f, const char *name)
{
  f = reinterpret_cast<F>(dlsym(m_Handle, name));
  if (!f)
  {
    throw std::runtime_error(std::string("JewelCore: ") + m_Lib + " lacks " + name);
  }
}

void JewelCore::load(const std::string &lib)
{
  // RTLD_DEEPBIND: the core resolves its own undefined symbols (pdfset_,
  // structm_ from LHAPDF 6) in its dependencies first. Without it they bind
  // to whatever is already global -- in a G4 job that is LHAPDF 5.9.1,
  // loaded by libPythia6 for the G4 decayer, and JEWEL stops in PDFSET.
  m_Handle = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL | RTLD_DEEPBIND);
  if (!m_Handle)
  {
    throw std::runtime_error("JewelCore: dlopen(" + lib + ") failed: " + dlerror());
  }
  sym(f_abi, "jw_abi_version");
  sym(f_err, "jw_last_error");
  sym(f_init, "jw_init");
  sym(f_chan, "jw_set_channel");
  sym(f_next, "jw_next");
  sym(f_text, "jw_event_text");
  sym(f_info, "jw_info");
  sym(f_counts, "jw_counts");
  sym(f_config, "jw_config");
  sym(f_finish, "jw_finish");
  if (f_abi() != kAbiVersion)
  {
    throw std::runtime_error("JewelCore: " + lib + " has ABI " + std::to_string(f_abi()) +
                             ", this code expects " + std::to_string(kAbiVersion) + "; rebuild the package");
  }
  m_Lib = lib;
}

void JewelCore::init(const std::string &card, const std::string &log, int seed, const std::string &hepmc,
                     const std::string &evinfo)
{
  if (!m_Handle)
  {
    throw std::runtime_error("JewelCore::init before load");
  }
  if (f_init(card.c_str(), log.c_str(), hepmc.c_str(), evinfo.c_str(), seed) != 0)
  {
    throw std::runtime_error(std::string("JewelCore: jw_init: ") + f_err());
  }
  m_Seed = seed;
  m_Config = readConfig();
  m_Rng.seed(static_cast<unsigned long long>(seed) * 2654435761ULL + 12345ULL);
  const double A = m_Config.mass;
  const double Z = m_Config.nproton;
  const double N = A - Z;
  m_Prob = {Z * Z / (A * A), Z * N / (A * A), N * Z / (A * A), N * N / (A * A)};
  m_Init = true;
}

const char *JewelCore::next()
{
  if (!m_Init)
  {
    throw std::runtime_error("JewelCore::next before init");
  }
  if (m_Channel == 0 || (m_Config.iso == 0 && m_InBlock >= m_Block))
  {
    chooseChannel();
  }
  long discards = 0;
  while (true)
  {
    ++m_Pass;
    if (f_next(static_cast<int>(m_Pass)) == 1)
    {
      break;
    }
    ++m_NDiscardedPasses;
    if (m_MaxDiscards > 0 && ++discards >= m_MaxDiscards)
    {
      throw std::runtime_error("JewelCore: " + std::to_string(discards) + " consecutive passes discarded by JEWEL");
    }
  }
  m_Last = readInfo();
  ++m_InBlock;
  ++m_SessN;
  ++m_NAccepted;
  m_SumW += m_Last.weight;
  m_ChanN[m_Channel - 1] += 1;
  m_ChanW[m_Channel - 1] += m_Last.weight;
  m_SessXsec = m_Last.xsec_mb;
  return f_text();
}

std::string JewelCore::event_stream() const
{
  return std::string("HepMC::Version 2.06.05\nHepMC::IO_GenEvent-START_EVENT_LISTING\n") + f_text() +
         "HepMC::IO_GenEvent-END_EVENT_LISTING\n";
}

void JewelCore::finish()
{
  if (m_Init && !m_Finished)
  {
    commitSession();
    f_finish();
    m_Finished = true;
  }
}

double JewelCore::xsec_mb() const
{
  const double sn = m_XsecN + m_SessXsec * m_SessN;
  const double n = m_XsecCount + m_SessN;
  return n > 0 ? sn / n : 0.0;
}

JewelCore::Counts JewelCore::counts() const
{
  double d[16] = {0};
  f_counts(d);
  Counts c;
  c.ngood = d[0];
  c.ndisc = d[1];
  c.nstrange = d[2];
  c.nevdone = d[3];
  c.sumw = d[4];
  c.wdisc = d[5];
  c.ntotspliti = d[6];
  c.noverspliti = d[7];
  c.ntotpdf = d[8];
  c.noverpdf = d[9];
  c.ntotxsec = d[10];
  c.noverxsec = d[11];
  c.ntotsuda = d[12];
  c.noversuda = d[13];
  return c;
}

JewelCore::Config JewelCore::readConfig() const
{
  int ic[16] = {0};
  double dc[8] = {0};
  f_config(ic, dc);
  Config c;
  c.mass = ic[0];
  c.nproton = ic[1];
  c.iso = ic[2];
  // PROCESS packed as four characters, most significant byte first
  for (const unsigned int shift : {24U, 16U, 8U, 0U})
  {
    c.process += static_cast<char>((static_cast<unsigned int>(ic[3]) >> shift) & 0xffU);
  }
  c.weighted = ic[4];
  c.keeprecoils = ic[5];
  c.writescatcen = ic[6];
  c.writedummies = ic[7];
  c.dosubtraction = ic[8];
  c.hadro = ic[9];
  c.shorthepmc = ic[10];
  c.isrscat = ic[11];
  c.mpi = ic[12];
  c.nf = ic[13];
  c.sqrts = dc[0];
  c.ptmin = dc[1];
  c.ptmax = dc[2];
  c.etamax = dc[3];
  c.q0 = dc[4];
  c.lqcd = dc[5];
  c.lps = dc[6];
  return c;
}

JewelCore::EventInfo JewelCore::readInfo() const
{
  double d[16] = {0};
  f_info(d);
  EventInfo e;
  e.weight = d[0];
  e.xsec_mb = d[1];
  e.b = d[2];
  e.centrality = d[3];
  e.isub = static_cast<int>(d[4]);
  e.id1 = static_cast<int>(d[5]);
  e.id2 = static_cast<int>(d[6]);
  e.x1 = d[7];
  e.x2 = d[8];
  e.pthat = d[9];
  e.Q2 = d[10];
  e.ngood = static_cast<long>(d[11]);
  e.nscatcen = static_cast<int>(d[12]);
  e.channel = static_cast<int>(d[13]);
  return e;
}

void JewelCore::commitSession()
{
  m_XsecN += m_SessXsec * m_SessN;
  m_XsecCount += m_SessN;
  m_SessN = 0;
}

void JewelCore::chooseChannel()
{
  int ich = m_Config.iso;  // 1..4 fixed
  if (ich == 0)
  {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double r = u(m_Rng);
    double acc = 0;
    ich = 4;
    for (int i = 0; i < 4; ++i)
    {
      acc += m_Prob[i];
      if (r < acc)
      {
        ich = i + 1;
        break;
      }
    }
  }
  if (ich < 1 || ich > 4)
  {
    throw std::runtime_error("JewelCore: JEWEL reports an unknown ISOCHANNEL");
  }
  m_InBlock = 0;
  if (ich == m_Channel)
  {
    return;  // same channel drawn again: PYTHIA keeps running
  }
  commitSession();
  if (f_chan(ich) != 0)
  {
    throw std::runtime_error(std::string("JewelCore: ") + f_err());
  }
  m_Channel = ich;
  ++m_Switches;
}
