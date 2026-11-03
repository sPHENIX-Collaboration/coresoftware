#include "JewelRun.h"

#include "JewelCore.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

JewelCard::Result startJewel(JewelCore &core, const JewelRunConfig &c, std::ostream &os)
{
  JewelCard::Options opt;
  opt.medium = c.medium;
  opt.allowTableGeneration = c.allowTableGeneration;
  JewelCard::Result r = JewelCard::validate(c.card, opt);
  for (const std::string &n : r.notes)
  {
    os << "PHJewel card note: " << n << "\n";
  }
  if (!r.ok)
  {
    std::string msg = "PHJewel: card " + c.card + " (medium " + c.medium + ") refused:";
    for (const std::string &e : r.errors)
    {
      msg += "\n  - " + e;
    }
    throw std::runtime_error(msg);
  }
  std::string err;
  const std::string canon = JewelCard::writeCanonical(r, c.stem, JewelCard::absPath(c.card, "."), c.medium, err);
  if (canon.empty())
  {
    throw std::runtime_error("PHJewel: " + err);
  }
  os << "PHJewel: canonical card " << canon << " (hash " << r.hash << ")\n";

  core.load(c.coreLib.empty() ? "libjewel_" + c.medium + ".so" : c.coreLib);
  core.set_isospin_block(c.isospinBlock);
  core.init(canon, c.stem + ".log", c.seed, c.fortranHepMC, c.evinfo);

  // JEWEL never reports a value it failed to read; compare its state with ours
  const JewelCore::Config &j = core.config();
  const std::map<std::string, std::string> &m = r.main;
  std::string bad;
  auto chkI = [&](const char *k, long have)
  {
    if (have != JewelCard::I(m, k))
    {
      bad += std::string(" ") + k + "(card " + m.at(k) + ", JEWEL " + std::to_string(have) + ")";
    }
  };
  auto chkR = [&](const char *k, double have)
  {
    const double want = JewelCard::R(m, k);
    // JEWEL's defaults are REAL*4 literals: agree to single precision
    if (std::fabs(have - want) > 1e-6 * std::max(1.0, std::fabs(want)))
    {
      bad += std::string(" ") + k + "(card " + m.at(k) + ", JEWEL " + std::to_string(have) + ")";
    }
  };
  auto chkL = [&](const char *k, bool have)
  {
    if (have != JewelCard::L(m, k))
    {
      bad += std::string(" ") + k + "(card " + m.at(k) + ", JEWEL " + (have ? "T" : "F") + ")";
    }
  };
  chkI("MASS", j.mass);
  chkI("NPROTON", j.nproton);
  chkI("NF", j.nf);
  chkR("SQRTS", j.sqrts);
  chkR("PTMIN", j.ptmin);
  chkR("PTMAX", j.ptmax);
  chkR("ETAMAX", j.etamax);
  chkR("Q0", j.q0);
  chkR("LAMBDAQCD", j.lqcd);
  chkR("LAMBDAQCDPS", j.lps);
  chkL("WEIGHTED", j.weighted);
  chkL("KEEPRECOILS", j.keeprecoils);
  chkL("WRITESCATCEN", j.writescatcen);
  chkL("WRITEDUMMIES", j.writedummies);
  chkL("DOSUBTRACTION", j.dosubtraction);
  chkL("HADRO", j.hadro);
  chkL("SHORTHEPMC", j.shorthepmc);
  chkL("ISRSCAT", j.isrscat);
  chkL("MPI", j.mpi);
  if (j.process != m.at("PROCESS"))
  {
    bad += " PROCESS(card " + m.at("PROCESS") + ", JEWEL " + j.process + ")";
  }
  static const char *isoName[] = {"XX", "PP", "PN", "NP", "NN"};
  if (j.iso < 0 || j.iso > 4 || m.at("ISOCHANNEL") != isoName[j.iso])
  {
    bad += " ISOCHANNEL(card " + m.at("ISOCHANNEL") + ", JEWEL code " + std::to_string(j.iso) + ")";
  }
  if (!bad.empty())
  {
    throw std::runtime_error("PHJewel: JEWEL read the canonical card differently from PHJewel:" + bad);
  }
  os << "PHJewel: JEWEL " << j.process << " sqrt(s)=" << j.sqrts << " A=" << j.mass << " Z=" << j.nproton
     << " iso=" << m.at("ISOCHANNEL") << " pThat>=" << j.ptmin << " medium=" << c.medium
     << " recoils=" << (j.keeprecoils ? "on" : "off") << " seed=" << c.seed << " core=" << core.library() << "\n";
  return r;
}
