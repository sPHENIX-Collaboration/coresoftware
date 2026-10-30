#include "JewelCard.h"

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <sstream>

#include <sys/stat.h>
#include <unistd.h>

namespace JewelCard
{
  namespace
  {
    enum class Type
    {
      Int,
      Real,
      Logical,
      Choice,
      Path
    };

    enum class Role
    {
      Physics,  // passed to JEWEL
      Owned,    // set by PHJewel; a value in the card is ignored (noted)
      Table,    // input table file; must exist unless generation is allowed
      Medium    // the medium card (main card only)
    };

    struct Key
    {
      const char *name;
      Type type;
      Role role;
      bool required;        // no default: must be in the card
      const char *def;      // JEWEL's built-in default ("" if required/owned)
      const char *choices;  // space-separated allowed values (Choice / Int sets)
      const char *doc;
    };

    // ---- the main card: every key JEWEL 2.6.0 init() reads -------------------
    const std::vector<Key> &mainKeys()
    {
      static const std::vector<Key> k = {
          // what the sample is -- no defaults
          {"PROCESS", Type::Choice, Role::Physics, true, "", "PPJJ PPYJ PPYQ PPYG PPZJ PPZQ PPZG PPWJ PPWQ PPWG PPDY", "hard process"},
          {"SQRTS", Type::Real, Role::Physics, true, "", "", "sqrt(s_NN) [GeV]"},
          {"PTMIN", Type::Real, Role::Physics, true, "", "", "pT-hat minimum [GeV], >= 3"},
          {"PTMAX", Type::Real, Role::Physics, true, "", "", "pT-hat maximum [GeV], -1 = none"},
          {"ETAMAX", Type::Real, Role::Physics, true, "", "", "max |eta| of ME partons and of the medium"},
          {"PDFSET", Type::Int, Role::Physics, true, "", "", "LHAPDF ID of the proton PDF (13100 = CT14nlo)"},
          {"MASS", Type::Int, Role::Physics, true, "", "", "nucleus A (1 for p+p)"},
          {"NPROTON", Type::Int, Role::Physics, true, "", "", "nucleus Z (1 for p+p)"},
          {"ISOCHANNEL", Type::Choice, Role::Physics, true, "", "XX PP PN NP NN", "XX = mix by Z/A, else one nucleon pair"},
          {"WEIGHTED", Type::Logical, Role::Physics, true, "", "", "pT-hat weighted generation"},
          {"KEEPRECOILS", Type::Logical, Role::Physics, true, "", "", "keep medium recoil partons"},
          {"SPLITINTFILE", Type::Path, Role::Table, true, "", "", "splitting-integral table"},
          {"PDFFILE", Type::Path, Role::Table, true, "", "", "medium-scattering PDF table (not the proton PDF)"},
          {"XSECFILE", Type::Path, Role::Table, true, "", "", "elastic cross-section table"},
          {"MEDIUMPARAMS", Type::Path, Role::Medium, false, "", "", "medium card (required unless medium = vac)"},
          // generation details -- JEWEL defaults, written out explicitly
          {"WEXPO", Type::Real, Role::Physics, false, "5.", "", "weight exponent (required if WEIGHTED T)"},
          {"NF", Type::Int, Role::Physics, false, "3", "3 4 5", "active flavours; must equal the medium NF"},
          {"LAMBDAQCD", Type::Real, Role::Physics, false, "0.4", "", "Lambda_QCD [GeV]"},
          {"LAMBDAQCDPS", Type::Real, Role::Physics, false, "0.4", "", "Lambda_QCD of the PYTHIA ISR [GeV]"},
          {"Q0", Type::Real, Role::Physics, false, "1.5", "", "shower cut-off [GeV]"},
          {"CHANNEL", Type::Choice, Role::Physics, false, "MUON", "MUON ELEC", "Z/W decay channel"},
          {"PDFALPHAS", Type::Real, Role::Physics, false, "0.118", "", "alpha_s(MZ) fallback"},
          {"XSECSCALEFAC", Type::Real, Role::Physics, false, "1.", "", "scale factor on the elastic cross sections"},
          {"MPI", Type::Logical, Role::Physics, false, "F", "", "PYTHIA multiple interactions"},
          {"MPIFSR", Type::Logical, Role::Physics, false, "T", "", "FSR off MPI partons"},
          {"ANGORD", Type::Int, Role::Physics, false, "1", "0 1 2 3 4 5", "angular-ordering mode"},
          {"AOFIRSTSPLIT", Type::Logical, Role::Physics, false, "F", "", "angular ordering for the first splitting"},
          {"SCATRECOIL", Type::Logical, Role::Physics, false, "F", "", "let recoils rescatter"},
          {"HADRO", Type::Logical, Role::Physics, false, "T", "", "hadronise (must be T in PHJewel)"},
          {"HADROTYPE", Type::Int, Role::Physics, false, "0", "0 1", "string assignment for hadronisation"},
          {"PI0DEC", Type::Logical, Role::Physics, false, "T", "", "decay pi0 in the generator"},
          {"ISRSCAT", Type::Logical, Role::Physics, false, "T", "", "ISR partons scatter in the medium"},
          {"SHORTHEPMC", Type::Logical, Role::Physics, false, "T", "", "final state only (F = full history)"},
          {"COMPRESS", Type::Logical, Role::Physics, false, "T", "", "compress the event record"},
          {"WRITESCATCEN", Type::Logical, Role::Physics, false, "F", "", "write thermal momenta (status 3); required T with recoils"},
          {"WRITEDUMMIES", Type::Logical, Role::Physics, false, "F", "", "dummy status-1 pi0s; refused in PHJewel"},
          {"DOSUBTRACTION", Type::Logical, Role::Physics, false, "F", "", "4MomSub inside JEWEL; refused in PHJewel"},
          {"RSUB", Type::Real, Role::Physics, false, "1.", "", "subtraction radius (unused without DOSUBTRACTION)"},
          {"PUTONSHELL", Type::Int, Role::Physics, false, "0", "0 1 2", "subtraction mass treatment (unused)"},
          {"RECSOFTCUT", Type::Real, Role::Physics, false, "0.", "", "recoil soft cut [GeV]"},
          {"RECHARDCUT", Type::Real, Role::Physics, false, "5.", "", "recoil hard cut [GeV]"},
          {"KINMODE", Type::Int, Role::Physics, false, "1", "0 1 2 3", "scattering kinematics mode"},
          {"RECMODE", Type::Int, Role::Physics, false, "0", "0 1 2 3", "recoil mode (2 switches recoils off)"},
          {"COHLENGTHFAC", Type::Real, Role::Physics, false, "1.", "", "coherence-length factor"},
          {"PCOHREJ", Type::Real, Role::Physics, false, "1.", "", "coherence rejection probability"},
          {"COHSCAT", Type::Logical, Role::Physics, false, "T", "", "coherent scattering"},
          // owned by PHJewel
          {"NEVENT", Type::Int, Role::Owned, false, "", "", "Fun4All decides (se->run(n))"},
          {"NJOB", Type::Int, Role::Owned, false, "", "", "seed comes from PHRandomSeed() or set_seed()"},
          {"LOGFILE", Type::Path, Role::Owned, false, "", "", "PHJewel::set_log_file"},
          {"HEPMCFILE", Type::Path, Role::Owned, false, "", "", "events go to the node tree"},
          {"EVINFOFILE", Type::Path, Role::Owned, false, "", "", "PHJewel::set_evinfo_file"},
          {"READRAN", Type::Logical, Role::Owned, false, "", "", "never read in.ran"},
      };
      return k;
    }

    // ---- medium cards ---------------------------------------------------------
    std::vector<Key> mediumKeys(const std::string &medium)
    {
      if (medium == "simple")
      {
        return {
            {"TAUI", Type::Real, Role::Physics, true, "", "", "initial proper time [fm/c]"},
            {"TI", Type::Real, Role::Physics, true, "", "", "initial temperature [GeV]"},
            {"TC", Type::Real, Role::Physics, true, "", "", "critical temperature [GeV]"},
            {"CENTRMIN", Type::Real, Role::Physics, true, "", "", "centrality class lower edge [%]"},
            {"CENTRMAX", Type::Real, Role::Physics, true, "", "", "centrality class upper edge [%]"},
            {"SIGMANN", Type::Real, Role::Physics, true, "", "", "inelastic NN cross section [fm^2] (4.2 = 42 mb)"},
            {"NF", Type::Int, Role::Physics, true, "", "3 4 5", "flavours in the medium; must equal main NF"},
            {"MDFACTOR", Type::Real, Role::Physics, true, "", "", "Debye mass factor"},
            {"MDSCALEFAC", Type::Real, Role::Physics, true, "", "", "Debye mass scale factor"},
            {"WOODSSAXON", Type::Logical, Role::Physics, false, "T", "", "Woods-Saxon (F = hard sphere)"},
            {"N0", Type::Real, Role::Physics, false, "0.17", "", "nuclear density [fm^-3]"},
            {"D", Type::Real, Role::Physics, false, "0.54", "", "Woods-Saxon diffuseness [fm]"},
        };
      }
      if (medium == "brick")
      {
        return {
            {"TBRICK", Type::Real, Role::Physics, true, "", "", "brick temperature [GeV]"},
            {"LBRICK", Type::Real, Role::Physics, true, "", "", "brick length [fm]"},
            {"NF", Type::Int, Role::Physics, true, "", "3 4 5", "flavours in the medium; must equal main NF"},
            {"MDFACTOR", Type::Real, Role::Physics, true, "", "", "Debye mass factor"},
            {"MDSCALEFAC", Type::Real, Role::Physics, true, "", "", "Debye mass scale factor"},
        };
      }
      return {};
    }

    // ---- helpers --------------------------------------------------------------
    std::string upper(std::string s)
    {
      for (char &c : s)
      {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      }
      return s;
    }

    std::vector<std::string> words(const std::string &s)
    {
      std::vector<std::string> w;
      std::istringstream is(s);
      std::string x;
      while (is >> x)
      {
        w.push_back(x);
      }
      return w;
    }

    // Fortran list-directed integer: optional sign, digits, nothing else.
    bool parseInt(const std::string &s, long &v)
    {
      if (s.empty())
      {
        return false;
      }
      size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
      if (i >= s.size())
      {
        return false;
      }
      for (size_t j = i; j < s.size(); ++j)
      {
        if (!std::isdigit(static_cast<unsigned char>(s[j])))
        {
          return false;
        }
      }
      v = std::strtol(s.c_str(), nullptr, 10);
      return true;
    }

    // Fortran real: [sign] digits [. digits] [(e|d) [sign] digits]; "200." and
    // "1.d0" are fine, "1,5" and "20GeV" are not.
    bool parseReal(const std::string &s, double &v)
    {
      std::string t = s;
      size_t i = 0;
      size_t n = t.size();
      if (i < n && (t[i] == '+' || t[i] == '-'))
      {
        ++i;
      }
      size_t d0 = i;
      while (i < n && std::isdigit(static_cast<unsigned char>(t[i])))
      {
        ++i;
      }
      size_t nd = i - d0;
      if (i < n && t[i] == '.')
      {
        ++i;
        size_t f0 = i;
        while (i < n && std::isdigit(static_cast<unsigned char>(t[i])))
        {
          ++i;
        }
        nd += i - f0;
      }
      if (nd == 0)
      {
        return false;
      }
      if (i < n && (t[i] == 'e' || t[i] == 'E' || t[i] == 'd' || t[i] == 'D'))
      {
        t[i] = 'e';
        ++i;
        if (i < n && (t[i] == '+' || t[i] == '-'))
        {
          ++i;
        }
        size_t e0 = i;
        while (i < n && std::isdigit(static_cast<unsigned char>(t[i])))
        {
          ++i;
        }
        if (i == e0)
        {
          return false;
        }
      }
      if (i != n)
      {
        return false;
      }
      v = std::strtod(t.c_str(), nullptr);
      return true;
    }

    bool parseLogical(const std::string &s, bool &v)
    {
      const std::string u = upper(s);
      if (u == "T" || u == ".T." || u == "TRUE" || u == ".TRUE.")
      {
        v = true;
        return true;
      }
      if (u == "F" || u == ".F." || u == "FALSE" || u == ".FALSE.")
      {
        v = false;
        return true;
      }
      return false;
    }

    bool fileExists(const std::string &p)
    {
      struct stat st{};
      return stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode);
    }

    std::string dirName(const std::string &p)
    {
      const size_t s = p.find_last_of('/');
      if (s == std::string::npos)
      {
        return ".";
      }
      if (s == 0)
      {
        return "/";
      }
      return p.substr(0, s);
    }

    std::string fmtReal(double v)
    {
      // same text as printf "%.10g"
      std::ostringstream os;
      os << std::setprecision(10) << v;
      std::string s = os.str();
      // keep it a Fortran real literal: "200" -> "200."
      if (s.find_first_of(".eEnN") == std::string::npos)
      {
        s += ".";
      }
      return s;
    }

    uint64_t fnv1a(const std::string &s)
    {
      uint64_t h = 1469598103934665603ULL;
      for (unsigned char c : s)
      {
        h ^= c;
        h *= 1099511628211ULL;
      }
      return h;
    }

    // ---- one parsed card ------------------------------------------------------
    struct Entry
    {
      std::string raw;    // value as written
      std::string canon;  // canonical value (what JEWEL is given)
      int line = 0;
    };

    struct Parsed
    {
      std::string path;
      std::map<std::string, Entry> entries;
    };

    // Read one card with the grammar both JEWEL and PHJewel agree on.
    //   maxLen : line length standalone JEWEL reads without truncation
    //            (main card CHARACTER*80, medium card CHARACTER*100)
    bool readCard(const std::string &path, const std::vector<Key> &keys, size_t maxLen,
                  Parsed &out, std::vector<std::string> &errors)
    {
      out = Parsed();
      out.path = path;
      std::ifstream in(path.c_str());
      if (!in)
      {
        errors.push_back(path + ": cannot open");
        return false;
      }
      std::map<std::string, const Key *> byName;
      for (const Key &k : keys)
      {
        byName[k.name] = &k;
      }

      const size_t before = errors.size();
      std::string line;
      int ln = 0;
      while (std::getline(in, line))
      {
        ++ln;
        const std::string where = path + ":" + std::to_string(ln) + ": ";
        if (ln > 1000)
        {
          errors.push_back(where + "JEWEL reads only the first 1000 lines");
          break;
        }
        if (!line.empty() && line[line.size() - 1] == '\r')
        {
          errors.push_back(where + "DOS line ending (carriage return); JEWEL would read it into the value");
          continue;
        }
        if (line.find('\t') != std::string::npos)
        {
          errors.push_back(where + "tab character; JEWEL splits key and value at the first SPACE only");
          continue;
        }
        if (line.find_first_not_of(' ') == std::string::npos)
        {
          continue;  // blank
        }
        if (line[0] == '#')
        {
          continue;  // comment: column 1 only
        }
        if (line.size() > maxLen)
        {
          errors.push_back(where + "line longer than " + std::to_string(maxLen) +
                           " characters; standalone JEWEL silently truncates it (use a shorter or relative path)");
          continue;
        }
        if (line[0] == ' ')
        {
          errors.push_back(where +
                           "leading space; JEWEL would read an empty key and ignore the line "
                           "(comments need '#' in column 1)");
          continue;
        }
        const size_t sp = line.find(' ');
        const std::string name = line.substr(0, sp);
        const std::string rest = (sp == std::string::npos) ? std::string() : line.substr(sp + 1);
        std::map<std::string, const Key *>::const_iterator it = byName.find(name);
        if (it == byName.end())
        {
          std::string hint;
          if (byName.contains(upper(name)))
          {
            hint = " (keys are case sensitive: " + upper(name) + ")";
          }
          std::string msg = where;
          msg += "unknown key '";
          msg += name;
          msg += "'";
          msg += hint;
          msg += "; JEWEL would ignore it and keep its default";
          errors.push_back(msg);
          continue;
        }
        const Key &k = *it->second;
        if (out.entries.contains(name))
        {
          std::string msg = where;
          msg += "duplicate key ";
          msg += name;
          msg += " (first on line ";
          msg += std::to_string(out.entries[name].line);
          msg += "); JEWEL would use the last one";
          errors.push_back(msg);
          continue;
        }
        Entry e;
        e.line = ln;
        if (k.type == Type::Path)
        {
          // JEWEL reads paths with '(a)': everything after the first space,
          // leading blanks and any trailing comment included.
          if (rest.empty() || rest[0] == ' ')
          {
            errors.push_back(where + name + ": exactly one space between key and path");
            continue;
          }
          const std::vector<std::string> w = words(rest);
          if (w.size() != 1)
          {
            std::string msg = where;
            msg += name;
            msg += ": path followed by more text ('";
            msg += rest;
            msg += "'); JEWEL would take all of it as the file name";
            errors.push_back(msg);
            continue;
          }
          e.raw = w[0];
          e.canon = w[0];  // made absolute by the caller
        }
        else
        {
          // list-directed read: first token is the value; a trailing '#...'
          // comment is harmless to JEWEL and allowed, anything else is not
          const std::vector<std::string> w = words(rest);
          if (w.empty())
          {
            errors.push_back(where + name + ": no value");
            continue;
          }
          if (w.size() > 1 && w[1][0] != '#')
          {
            std::string msg = where;
            msg += name;
            msg += ": unexpected text after the value: '";
            msg += rest;
            msg += "'";
            errors.push_back(msg);
            continue;
          }
          e.raw = w[0];
          bool ok = true;
          if (k.type == Type::Int)
          {
            long v = 0;
            ok = parseInt(w[0], v);
            if (!ok)
            {
              errors.push_back(where + name + " = '" + w[0] + "' is not an integer" +
                               (w[0].find('.') != std::string::npos ? " (no decimal point: JEWEL rejects it and keeps its default)" : ""));
            }
            else
            {
              e.canon = std::to_string(v);
            }
          }
          else if (k.type == Type::Real)
          {
            double v = 0;
            ok = parseReal(w[0], v);
            if (!ok)
            {
              errors.push_back(where + name + " = '" + w[0] + "' is not a number");
            }
            else
            {
              e.canon = fmtReal(v);
            }
          }
          else if (k.type == Type::Logical)
          {
            bool v = false;
            ok = parseLogical(w[0], v);
            if (!ok)
            {
              errors.push_back(where + name + " = '" + w[0] + "' is not a logical (T or F)");
            }
            else
            {
              e.canon = v ? "T" : "F";
            }
          }
          else if (k.type == Type::Choice)
          {
            e.canon = w[0];
          }
          if (!ok)
          {
            continue;
          }
          if (k.choices[0] != '\0')
          {
            const std::vector<std::string> allowed = words(k.choices);
            bool found = false;
            for (const std::string &a : allowed)
            {
              if (a == e.canon)
              {
                found = true;
              }
            }
            if (!found)
            {
              errors.push_back(where + name + " = '" + w[0] + "' not one of: " + k.choices);
              continue;
            }
          }
        }
        out.entries[name] = e;
      }
      return errors.size() == before;
    }
  }  // namespace

  std::string absPath(const std::string &p, const std::string &base)
  {
    std::string r;
    if (!p.empty() && p[0] == '/')
    {
      r = p;
    }
    else
    {
      std::string b = base;
      if (b.empty() || b[0] != '/')
      {
        char cwd[4096];
        if (!getcwd(cwd, sizeof(cwd)))
        {
          return p;
        }
        b = (b.empty() || b == ".") ? std::string(cwd) : std::string(cwd) + "/" + b;
      }
      r = b + "/" + p;
    }
    // collapse "/./" and "a/../" so the canonical card reads cleanly
    std::vector<std::string> parts;
    std::string seg;
    std::istringstream is(r);
    while (std::getline(is, seg, '/'))
    {
      if (seg.empty() || seg == ".")
      {
        continue;
      }
      if (seg == ".." && !parts.empty() && parts.back() != "..")
      {
        parts.pop_back();
      }
      else
      {
        parts.push_back(seg);
      }
    }
    std::string out;
    for (const std::string &x : parts)
    {
      out += "/" + x;
    }
    return out.empty() ? "/" : out;
  }

  double R(const std::map<std::string, std::string> &m, const char *k) { return std::strtod(m.at(k).c_str(), nullptr); }
  long I(const std::map<std::string, std::string> &m, const char *k) { return std::strtol(m.at(k).c_str(), nullptr, 10); }
  bool L(const std::map<std::string, std::string> &m, const char *k) { return m.at(k) == "T"; }

  Result validate(const std::string &cardPath, const Options &opt)
  {
    Result r;
    if (opt.medium != "simple" && opt.medium != "brick" && opt.medium != "vac")
    {
      r.errors.push_back("medium '" + opt.medium + "' unknown: simple | brick | vac");
      return r;
    }
    const std::string card = absPath(cardPath, ".");
    const std::string base = dirName(card);

    Parsed pm;
    readCard(card, mainKeys(), 80, pm, r.errors);

    // owned keys: note and drop
    for (const Key &k : mainKeys())
    {
      if (k.role == Role::Owned && pm.entries.contains(k.name))
      {
        r.notes.push_back(std::string(k.name) + " " + pm.entries[k.name].raw + " in the card is ignored: " + k.doc);
        pm.entries.erase(k.name);
      }
    }
    // required keys, then defaults made explicit
    for (const Key &k : mainKeys())
    {
      if (k.role == Role::Owned)
      {
        continue;
      }
      if (pm.entries.contains(k.name))
      {
        r.main[k.name] = pm.entries[k.name].canon;
        r.explicitMain.insert(k.name);
      }
      else if (k.required)
      {
        r.errors.push_back(card + ": required key " + k.name + " missing (" + k.doc + "); " +
                           (k.role == Role::Table ? "JEWEL would look for its default file name in the working directory"
                                                  : "JEWEL's built-in default is for LHC Pb+Pb"));
      }
      else if (k.role != Role::Medium)
      {
        r.main[k.name] = k.def;
      }
    }
    if (!pm.path.empty() && r.main.empty() && !r.errors.empty())
    {
      return r;  // unreadable card
    }

    // ---- physics checks on the main card ----
    // Each check runs when the keys it needs are there, so one pass reports
    // every problem in the card.
    std::vector<std::string> &E = r.errors;
    const std::map<std::string, std::string> &m = r.main;
    auto has = [&m](std::initializer_list<const char *> ks)
    {
      for (const char *k : ks)
      {
        if (!m.contains(k))
        {
          return false;
        }
      }
      return true;
    };
    if (has({"SQRTS"}) && (R(m, "SQRTS") <= 0))
    {
      E.emplace_back("SQRTS must be > 0");
    }
    if (has({"PTMIN"}) && (R(m, "PTMIN") < 3.0))
    {
      E.emplace_back("PTMIN < 3 GeV: JEWEL silently raises it to 3");
    }
    if (has({"PTMAX", "PTMIN"}) && (R(m, "PTMAX") != -1.0 && R(m, "PTMAX") <= R(m, "PTMIN")))
    {
      E.emplace_back("PTMAX must be -1 (no maximum) or larger than PTMIN");
    }
    if (has({"ETAMAX"}) && (R(m, "ETAMAX") <= 0))
    {
      E.emplace_back("ETAMAX must be > 0");
    }
    if (has({"MASS"}) && (I(m, "MASS") < 1))
    {
      E.emplace_back("MASS must be >= 1");
    }
    if (has({"NPROTON", "MASS"}) && (I(m, "NPROTON") < 0 || I(m, "NPROTON") > I(m, "MASS")))
    {
      E.emplace_back("NPROTON must be in 0..MASS");
    }
    if (has({"PDFSET"}) && (I(m, "PDFSET") <= 0))
    {
      E.emplace_back("PDFSET must be a positive LHAPDF ID");
    }
    if (R(m, "LAMBDAQCD") <= 0 || R(m, "LAMBDAQCDPS") <= 0 || R(m, "Q0") <= 0)
    {
      E.emplace_back("LAMBDAQCD, LAMBDAQCDPS and Q0 must be > 0");
    }
    if (R(m, "XSECSCALEFAC") <= 0)
    {
      E.emplace_back("XSECSCALEFAC must be > 0");
    }
    if (has({"WEIGHTED"}) && (L(m, "WEIGHTED") && !pm.entries.contains("WEXPO")))
    {
      E.emplace_back("WEIGHTED T needs an explicit WEXPO (JEWEL's 5. is the usual choice; write it down)");
    }
    if (has({"WEIGHTED"}) && (L(m, "WEIGHTED") && R(m, "WEXPO") <= 0))
    {
      E.emplace_back("WEXPO must be > 0");
    }
    if (!L(m, "HADRO"))
    {
      E.emplace_back("HADRO F leaves partons as final-state particles; Geant cannot transport them");
    }
    if (L(m, "WRITEDUMMIES"))
    {
      E.emplace_back(
          "WRITEDUMMIES T adds one status-1 pi0 placeholder (E = 1e-6 GeV) per thermal parton, hundreds "
          "per event, which Geant and any truth-particle count treat as real particles; recoil subtraction "
          "needs only WRITESCATCEN T");
    }
    if (L(m, "DOSUBTRACTION"))
    {
      E.emplace_back(
          "DOSUBTRACTION T rewrites the final state inside JEWEL; a detector simulation needs the "
          "unsubtracted event (subtract at analysis level from the status-3 thermal momenta)");
    }
    if (has({"KEEPRECOILS"}) && (L(m, "KEEPRECOILS") && !L(m, "WRITESCATCEN")))
    {
      E.emplace_back(
          "KEEPRECOILS T needs WRITESCATCEN T: without the thermal momenta (status 3) the recoil "
          "contribution can never be subtracted");
    }
    if (has({"KEEPRECOILS"}) && (I(m, "RECMODE") == 2 && L(m, "KEEPRECOILS")))
    {
      E.emplace_back("RECMODE 2 makes JEWEL switch KEEPRECOILS off; the card says T");
    }
    const std::string iso = m.contains("ISOCHANNEL") ? m.at("ISOCHANNEL") : "";
    if (has({"MASS", "NPROTON", "ISOCHANNEL"}) && I(m, "MASS") == 1 && iso == "XX")
    {
      r.notes.push_back("MASS 1 with ISOCHANNEL XX: every event is the p+p channel (Z/A = " +
                        m.at("NPROTON") + "/1)");
    }
    if (has({"MASS", "ISOCHANNEL"}) && I(m, "MASS") > 1 && iso != "XX")
    {
      r.notes.push_back("ISOCHANNEL " + iso + " with MASS " + m.at("MASS") +
                        ": one nucleon pair only, no isospin mixing");
    }
    if (L(m, "MPI"))
    {
      r.notes.emplace_back("MPI T: PYTHIA 6 underlying event; JEWEL's MPI partons are not quenched");
    }
    if (pm.entries.contains("RSUB") || pm.entries.contains("PUTONSHELL"))
    {
      r.notes.emplace_back("RSUB / PUTONSHELL have no effect without DOSUBTRACTION");
    }

    // ---- tables ----
    for (const char *t : {"SPLITINTFILE", "PDFFILE", "XSECFILE"})
    {
      if (!m.contains(t))
      {
        continue;
      }
      const std::string p = absPath(m.at(t), base);
      r.main[t] = p;
      if (!fileExists(p))
      {
        if (opt.allowTableGeneration)
        {
          r.notes.push_back(std::string(t) + " " + p +
                            " does not exist: JEWEL will integrate and write it "
                            "(slow; never in parallel jobs sharing the path)");
        }
        else
        {
          E.push_back(std::string(t) + " " + p +
                      " does not exist; JEWEL would spend a long time integrating it "
                      "and write it into a shared area (PHJewel::allow_table_generation to permit)");
        }
      }
    }
    r.notes.emplace_back(
        "tables are only as good as the parameters they were built with (Q0, LAMBDAQCD, NF, "
        "the medium Debye-mass range and the energy reach); JEWEL counts extrapolations, "
        "PHJewel prints them at End");

    // ---- medium card ----
    std::string medIn;
    if (pm.entries.contains("MEDIUMPARAMS"))
    {
      medIn = absPath(pm.entries["MEDIUMPARAMS"].raw, base);
    }
    if (opt.medium == "vac")
    {
      if (!medIn.empty())
      {
        r.notes.emplace_back("MEDIUMPARAMS is ignored by the vacuum medium");
      }
    }
    else
    {
      if (medIn.empty())
      {
        E.push_back("MEDIUMPARAMS missing: without it the " + opt.medium +
                    " medium silently runs with built-in defaults");
      }
      else if (!fileExists(medIn))
      {
        E.push_back("MEDIUMPARAMS " + medIn +
                    " does not exist: JEWEL would print a notice to its log "
                    "and run with built-in defaults");
      }
      else
      {
        Parsed pd;
        const std::vector<Key> mk = mediumKeys(opt.medium);
        readCard(medIn, mk, 100, pd, E);
        for (const Key &k : mk)
        {
          if (pd.entries.contains(k.name))
          {
            r.medium[k.name] = pd.entries[k.name].canon;
            r.explicitMedium.insert(k.name);
          }
          else if (k.required)
          {
            E.push_back(medIn + ": required key " + k.name + " missing (" + k.doc + ")");
          }
          else
          {
            r.medium[k.name] = k.def;
          }
        }
        r.mediumPath = medIn;
      }
    }
    // medium checks run whenever the medium card itself was read completely
    const std::map<std::string, std::string> &d = r.medium;
    const bool medOK = opt.medium != "vac" && d.size() == mediumKeys(opt.medium).size();
    if (medOK)
    {
      if (d.contains("NF") && I(d, "NF") != I(m, "NF"))
      {
        E.push_back("NF " + d.at("NF") + " in the medium card but " + m.at("NF") + " in the main card");
      }
      if (R(d, "MDFACTOR") <= 0 || R(d, "MDSCALEFAC") <= 0)
      {
        E.emplace_back("MDFACTOR and MDSCALEFAC must be > 0");
      }
    }
    if (medOK && opt.medium == "simple")
    {
      if (m.contains("MASS") && I(m, "MASS") < 2)
      {
        E.emplace_back("the simple medium needs a nucleus: MASS >= 2");
      }
      if (R(d, "TAUI") <= 0)
      {
        E.emplace_back("TAUI must be > 0");
      }
      if (!(R(d, "TI") > R(d, "TC") && R(d, "TC") > 0))
      {
        E.emplace_back("need TI > TC > 0 (GeV)");
      }
      if (R(d, "TI") > 2.0)
      {
        E.push_back("TI " + d.at("TI") + " looks like MeV; JEWEL wants GeV");
      }
      if (!(R(d, "CENTRMIN") >= 0 && R(d, "CENTRMIN") < R(d, "CENTRMAX") && R(d, "CENTRMAX") <= 100))
      {
        E.emplace_back("need 0 <= CENTRMIN < CENTRMAX <= 100");
      }
      if (R(d, "SIGMANN") <= 0)
      {
        E.emplace_back("SIGMANN must be > 0");
      }
      if (R(d, "SIGMANN") > 15.0)
      {
        E.push_back("SIGMANN " + d.at("SIGMANN") + " looks like mb; JEWEL wants fm^2 (42 mb = 4.2)");
      }
      if (R(d, "N0") <= 0 || R(d, "D") <= 0)
      {
        E.emplace_back("N0 and D must be > 0");
      }
    }
    if (medOK && opt.medium == "brick")
    {
      if (R(d, "TBRICK") <= 0 || R(d, "TBRICK") > 2.0)
      {
        E.emplace_back("TBRICK must be in (0, 2] GeV");
      }
      if (R(d, "LBRICK") <= 0)
      {
        E.emplace_back("LBRICK must be > 0 (fm)");
      }
    }
    if (!E.empty())
    {
      return r;
    }

    // ---- canonical cards ----
    // Keys the user wrote are passed to JEWEL; the others are listed as
    // comments with JEWEL's built-in value but NOT passed: several built-in
    // defaults are single-precision literals (lqcd = 0.4 is 0.4000000059),
    // so writing "LAMBDAQCD 0.4" would change the physics in the 8th digit
    // and break event-for-event identity with the standalone binary.
    std::ostringstream med;
    if (opt.medium != "vac")
    {
      for (const Key &k : mediumKeys(opt.medium))
      {
        if (!r.explicitMedium.contains(k.name))
        {
          med << "#default ";
        }
        med << k.name << " " << d.at(k.name) << "\n";
      }
    }
    r.mediumCanon = med.str();

    std::ostringstream mc;
    for (const Key &k : mainKeys())
    {
      if (k.role == Role::Owned || k.role == Role::Medium)
      {
        continue;
      }
      if (!r.explicitMain.contains(k.name))
      {
        mc << "#default ";
      }
      mc << k.name << " " << r.main[k.name] << "\n";
    }
    r.mainCanon = mc.str();  // MEDIUMPARAMS added by the writer (path of the copy)

    // 16 hex digits, as printf "%016llx"
    std::ostringstream hx;
    hx << std::hex << std::setw(16) << std::setfill('0')
       << fnv1a("medium=" + opt.medium + "\n" + r.mainCanon + "--\n" + r.mediumCanon);
    r.hash = hx.str();
    r.ok = true;
    return r;
  }

  std::string writeCanonical(const Result &r, const std::string &stem, const std::string &source,
                             const std::string &medium, std::string &err)
  {
    std::string mainPath = absPath(stem + ".card.dat", ".");
    const std::string medPath = absPath(stem + ".medium.dat", ".");
    if (medium != "vac")
    {
      std::ofstream o(medPath.c_str());
      if (!o)
      {
        err = "cannot write " + medPath;
        return "";
      }
      o << "# canonical JEWEL medium card (" << medium << ") written by PHJewel from\n# "
        << r.mediumPath << "\n# card hash " << r.hash << "\n"
        << r.mediumCanon;
    }
    std::ofstream o(mainPath.c_str());
    if (!o)
    {
      err = "cannot write " + mainPath;
      return "";
    }
    o << "# canonical JEWEL card written by PHJewel from\n# " << source
      << "\n# medium " << medium << ", card hash " << r.hash
      << "\n# every parameter is explicit; NEVENT/NJOB/LOGFILE/HEPMCFILE/EVINFOFILE are set by PHJewel\n"
      << r.mainCanon;
    if (medium != "vac")
    {
      o << "MEDIUMPARAMS " << medPath << "\n";
    }
    return mainPath;
  }
}  // namespace JewelCard
