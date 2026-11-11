#include "TpcQAInternal.h"

namespace TpcCalibrationQADetail
{
  bool sameAxis(const TAxis *a, const TAxis *b)
  {
    if (a->GetNbins() != b->GetNbins())
    {
      return false;
    }
    for (int i = 1; i <= a->GetNbins() + 1; ++i)
    {
      if (std::abs(a->GetBinLowEdge(i) - b->GetBinLowEdge(i)) > 1e-10)
      {
        return false;
      }
    }
    return true;
  }

  bool sameBins(const TH1 *a, const TH1 *b)
  {
    return a->GetDimension() == b->GetDimension() && sameAxis(a->GetXaxis(), b->GetXaxis()) &&
           (a->GetDimension() < 2 || sameAxis(a->GetYaxis(), b->GetYaxis())) &&
           (a->GetDimension() < 3 || sameAxis(a->GetZaxis(), b->GetZaxis()));
  }

  std::vector<std::string> wanted()
  {
    std::vector<std::string> names = {"nrecotracks1d", "pt_pos", "pt_neg",
                                      "nrecotracks_pos", "nrecotracks_neg", "tpc_only_crossing",
                                      "tpc_only_status"};
    for (const std::string side : {"south", "north"})
    {
      names.push_back("tpc_only_" + side + "_rdca_qoverpt");
      for (const std::string q : {"pos", "neg"})
      {
        names.push_back(std::string("dcaxyorigin_phi_").append(side).append("_").append(q));
        for (const std::string obs : {"rdca_phi_pt", "rdca_eta_pt", "pt_phi_eta"})
        {
          names.push_back(
              std::string("tpc_only_").append(side).append("_").append(q).append("_").append(obs));
        }
      }
    }
    return names;
  }

  Hists load(const Sample &sample, const Config &cfg, std::ostream &availability)
  {
    std::ifstream in(sample.list);
    if (!in)
    {
      throw std::runtime_error("Cannot read file list: " + sample.list);
    }
    std::set<std::string> seen;
    std::map<std::string, int> present;
    Hists sums;
    int files = 0;
    std::string path;
    while (std::getline(in, path))
    {
      if (path.empty() || path[0] == '#')
      {
        continue;
      }
      if (!seen.insert(path).second)
      {
        throw std::runtime_error("Duplicate input: " + path);
      }
      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
      {
        throw std::runtime_error("Bad ROOT input: " + path);
      }
      ++files;
      for (const auto &name : wanted())
      {
        auto *h = dynamic_cast<TH1 *>(file->Get((cfg.prefix + name).c_str()));
        if (!h)
        {
          continue;
        }
        int dim = 1;
        if (name.find("rdca_qoverpt") != std::string::npos || name.find("dcaxyorigin_phi") == 0 ||
            name == "nrecotracks_pos" || name == "nrecotracks_neg")
        {
          dim = 2;
        }
        else if (name.find("_pt_phi_eta") != std::string::npos ||
                 name.find("_rdca_phi_pt") != std::string::npos ||
                 name.find("_rdca_eta_pt") != std::string::npos)
        {
          dim = 3;
        }
        if (h->GetDimension() != dim)
        {
          throw std::runtime_error("Wrong histogram dimension: " + name);
        }
        ++present[name];
        auto it = sums.find(name);
        if (it == sums.end())
        {
          Hist copy(static_cast<TH1 *>(h->Clone(name.c_str())));
          copy->SetDirectory(nullptr);
          sums.emplace(name, std::move(copy));
        }
        else
        {
          if (!sameBins(it->second.get(), h))
          {
            throw std::runtime_error("Incompatible segment binning: " + name);
          }
          if (!it->second->Add(h))
          {
            throw std::runtime_error("Cannot merge: " + name);
          }
        }
      }
    }
    if (!files)
    {
      throw std::runtime_error("Empty input list: " + sample.list);
    }
    for (const auto &name : wanted())
    {
      std::string state = "absent";
      if (present[name] == files)
      {
        state = "complete";
      }
      else if (present[name])
      {
        state = "partial_SKIPPED";
      }
      availability << sample.label << '\t' << name << '\t' << present[name] << '\t' << files << '\t'
                   << state << '\n';
      if (present[name] != files)
      {
        sums.erase(name);
      }
    }
    return sums;
  }
}  // namespace TpcCalibrationQADetail
