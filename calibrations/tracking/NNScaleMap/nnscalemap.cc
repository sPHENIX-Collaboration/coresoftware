#include "nnscalemap.h"

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrack_v4.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

namespace
{
// lower bracketing index i such that grid[i] <= x <= grid[i+1] (x already clamped to [grid.front(), grid.back()]),
// plus the fractional distance to the next grid point.
size_t bracketIndex(const std::vector<double>& grid, double x, double& frac)
{
  auto it = std::upper_bound(grid.begin(), grid.end(), x);
  size_t i = (it == grid.begin()) ? 0 : std::distance(grid.begin(), it) - 1;
  i = std::min(i, grid.size() - 2);
  const double lo = grid[i];
  const double hi = grid[i + 1];
  frac = (hi > lo) ? (x - lo) / (hi - lo) : 0.0;
  return i;
}

std::vector<double> uniqueSorted(std::vector<double> v)
{
  std::sort(v.begin(), v.end());
  std::vector<double> out;
  for (double x : v)
  {
    if (out.empty() || std::abs(x - out.back()) > 1e-6 * std::max(1.0, std::abs(x)))
    {
      out.push_back(x);
    }
  }
  return out;
}

// returns -1 if x is not found on grid within tolerance
long findIndex(const std::vector<double>& grid, double x)
{
  auto it = std::lower_bound(grid.begin(), grid.end(), x);
  if (it != grid.begin() && (it == grid.end() || std::abs(*it - x) > std::abs(*std::prev(it) - x)))
  {
    --it;
  }
  if (it == grid.end())
  {
    return -1;
  }
  const long idx = std::distance(grid.begin(), it);
  if (std::abs(grid[idx] - x) > 1e-6 * std::max(1.0, std::abs(x)))
  {
    return -1;
  }
  return idx;
}

std::string trim(const std::string& s)
{
  const size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos)
  {
    return "";
  }
  const size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

std::vector<std::string> splitCsvLine(const std::string& line)
{
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string tok;
  while (std::getline(ss, tok, ','))
  {
    out.push_back(trim(tok));
  }
  return out;
}
}  // namespace

NNScaleMap::NNScaleMap(const std::string& name)
  : SubsysReco(name)
{
}

int NNScaleMap::Init(PHCompositeNode* /*topNode*/)
{
  if (!m_useCDB && m_lookupFile.empty())
  {
    std::cout << PHWHERE << " no kappa lookup source configured -- call setKappaLookupFile(path)"
               << " and/or setUseCDB(true)." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int NNScaleMap::InitRun(PHCompositeNode* topNode)
{
  const std::string path = ResolveLookupFile();
  if (path.empty())
  {
    std::cout << PHWHERE << " could not resolve a kappa lookup file for this run"
               << (m_useCDB ? (" (CDB domain '" + m_cdbDomainName + "' and no local fallback set)") : "")
               << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (path != m_loadedFile)
  {
    if (!LoadLookupTable(path))
    {
      std::cout << PHWHERE << " failed to load kappa lookup table from '" << path << "'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    m_loadedFile = path;
  }
  else if (Verbosity() > 0)
  {
    std::cout << "NNScaleMap: '" << path << "' already loaded, reusing it for this run" << std::endl;
  }

  const int ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int NNScaleMap::process_event(PHCompositeNode* topNode)
{
  const int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  m_outputTrackMap->Reset();

  for (const auto& entry : *m_inputTrackMap)
  {
    const SvtxTrack* track = entry.second;
    if (!track)
    {
      continue;
    }
    ++m_nTracksSeen;

    const float px = track->get_px();
    const float py = track->get_py();
    const float pz = track->get_pz();
    const float pt = std::sqrt(px * px + py * py);

    SvtxTrack_v4 corrected(*track);

    if (pt <= 0.f)
    {
      ++m_nTracksSkippedZeroPt;
    }
    else
    {
      const float eta = std::asinh(pz / pt);
      const float phi = std::atan2(py, px);
      const float kappa = GetKappa(track->get_charge(), pt, eta, phi);

      corrected.set_px(kappa * px);
      corrected.set_py(kappa * py);
      corrected.set_pz(kappa * pz);

      if (m_scaleCovariance)
      {
        // momentum components live at covariance indices 3,4,5 (x,y,z,px,py,pz);
        // scaling px,py,pz by kappa scales those rows/columns of the covariance
        // by kappa for each momentum index involved (kappa^2 if both are).
        for (int i = 0; i < 6; ++i)
        {
          for (int j = i; j < 6; ++j)
          {
            const bool iMom = (i >= 3);
            const bool jMom = (j >= 3);
            if (!iMom && !jMom)
            {
              continue;
            }
            const float scale = (iMom && jMom) ? kappa * kappa : kappa;
            corrected.set_error(i, j, scale * track->get_error(i, j));
          }
        }
      }

      ++m_nTracksCorrected;
    }

    m_outputTrackMap->insertWithKey(&corrected, track->get_id());
  }

  ++m_nEvents;
  return Fun4AllReturnCodes::EVENT_OK;
}

int NNScaleMap::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "NNScaleMap::End - processed " << m_nEvents << " events, "
             << m_nTracksSeen << " tracks seen, " << m_nTracksCorrected << " corrected, "
             << m_nTracksSkippedZeroPt << " skipped (pT<=0)" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void NNScaleMap::Print(const std::string& what) const
{
  std::cout << "NNScaleMap::Print - " << what
             << ", useCDB=" << m_useCDB
             << ", cdbDomain=" << m_cdbDomainName
             << ", lookupFile=" << m_lookupFile
             << ", loaded=" << m_loadedFile
             << ", input=" << m_inputTrackMapName
             << ", output=" << m_outputTrackMapName
             << ", grid=" << m_ptEdges.size() << "(pT) x" << m_etaEdges.size()
             << "(eta) x" << m_phiEdges.size() << "(phi)"
             << std::endl;
}

std::string NNScaleMap::ResolveLookupFile() const
{
  if (m_useCDB)
  {
    // CDBInterface::getUrl resolves the payload for this run's CDB_GLOBALTAG/TIMESTAMP
    // (both must already be set on recoConsts by the macro, same as any other CDB
    // payload). If the domain has no entry (and no "<domain>_default" entry) for this
    // run, it falls back to returning m_lookupFile verbatim -- so a local override
    // still works even with CDB enabled.
    return CDBInterface::instance()->getUrl(m_cdbDomainName, m_lookupFile);
  }
  return m_lookupFile;
}

int NNScaleMap::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator iter_dst(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  auto* existing = findNode::getClass<SvtxTrackMap>(topNode, m_outputTrackMapName);
  if (existing)
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << " node " << m_outputTrackMapName << " already exists, reusing it." << std::endl;
    }
    m_outputTrackMap = existing;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_outputTrackMap = new SvtxTrackMap_v2();
  auto* node = new PHIODataNode<PHObject>(m_outputTrackMap, m_outputTrackMapName, "PHObject");
  svtxNode->addNode(node);

  if (Verbosity() > 0)
  {
    std::cout << "NNScaleMap: added " << m_outputTrackMapName << " node under DST/SVTX" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int NNScaleMap::GetNodes(PHCompositeNode* topNode)
{
  m_inputTrackMap = findNode::getClass<SvtxTrackMap>(topNode, m_inputTrackMapName);
  if (!m_inputTrackMap)
  {
    std::cout << PHWHERE << " missing required node " << m_inputTrackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (!m_outputTrackMap)
  {
    m_outputTrackMap = findNode::getClass<SvtxTrackMap>(topNode, m_outputTrackMapName);
  }
  if (!m_outputTrackMap)
  {
    std::cout << PHWHERE << " missing output node " << m_outputTrackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool NNScaleMap::LoadLookupTable(const std::string& path)
{
  std::ifstream fin(path);
  if (!fin.is_open())
  {
    return false;
  }

  // The CSV may start with '#'-prefixed convention/provenance comment lines,
  // then a header row naming its columns (e.g.
  // "q,pT,curv,eta,phi,kappa,kappa_curv,eps,delta" for a curvature-space
  // training run, or plain "q,pT,eta,phi,kappa"), then data rows. Only
  // q,pT,eta,phi,kappa are needed for deployment; look them up by name
  // instead of assuming a fixed column count/order so a future column
  // reshuffle or addition doesn't silently misread the grid.
  std::map<std::string, int> colIndex;
  std::string line;
  bool haveHeader = false;
  while (!haveHeader && std::getline(fin, line))
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    const std::vector<std::string> cols = splitCsvLine(line);
    for (size_t i = 0; i < cols.size(); ++i)
    {
      colIndex[cols[i]] = static_cast<int>(i);
    }
    haveHeader = true;
  }

  if (!haveHeader)
  {
    std::cout << PHWHERE << " kappa lookup file '" << path << "' has no header row" << std::endl;
    return false;
  }

  const char* required[] = {"q", "pT", "eta", "phi", "kappa"};
  for (const char* name : required)
  {
    if (colIndex.find(name) == colIndex.end())
    {
      std::cout << PHWHERE << " kappa lookup file '" << path
                 << "' header is missing required column '" << name << "'" << std::endl;
      return false;
    }
  }
  const int iQ = colIndex["q"];
  const int iPt = colIndex["pT"];
  const int iEta = colIndex["eta"];
  const int iPhi = colIndex["phi"];
  const int iKappa = colIndex["kappa"];
  const int nColsNeeded = 1 + std::max({iQ, iPt, iEta, iPhi, iKappa});

  std::vector<double> qs;
  std::vector<double> pts;
  std::vector<double> etas;
  std::vector<double> phis;
  std::vector<float> kappas;

  while (std::getline(fin, line))
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }

    const std::vector<std::string> cols = splitCsvLine(line);
    if (static_cast<int>(cols.size()) < nColsNeeded)
    {
      continue;
    }

    try
    {
      const double q = std::stod(cols[iQ]);
      const double pt = std::stod(cols[iPt]);
      const double eta = std::stod(cols[iEta]);
      const double phi = std::stod(cols[iPhi]);
      const double kappa = std::stod(cols[iKappa]);
      qs.push_back(q);
      pts.push_back(pt);
      etas.push_back(eta);
      phis.push_back(phi);
      kappas.push_back(static_cast<float>(kappa));
    }
    catch (const std::exception&)
    {
      continue;  // malformed data line
    }
  }

  if (qs.empty())
  {
    std::cout << PHWHERE << " kappa lookup file '" << path << "' contained no data rows" << std::endl;
    return false;
  }

  m_ptEdges = uniqueSorted(pts);
  m_etaEdges = uniqueSorted(etas);
  m_phiEdges = uniqueSorted(phis);

  const size_t nPt = m_ptEdges.size();
  const size_t nEta = m_etaEdges.size();
  const size_t nPhi = m_phiEdges.size();

  if (nPt < 2 || nEta < 2 || nPhi < 2)
  {
    std::cout << PHWHERE << " kappa lookup grid too small to interpolate (nPt=" << nPt
               << " nEta=" << nEta << " nPhi=" << nPhi << ")" << std::endl;
    return false;
  }

  m_kappaGrid[0].assign(nPt * nEta * nPhi, 1.0f);
  m_kappaGrid[1].assign(nPt * nEta * nPhi, 1.0f);

  size_t nBad = 0;
  for (size_t r = 0; r < qs.size(); ++r)
  {
    const int qIdx = (qs[r] > 0) ? 1 : 0;
    const long ip = findIndex(m_ptEdges, pts[r]);
    const long ie = findIndex(m_etaEdges, etas[r]);
    const long iph = findIndex(m_phiEdges, phis[r]);
    if (ip < 0 || ie < 0 || iph < 0)
    {
      ++nBad;
      continue;
    }
    const size_t idx = (static_cast<size_t>(ip) * nEta + static_cast<size_t>(ie)) * nPhi + static_cast<size_t>(iph);
    m_kappaGrid[qIdx][idx] = kappas[r];
  }

  if (nBad > 0)
  {
    std::cout << PHWHERE << " warning: " << nBad << " rows in '" << path
               << "' could not be placed on the (pT,eta,phi) grid" << std::endl;
  }

  std::cout << "NNScaleMap: loaded kappa lookup grid " << nPt << "(pT) x " << nEta
             << "(eta) x " << nPhi << "(phi), " << qs.size() << " rows from " << path << std::endl;

  return true;
}

float NNScaleMap::GetKappa(int charge, float pt, float eta, float phi) const
{
  if (m_ptEdges.empty())
  {
    return 1.0f;
  }

  const int qIdx = (charge > 0) ? 1 : 0;

  const double ptC = std::clamp<double>(pt, m_ptEdges.front(), m_ptEdges.back());
  const double etaC = std::clamp<double>(eta, m_etaEdges.front(), m_etaEdges.back());
  const double phiC = std::clamp<double>(phi, m_phiEdges.front(), m_phiEdges.back());

  double fpt = 0.0;
  double feta = 0.0;
  double fphi = 0.0;
  const size_t ipt = bracketIndex(m_ptEdges, ptC, fpt);
  const size_t ieta = bracketIndex(m_etaEdges, etaC, feta);
  const size_t iphi = bracketIndex(m_phiEdges, phiC, fphi);

  const size_t nEta = m_etaEdges.size();
  const size_t nPhi = m_phiEdges.size();
  const std::vector<float>& grid = m_kappaGrid[qIdx];

  auto at = [&](size_t ip, size_t ie, size_t iph) -> float
  {
    return grid[(ip * nEta + ie) * nPhi + iph];
  };

  const float c000 = at(ipt, ieta, iphi);
  const float c001 = at(ipt, ieta, iphi + 1);
  const float c010 = at(ipt, ieta + 1, iphi);
  const float c011 = at(ipt, ieta + 1, iphi + 1);
  const float c100 = at(ipt + 1, ieta, iphi);
  const float c101 = at(ipt + 1, ieta, iphi + 1);
  const float c110 = at(ipt + 1, ieta + 1, iphi);
  const float c111 = at(ipt + 1, ieta + 1, iphi + 1);

  const float c00 = static_cast<float>(c000 * (1.0 - fphi) + c001 * fphi);
  const float c01 = static_cast<float>(c010 * (1.0 - fphi) + c011 * fphi);
  const float c10 = static_cast<float>(c100 * (1.0 - fphi) + c101 * fphi);
  const float c11 = static_cast<float>(c110 * (1.0 - fphi) + c111 * fphi);

  const float c0 = static_cast<float>(c00 * (1.0 - feta) + c01 * feta);
  const float c1 = static_cast<float>(c10 * (1.0 - feta) + c11 * feta);

  return static_cast<float>(c0 * (1.0 - fpt) + c1 * fpt);
}
