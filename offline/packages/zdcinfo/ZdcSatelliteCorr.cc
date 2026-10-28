#include "ZdcSatelliteCorr.h"

#include "Zdcinfo.h"
#include "Zdcinfov2.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>

#include <mbd/MbdOut.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <cmath>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>

ZdcSatelliteCorr::ZdcSatelliteCorr(const std::string &name)
  : SubsysReco(name)
{
}

ZdcSatelliteCorr::~ZdcSatelliteCorr()
{
  delete m_cdb;
}

int ZdcSatelliteCorr::InitRun(PHCompositeNode *topNode)
{
  delete m_cdb;
  m_cdb = nullptr;
  m_have = false;
  m_gain = {1., 1.};
  m_fm5 = {1., 1.};
  m_fp5 = {1., 1.};
  m_status = {1, 1};
  m_nt0bin = 0;

  int run = 0;
  EventHeader *eh = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (eh)
  {
    run = eh->get_RunNumber();
  }
  if (run <= 0)
  {
    run = static_cast<int>(recoConsts::instance()->get_uint64Flag("TIMESTAMP"));
  }

  std::string url;
  if (!m_calibFile.empty())
  {
    if (!std::filesystem::exists(m_calibFile))
    {
      std::cout << Name() << ": local calibration file " << m_calibFile << " not found, abort run" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    url = m_calibFile;
    std::cout << Name() << ": TEST, constants from the local file " << url << ", not from the CDB" << std::endl;
  }
  else
  {
    url = CDBInterface::instance()->getUrl(m_calibName);
  }
  if (url.empty())
  {
    std::cout << Name() << ": no " << m_calibName << " payload for this run, ZDC energies are not corrected" << std::endl;
    return createNodes(topNode);
  }
  m_cdb = new CDBTTree(url);
  m_cdb->LoadCalibrations();
  m_have = true;

  if (!m_calibFile.empty())
  {
    const int rfirst = m_cdb->GetSingleIntValue("run_first");
    const int rlast = m_cdb->GetSingleIntValue("run_last");
    if (run < rfirst || run > rlast)
    {
      std::cout << Name() << ": WARNING run " << run << " is outside the file's run range " << rfirst << "-" << rlast
                << ", its constants are applied anyway" << std::endl;
    }
  }

  for (int arm = 0; arm < 2; ++arm)
  {
    const float g = m_cdb->GetFloatValue(arm, "gain");
    const float fm = m_cdb->GetFloatValue(arm, "f_m5");
    const float fp = m_cdb->GetFloatValue(arm, "f_p5");
    const float st = m_cdb->GetFloatValue(arm, "status");
    m_gain[arm] = (std::isfinite(g) && g > 0) ? g : 1.;
    m_fm5[arm] = (std::isfinite(fm) && fm > 0) ? fm : 1.;
    m_fp5[arm] = (std::isfinite(fp) && fp > 0) ? fp : 1.;
    m_status[arm] = (std::isfinite(st) && st < 0.5) ? 0 : 1;
  }

  // t0-binned f
  if (m_useT0Bins)
  {
    m_nt0bin = m_cdb->GetSingleIntValue("nt0bin");
    m_t0lo = m_cdb->GetSingleFloatValue("t0bin_lo");
    m_t0w = m_cdb->GetSingleFloatValue("t0bin_w");
    if (m_nt0bin <= 0 || !(m_t0w > 0))
    {
      m_nt0bin = 0;
    }
    for (int arm = 0; arm < 2 && m_nt0bin > 0; ++arm)
    {
      m_fb[arm].assign(m_nt0bin, 1.);
      for (int k = 0; k < m_nt0bin; ++k)
      {
        const float v = m_cdb->GetFloatValue(arm, std::format("fb{:02d}", k));
        if (std::isfinite(v) && v > 0)
        {
          m_fb[arm][k] = v;
        }
      }
    }
  }

  // per-run entries: main-bunch scale
  const int nruns = m_cdb->GetSingleIntValue("nruns_t0");
  for (int i = 0; i < nruns; ++i)
  {
    if (m_cdb->GetIntValue(100 + i, "run") != run)
    {
      continue;
    }
    if (m_perRunGain)
    {
      const std::array<float, 2> gr{m_cdb->GetFloatValue(100 + i, "gain_run_s"), m_cdb->GetFloatValue(100 + i, "gain_run_n")};
      for (int arm = 0; arm < 2; ++arm)
      {
        if (std::isfinite(gr[arm]) && gr[arm] > 0 && m_status[arm] == 1)
        {
          m_gain[arm] = gr[arm];
        }
      }
    }
    break;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": " << url << " run " << run << ", t0 bins " << m_nt0bin
              << ", S gain " << m_gain[0] << " f(-5) " << m_fm5[0] << " f(+5) " << m_fp5[0] << " status " << m_status[0]
              << ", N gain " << m_gain[1] << " f(-5) " << m_fm5[1] << " f(+5) " << m_fp5[1] << " status " << m_status[1] << std::endl;
  }
  return createNodes(topNode);
}

int ZdcSatelliteCorr::createNodes(PHCompositeNode *topNode)
{
  if (m_overwrite || findNode::getClass<Zdcinfo>(topNode, "ZdcinfoCorr"))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << Name() << ": no DST node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator dstIter(dstNode);
  PHCompositeNode *zdcNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "ZDC"));
  if (!zdcNode)
  {
    zdcNode = new PHCompositeNode("ZDC");
    dstNode->addNode(zdcNode);
  }
  zdcNode->addNode(new PHIODataNode<PHObject>(new Zdcinfov2(), "ZdcinfoCorr", "PHObject"));
  return Fun4AllReturnCodes::EVENT_OK;
}

float ZdcSatelliteCorr::factor(const int arm, const float t0) const
{
  // returns G x f for this arm and (re-centered) t0, or 0 when the event is not corrected
  if (!m_have || m_status[arm] == 0 || !std::isfinite(t0) || std::fabs(t0) >= m_edge2)
  {
    return 0.;
  }
  float f = 1.;
  if (t0 < -m_edge1)
  {
    f = m_fm5[arm];
  }
  else if (t0 >= m_edge1)
  {
    f = m_fp5[arm];
  }
  if (m_nt0bin > 0)
  {
    const int k = static_cast<int>(std::floor((t0 - m_t0lo) / m_t0w));
    if (k >= 0 && k < m_nt0bin)
    {
      f = m_fb[arm][k];
    }
  }
  return m_gain[arm] * f;
}

int ZdcSatelliteCorr::process_event(PHCompositeNode *topNode)
{
  ++m_nevt;
  Zdcinfo *in = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  if (!in)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  Zdcinfo *out = m_overwrite ? in : findNode::getClass<Zdcinfo>(topNode, "ZdcinfoCorr");
  if (!out)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  MbdOut *mbdout = findNode::getClass<MbdOut>(topNode, "MbdOut");
  const float t0 = (mbdout && mbdout->get_npmt(0) > 0 && mbdout->get_npmt(1) > 0) ? mbdout->get_t0() : std::numeric_limits<float>::quiet_NaN();
  bool corrected = false;
  for (int arm = 0; arm < 2; ++arm)
  {
    const float e = in->get_zdc_energy(arm);
    const float gf = factor(arm, t0);
    if (gf > 0)
    {
      out->set_zdc_energy(arm, e / gf);
      corrected = true;
    }
    else if (!m_overwrite)
    {
      out->set_zdc_energy(arm, e);
    }
  }
  if (corrected)
  {
    ++m_ncorr;
  }
  if (!m_overwrite)
  {
    out->set_zvertex(in->get_zvertex());
    for (int arm = 0; arm < 2; ++arm)
    {
      out->set_radius(arm, in->get_radius(arm));
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int ZdcSatelliteCorr::End(PHCompositeNode * /*topNode*/)
{
  std::cout << Name() << ": corrected " << m_ncorr << " of " << m_nevt << " events" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
