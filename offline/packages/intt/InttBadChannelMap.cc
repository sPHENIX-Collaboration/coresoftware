#include "InttBadChannelMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <phool/phool.h>

#include <Rtypes.h> // For Int_t, Long64_t, etc

#include <filesystem>  // for exists
#include <limits>

int InttBadChannelMap::Load(std::string const& name)
{
  if (name.empty())
  {
    std::cout
      << PHWHERE << "\n"
      << "\tArgument 'name' is empty string\n"
      << std::flush;
    return 1;
  }

  std::string filename = name.find(".root") != std::string::npos ? name : CDBInterface::instance()->getUrl(name);

  if (filename.empty())
  {
    std::cout
      << PHWHERE << "\n"
      << "\tCalibration '" << name << "' not found in CDB\n"
      << std::flush;
    return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cout
      << PHWHERE << "\n"
      << "\tFile '" << filename << "' not does not exist\n"
      << std::flush;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();
  return v_LoadFromCDBTTree(cdbttree);
}

int InttBadChannelMap::v_LoadFromCDBTTree(CDBTTree& cdbttree)
{
  m_offline_set.clear();
  m_rawdata_set.clear();

  m_offline_loaded = true;
  m_rawdata_loaded = true;
  m_size = cdbttree.GetSingleIntValue("size");

  if (m_size == 0)
  {
    return 0;
  }

  // Check if the CDBTTree has branches corresponding to the offline convention
  m_offline_loaded = m_offline_loaded && (cdbttree.GetIntValue(0, "layer") != std::numeric_limits<int>::min());
  m_offline_loaded = m_offline_loaded && (cdbttree.GetIntValue(0, "ladder_phi") != std::numeric_limits<int>::min());
  m_offline_loaded = m_offline_loaded && (cdbttree.GetIntValue(0, "ladder_z") != std::numeric_limits<int>::min());
  m_offline_loaded = m_offline_loaded && (cdbttree.GetIntValue(0, "strip_z") != std::numeric_limits<int>::min());
  m_offline_loaded = m_offline_loaded && (cdbttree.GetIntValue(0, "strip_phi") != std::numeric_limits<int>::min());

  // Check if the CDBTTree has branches corresponding to the rawdata convention
  m_rawdata_loaded = m_rawdata_loaded && (cdbttree.GetIntValue(0, "felix_server") != std::numeric_limits<int>::min());
  m_rawdata_loaded = m_rawdata_loaded && (cdbttree.GetIntValue(0, "felix_channel") != std::numeric_limits<int>::min());
  m_rawdata_loaded = m_rawdata_loaded && (cdbttree.GetIntValue(0, "chip") != std::numeric_limits<int>::min());
  m_rawdata_loaded = m_rawdata_loaded && (cdbttree.GetIntValue(0, "channel") != std::numeric_limits<int>::min());

  if (!m_offline_loaded && !m_rawdata_loaded)
  {
    std::cout
      << PHWHERE << "\n"
      << "\tCDBTTree does not have expected calibrations\n"
      << "\tAvailable calibrations:\n"
      << std::flush;
    cdbttree.Print();

    return 1;
  }

  for (Long64_t n = 0; n < m_size; ++n)
  {
    if (m_offline_loaded)
    {
      m_offline_set.insert((struct InttNameSpace::Offline_s){
        .layer = cdbttree.GetIntValue(n, "layer"),
        .ladder_phi = cdbttree.GetIntValue(n, "ladder_phi"),
        .ladder_z = cdbttree.GetIntValue(n, "ladder_z"),
        .strip_x = cdbttree.GetIntValue(n, "strip_phi"),
        .strip_y = cdbttree.GetIntValue(n, "strip_z"),
      });
    }

    if (m_rawdata_loaded)
    {
      m_rawdata_set.insert((struct InttNameSpace::RawData_s){
        .felix_server = cdbttree.GetIntValue(n, "felix_server"),
        .felix_channel = cdbttree.GetIntValue(n, "felix_channel"),
        .chip = cdbttree.GetIntValue(n, "chip"),
        .channel = cdbttree.GetIntValue(n, "channel"),
      });
    }
  }

  return 0;
}

void InttBadChannelMap::identify(std::ostream& out) const
{
  out
    << PHWHERE << "\n"
    << "\tsize: " << m_size << "\n"
    << std::flush;
}

void InttBadChannelMap::Print(std::ostream& out) const
{
  if (!m_offline_loaded && !m_rawdata_loaded)
  {
    out
      << PHWHERE << "\n"
      << "\tNo channels loaded\n"
      << std::flush;
    return;
  }

  out
    << PHWHERE << "\n";
  if (m_rawdata_loaded)
  {
    out << "\tmasked channels (rawdata convention): " << m_rawdata_set.size() << "\n";
    for (auto const& raw : m_rawdata_set)
    {
      out << "\t" << raw.felix_server << " " << raw.felix_channel << " " << raw.chip << " " << raw.channel << "\n";
    }
    out << std::flush;
  }

  if (m_offline_loaded)
  {
    out << "\tmasked channels (offline convention): " << m_offline_set.size() << "\n";
    for (auto const& ofl : m_offline_set)
    {
      out << "\t" << ofl.layer << " " << ofl.ladder_phi << " " << ofl.ladder_z << " " << ofl.strip_y << " " << ofl.strip_x << "\n";
    }
    out << std::flush;
  }
}

bool InttBadChannelMap::IsBad(InttNameSpace::Online_s const& /*unused*/) const
{
  std::cout
    << PHWHERE << "\n"
    << "\tUnimplemented overload\n"
    << std::flush;
  return false;
}

bool InttBadChannelMap::IsBad(InttNameSpace::Offline_s const& ofl) const
{
  if(m_offline_loaded)
  {
    return m_offline_set.find(ofl) != m_offline_set.end();
  }

  std::cout
    << PHWHERE << "\n"
    << "\tMap was not loaded with offline indexing convention\n"
    << "\t(layer, ladder_phi/z, strip_z/phi)\n"
    << "\tUse a different overload\n"
    << std::flush;
  return false;
}

bool InttBadChannelMap::IsBad(InttNameSpace::RawData_s const& raw) const
{
  if(m_rawdata_loaded)
  {
    return m_rawdata_set.find(raw) != m_rawdata_set.end();
  }

  std::cout
    << PHWHERE << "\n"
    << "\tMap was not loaded with rawdata indexing convention\n"
    << "\t(felix_server/channel, chip, channel)\n"
    << "\tUse a different overload\n"
    << std::flush;
  return false;
}

