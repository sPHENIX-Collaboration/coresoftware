#include "InttBadChannelMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <phool/phool.h>

#include <Rtypes.h> // For Int_t, Long64_t, etc

#include <filesystem>  // for exists

InttBadChannelMap::~InttBadChannelMap()
{
  delete m_bad_channel_set;
}

int InttBadChannelMap::LoadFromFile(std::string const& filename)
{
  if (filename.empty())
  {
    std::cerr
      << PHWHERE << "\n"
      << "\tArgument 'filename' is empty string\n"
      << std::flush;
    return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cerr
      << PHWHERE << "\n"
      << "\tFile '" << filename << "' does not exist\n"
      << std::flush;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

int InttBadChannelMap::LoadFromCDB(std::string const& name)
{
  if (name.empty())
  {
    std::cerr
      << PHWHERE << "\n"
      << "\tArgument 'name' is empty string\n"
      << std::flush;
    return 1;
  }

  std::string database = CDBInterface::instance()->getUrl(name);
  CDBTTree cdbttree(database);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

int InttBadChannelMap::v_LoadFromCDBTTree(CDBTTree& cdbttree)
{
  delete m_bad_channel_set;
  m_bad_channel_set = new Set_t;

  m_bad_channel_set->clear();
  Long64_t N = cdbttree.GetSingleIntValue("size");
  for (Long64_t n = 0; n < N; ++n)
  {
    m_bad_channel_set->insert((struct InttMap::Offline_s){
        .layer = cdbttree.GetIntValue(n, "layer"),
        .ladder_phi = cdbttree.GetIntValue(n, "ladder_phi"),
        .ladder_z = cdbttree.GetIntValue(n, "ladder_z"),
        .strip_z = cdbttree.GetIntValue(n, "strip_z"),
        .strip_phi = cdbttree.GetIntValue(n, "strip_phi"),
    });
  }

  return 0;
}

void InttBadChannelMap::identify(std::ostream& out) const
{
  out
    << PHWHERE << "\n"
    << "\tsize: " << size() << "\n"
    << std::flush;
}

std::size_t InttBadChannelMap::size() const
{
  return m_bad_channel_set ? m_bad_channel_set->size() : 0;
}

bool InttBadChannelMap::IsBad(InttMap::Online_s const& /*unused*/) const
{
  std::cerr
    << PHWHERE << "\n"
    << "\tUnimplemented overload\n"
    << std::flush;
  return false;
}

bool InttBadChannelMap::IsBad(InttMap::Offline_s const& ofl) const
{
  if (!m_bad_channel_set)
  {
    return false;
  }
  return m_bad_channel_set->find(ofl) != m_bad_channel_set->end();
}

bool InttBadChannelMap::IsBad(InttMap::RawData_s const& /*unused*/) const
{
  std::cerr
    << PHWHERE << "\n"
    << "\tUnimplemented overload\n"
    << std::flush;
  return false;
}

