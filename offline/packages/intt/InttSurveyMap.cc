#include "InttSurveyMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <phool/phool.h>

#include <Rtypes.h> // For Int_t, Long64_t, etc

#include <filesystem>  // for exists
#include <utility>     // for pair

#include <boost/format.hpp>

InttSurveyMap::~InttSurveyMap()
{
  delete m_absolute_transforms;
}

int InttSurveyMap::LoadFromFile(std::string const& filename)
{
  if (filename.empty())
  {
    std::cerr
	  << PHWHERE << "\n"
      << "\tArgument 'filename' is empty string"
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

int InttSurveyMap::LoadFromCDB(
    std::string const& name)
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

int InttSurveyMap::v_LoadFromCDBTTree(CDBTTree& cdbttree)
{
  Eigen::Affine3d aff;
  InttMap::Offline_s ofl;

  delete m_absolute_transforms;
  m_absolute_transforms = new map_t;

  Int_t N = cdbttree.GetSingleIntValue("size");
  for (Int_t n = 0; n < N; ++n)
  {
    ofl.layer = cdbttree.GetIntValue(n, "layer");
    ofl.ladder_phi = cdbttree.GetIntValue(n, "ladder_phi");
    ofl.ladder_z = cdbttree.GetIntValue(n, "ladder_z");
    ofl.strip_z = cdbttree.GetIntValue(n, "strip_z");
    ofl.strip_phi = cdbttree.GetIntValue(n, "strip_phi");

    for (int i = 0; i < 16; ++i)
    {
      std::string boost_formatted = boost::str(boost::format("m_abs_%01d_%01d") % (i / 4) % (i % 4));
      aff.matrix()(i / 4, i % 4) = cdbttree.GetDoubleValue(n, boost_formatted);
    }
    m_absolute_transforms->insert({ofl, aff});
  }

  return 0;
}

int InttSurveyMap::GetStripTransform(key_t const& k, val_t& v) const
{
  val_t const* transform_ptr = nullptr;
  key_t ofl{
      k.layer,
      k.ladder_phi,
      k.ladder_z,
      InttMap::Wildcard,
      InttMap::Wildcard};

  transform_ptr = GetAbsoluteTransform(ofl);
  if (!transform_ptr)
  {
    return 1;
  }

  for (int i = 0; i < 16; ++i)
  {
    v.matrix()(i / 4, i % 4) = i / 4 == i % 4;
  }

  // strip_z determines local z
  //  10, 16 are double the max range strip_z takes for the sensor type (ladder_z % 2)
  //  100, 128 are the lengths of the sensor type (ladder_z % 2) in mm
  v.matrix()(2, 3) = (2.0 * k.strip_z + 1.0) / ((k.ladder_z % 2) ? 10.0 : 16.0) - 0.5;
  v.matrix()(2, 3) *= (k.ladder_z % 2) ? 100.0 : 128.0;

  // strip_phi determines the local x
  // 512 is double the max range strip_phi takes
  // 19.968 is the sensor width in mm
  v.matrix()(0, 3) = (2.0 * k.strip_phi + 1.0) / 512.0 - 0.5;
  v.matrix()(0, 3) *= 19.968;

  v = *transform_ptr * v;

  return 0;
}

int InttSurveyMap::GetSensorTransform(key_t const& k, val_t& v) const
{
  val_t const* transform_ptr = nullptr;
  key_t ofl{
      k.layer,
      k.ladder_phi,
      k.ladder_z,
      InttMap::Wildcard,
      InttMap::Wildcard};

  transform_ptr = GetAbsoluteTransform(ofl);
  if (!transform_ptr)
  {
    return 1;
  }

  v = *transform_ptr;

  return 0;
}

int InttSurveyMap::GetLadderTransform(key_t const& k, val_t& v) const
{
  val_t const* transform_ptr = nullptr;
  key_t ofl{
      k.layer,
      k.ladder_phi,
      InttMap::Wildcard,
      InttMap::Wildcard,
      InttMap::Wildcard};

  transform_ptr = GetAbsoluteTransform(ofl);
  if (!transform_ptr)
  {
    return 1;
  }

  v = *transform_ptr;

  return 0;
}

void InttSurveyMap::identify(std::ostream& out) const
{
  out
    << PHWHERE << "\n"
    << std::flush;
}

std::size_t InttSurveyMap::size() const
{
  return m_absolute_transforms ? m_absolute_transforms->size() : 0;
}

InttSurveyMap::val_t const* InttSurveyMap::GetAbsoluteTransform(key_t const& k) const
{
  if (!m_absolute_transforms)
  {
    return nullptr;
  }
  map_t::const_iterator itr = m_absolute_transforms->find(k);
  return itr != m_absolute_transforms->end() ? &itr->second : nullptr;
}

