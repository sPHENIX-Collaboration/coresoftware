#include "InttSurveyMap.h"

#include <cdbobjects/CDBTTree.h>

#include <boost/format.hpp>

InttSurveyMap::~InttSurveyMap()
{
  delete m_transforms;
}

int InttSurveyMap::GetStripTransform(
  key_t k,
  val_t& v) const
{
  val_t const* transform_ptr = nullptr;
  key_t ofl{
    k.layer,
    k.ladder_phi,
    k.ladder_z,
    InttMap::Wildcard,
    InttMap::Wildcard
  };

  transform_ptr = GetTransform(ofl);
  if(!transform_ptr)
  {
    return 1;
  }

  for(int i = 0; i < 16; ++i)
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

int InttSurveyMap::GetSensorTransform(
  key_t k,
  val_t& v) const
{
  val_t const* transform_ptr = nullptr;
  key_t ofl{
    k.layer,
    k.ladder_phi,
    k.ladder_z,
    InttMap::Wildcard,
    InttMap::Wildcard
  };

  transform_ptr = GetTransform(ofl);
  if(!transform_ptr)
  {
    return 1;
  }

  v = *transform_ptr;

  return 0;
}

int InttSurveyMap::GetLadderTransform(
  key_t k,
  val_t& v) const
{
  val_t const* transform_ptr = nullptr;
  k.ladder_z  = InttMap::Wildcard;
  k.strip_z   = InttMap::Wildcard;
  k.strip_phi = InttMap::Wildcard;

  transform_ptr = GetTransform(k);
  if(!transform_ptr)
  {
    return 1;
  }

  v = *transform_ptr;

  return 0;
}

InttSurveyMap::val_t const* InttSurveyMap::GetTransform(
    key_t const& k) const
{
  if (!m_transforms)
  {
    return nullptr;
  }
  map_t::const_iterator itr = m_transforms->find(k);
  return itr != m_transforms->end() ? &itr->second : nullptr;
}

int InttSurveyMap::LoadFromCDBTTree(
    CDBTTree& cdbttree)
{
  delete m_transforms;
  m_transforms = new map_t;

  Eigen::Affine3d aff;
  InttMap::Offline_s ofl;
  std::map<std::string, InttMap::field_t*> ofl_fields = {
    {"layer",      &ofl.layer},
    {"ladder_phi", &ofl.ladder_phi},
    {"ladder_z",   &ofl.ladder_z},
    {"strip_z",    &ofl.strip_z},
    {"strip_phi",  &ofl.strip_phi},
  };

  for (int n = 0, N = cdbttree.GetSingleIntValue("size"); n < N; ++n)
  {
    for(auto& p : ofl_fields)
    {
      if((*p.second = cdbttree.GetIntValue(n, p.first)) != std::numeric_limits<int>::min())
      {
        continue;
      }

      std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\tCDBTTree::GetIntValue returned std::numeric_limits<int>::min()\n"
                << "\tname, channel combination does not exist in \"cdbttree\"\n" << std::endl;
      return 1;
    }

    for (int i = 0; i < 16; ++i)
    {
      std::string boost_formatted = boost::str(boost::format("m_abs_%01d_%01d") %  (i/4) % (i%4));
      if((aff.matrix()(i / 4, i % 4) = cdbttree.GetDoubleValue(n, boost_formatted)) != std::numeric_limits<int>::min())
	  {
		  continue;
	  }

	  std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\tCDBTTree::GetIntValue returned std::numeric_limits<int>::min()\n"
                << "\tname, channel combination does not exist in \"cdbttree\"\n" << std::endl;
	  return 1;
    }

    m_transforms->insert({ofl, aff});
  }

  return 0;
}
