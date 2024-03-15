#include "InttSurveyMapv1.h"

#include <cdbobjects/CDBTTree.h>

#include <boost/format.hpp>

InttSurveyMapv1::~InttSurveyMapv1()
{
  delete m_absolute_transforms;
  delete m_relative_transforms;
}

void InttSurveyMapv1::identify(
    std::ostream& out) const
{
  out << "InttSurveyMapv1"
      << "\n"
      << "\tsize: " << size() << std::endl;
}

std::size_t InttSurveyMapv1::size() const
{
  if (!m_absolute_transforms)
  {
    return 0;
  }
  return m_absolute_transforms->size();
}

int InttSurveyMapv1::v_LoadFromCDBTTree(
    CDBTTree& cdbttree)
{
  Eigen::Affine3d aff;
  InttMap::Offline_s ofl;

  delete m_absolute_transforms;
  delete m_relative_transforms;

  m_absolute_transforms = new map_t;
  m_relative_transforms = new map_t;

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
      std::string boost_formatted = boost::str(boost::format("m_abs_%01d_%01d") %  (i/4) % (i%4));
      aff.matrix()(i / 4, i % 4) = cdbttree.GetDoubleValue(n, boost_formatted);
    }
    m_absolute_transforms->insert({ofl, aff});

    for (int i = 0; i < 16; ++i)
    {
      std::string boost_formatted = boost::str(boost::format("m_rel_%01d_%01d") %  (i/4) % (i%4));
      aff.matrix()(i / 4, i % 4) = cdbttree.GetDoubleValue(n,boost_formatted );
    }
    m_relative_transforms->insert({ofl, aff});
  }

  return 0;
}

InttSurveyMap::val_t const* InttSurveyMapv1::GetAbsoluteTransform(
    key_t const& k) const
{
  if (!m_absolute_transforms)
  {
    return nullptr;
  }
  map_t::const_iterator itr = m_absolute_transforms->find(k);
  return itr != m_absolute_transforms->end() ? &itr->second : nullptr;
}

InttSurveyMap::val_t const* InttSurveyMapv1::GetRelativeTransform(
    key_t const& k) const
{
  if (!m_relative_transforms)
  {
    return nullptr;
  }
  map_t::const_iterator itr = m_relative_transforms->find(k);
  return itr != m_relative_transforms->end() ? &itr->second : nullptr;
}
