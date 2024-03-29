#include "InttFeeMapv1.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <filesystem>  // for exists
#include <utility>     // for pair

InttFeeMapv1::~InttFeeMapv1()
{
  delete m_onl_to_raw;
  delete m_raw_to_onl;
}
void InttFeeMapv1::identify(
  std::ostream& out) const
{
  out << "InttFeeMapv1" << std::endl;
}

int InttFeeMapv1::Convert(
  InttMap::Online_s& onl,
  InttMap::Offline_s const& ofl) const
{
  if(!InttMap::IsValid(ofl))
  {
	  return 1;
  }
  int num_ldrs = ofl.layer < 5 ? 12 : 16;

  onl.lyr = ofl.layer - 3;
  onl.ldr = (7 * num_ldrs / 4 - ofl.ladder_phi) % num_ldrs;
  // onl.ldr = (7 * num_ldrs / 4 - ofl.ladder_phi + (ofl.layer % 2 ? num_ldrs - 1 : 0)) % num_ldrs;
  onl.arm = ofl.ladder_z / 2;
  switch(ofl.ladder_z)
  {
  case 1:
	  onl.chp = ofl.strip_z + 13 * (ofl.strip_phi < 128);
	  break;
  case 0:
	  onl.chp = ofl.strip_z + 13 * (ofl.strip_phi < 128) + 5;
	  break;
  case 2:
	  onl.chp = 12 - ofl.strip_z + 13 * (127 < ofl.strip_phi);
	  break;
  case 3:
	  onl.chp = 4 - ofl.strip_z + 13 * (127 < ofl.strip_phi);
	  break;
  default:
	  break;
  }
  onl.chn = ofl.strip_phi < 128 ? ofl.strip_phi : 255 - ofl.strip_phi;
  return 0;
}

int InttFeeMapv1::Convert(
  InttMap::Offline_s& ofl,
  InttMap::Online_s const& onl) const
{
  if(!InttMap::IsValid(onl))
  {
	  return 1;
  }
  int num_ldrs = onl.lyr < 2 ? 12 : 16;

  ofl.layer = onl.lyr + 3;
  ofl.ladder_phi = (7 * num_ldrs / 4 - onl.ldr) % num_ldrs;
  // ofl.ladder_phi = (7 * num_ldrs / 4 - onl.ldr + (onl.lyr % 2 ? 0 : num_ldrs - 1)) % num_ldrs;
  ofl.ladder_z = 2 * onl.arm + (onl.chp % 13 < 5);
  switch(ofl.ladder_z)
  {
  case 1:
	  ofl.strip_z = onl.chp % 13;
	  break;
  case 0:
	  ofl.strip_z = onl.chp % 13 - 5;
	  break;
  case 2:
	  ofl.strip_z = 12 - (onl.chp % 13);
	  break;
  case 3:
	  ofl.strip_z = 4 - (onl.chp % 13);
	  break;
  default:
	  break;
  }
  ofl.strip_phi = (onl.arm == (onl.chp / 13)) ? 255 - onl.chn : onl.chn;
  return 0;
}

int InttFeeMapv1::Convert(
  InttMap::RawData_s& raw,
  InttMap::Online_s const& onl) const
{
  if(!m_onl_to_raw)
  {
    return 1;
  }

  if(m_onl_to_raw->find(onl) != m_onl_to_raw->end())
  {
    raw = m_onl_to_raw->at(onl);
	raw.chp = onl.chp;
	raw.chn = onl.chn;

	return 0;
  }

  return 1;
}

int InttFeeMapv1::Convert(
  InttMap::Online_s& onl,
  InttMap::RawData_s const& raw) const
{
  if(!m_raw_to_onl)
  {
    return 1;
  }

  if(m_raw_to_onl->find(raw) != m_raw_to_onl->end())
  {
    onl = m_raw_to_onl->at(raw);
	onl.chp = raw.chp;
	onl.chn = raw.chn;

	return 0;
  }

  return 1;
}

int InttFeeMapv1::Convert(
  InttMap::RawData_s& raw,
  InttMap::Offline_s const& ofl) const
{
  InttMap::Online_s onl;

  if(Convert(onl, ofl))
  {
    return 1;
  }
  if(Convert(raw, onl))
  {
    return 1;
  }
  return 0;
}

int InttFeeMapv1::Convert(
  InttMap::Offline_s& ofl,
  InttMap::RawData_s const& raw) const
{
  InttMap::Online_s onl;

  if(Convert(onl, raw))
  {
    return 1;
  }
  if(Convert(ofl, onl))
  {
    return 1;
  }
  return 0;
}

int InttFeeMapv1::v_LoadFromCDBTTree(
    CDBTTree& cdbttree)
{
  InttMap::Online_s onl;
  onl.chp = InttMap::Wildcard;
  onl.chn = InttMap::Wildcard;

  InttMap::RawData_s raw;
  raw.chp = InttMap::Wildcard;
  raw.chn = InttMap::Wildcard;

  delete m_onl_to_raw;
  delete m_raw_to_onl;

  m_onl_to_raw = new std::map<InttMap::Online_s, InttMap::RawData_s, InttMap::OnlineWildcardComparator>;
  m_raw_to_onl = new std::map<InttMap::RawData_s, InttMap::Online_s, InttMap::RawDataWildcardComparator>;

  Int_t N = cdbttree.GetSingleIntValue("size");
  for(Int_t n = 0; n < N; ++n)
  {
	  onl.lyr = cdbttree.GetIntValue(n, "lyr");
	  onl.ldr = cdbttree.GetIntValue(n, "ldr");
	  onl.arm = cdbttree.GetIntValue(n, "arm");

	  raw.pid = cdbttree.GetIntValue(n, "pid");
	  raw.fee = cdbttree.GetIntValue(n, "fee");

	  m_onl_to_raw->insert({onl, raw});
	  m_raw_to_onl->insert({raw, onl});
  }

  return 0;
}
