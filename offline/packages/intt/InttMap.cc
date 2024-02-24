#include "InttMap.h"

bool InttMap::OnlineComparator::operator()(
    struct InttMap::Online_s const& lhs,
    struct InttMap::Online_s const& rhs) const
{
  if (lhs.lyr != rhs.lyr)
  {
    return lhs.lyr < rhs.lyr;
  }
  if (lhs.ldr != rhs.ldr)
  {
    return lhs.ldr < rhs.ldr;
  }
  if (lhs.arm != rhs.arm)
  {
    return lhs.arm < rhs.arm;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool InttMap::OnlineWildcardComparator::operator()(
    struct InttMap::Online_s const& lhs,
    struct InttMap::Online_s const& rhs) const
{
  if (lhs.lyr == Wildcard || rhs.lyr == Wildcard)
  {
    return false;
  }
  if (lhs.lyr != rhs.lyr)
  {
    return lhs.lyr < rhs.lyr;
  }
  // if(lhs.lyr != Wildcard && rhs.lyr == Wildcard && lhs.lyr != rhs.lyr)return lhs.lyr < rhs.lyr;

  if (lhs.ldr == Wildcard || rhs.ldr == Wildcard)
  {
    return false;
  }
  if (lhs.ldr != rhs.ldr)
  {
    return lhs.ldr < rhs.ldr;
  }
  // if(lhs.ldr != Wildcard && rhs.ldr != Wildcard && lhs.ldr != rhs.ldr)return lhs.ldr < rhs.ldr;

  if (lhs.arm == Wildcard || rhs.arm == Wildcard)
  {
    return false;
  }
  if (lhs.arm != rhs.arm)
  {
    return lhs.arm < rhs.arm;
  }
  // if(lhs.arm != Wildcard && rhs.arm != Wildcard && lhs.arm != rhs.arm)return lhs.arm < rhs.arm;

  if (lhs.chp == Wildcard || rhs.chp == Wildcard)
  {
    return false;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  // if(lhs.chp != Wildcard && rhs.chp != Wildcard && lhs.chp != rhs.chp)return lhs.chp < rhs.chp;

  if (lhs.chn == Wildcard || rhs.chn == Wildcard)
  {
    return false;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }
  // if(lhs.chn != Wildcard && rhs.chn != Wildcard && lhs.chn != rhs.chn)return lhs.chn < rhs.chn;

  return false;
}

bool InttMap::RawDataComparator::operator()(
    struct InttMap::RawData_s const& lhs,
    struct InttMap::RawData_s const& rhs) const
{
  if (lhs.pid != rhs.pid)
  {
    return lhs.pid < rhs.pid;
  }
  if (lhs.fee != rhs.fee)
  {
    return lhs.fee < rhs.fee;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool InttMap::RawDataWildcardComparator::operator()(
    struct InttMap::RawData_s const& lhs,
    struct InttMap::RawData_s const& rhs) const
{
  if (lhs.pid == Wildcard || rhs.pid == Wildcard)
  {
    return false;
  }
  if (lhs.pid != rhs.pid)
  {
    return lhs.pid < rhs.pid;
  }
  // if(lhs.pid != Wildcard && rhs.pid != Wildcard && lhs.pid != rhs.pid)return lhs.pid < rhs.pid;

  if (lhs.fee == Wildcard || rhs.fee == Wildcard)
  {
    return false;
  }
  if (lhs.fee != rhs.fee)
  {
    return lhs.fee < rhs.fee;
  }
  // if(lhs.fee != Wildcard && rhs.fee != Wildcard && lhs.fee != rhs.fee)return lhs.fee < rhs.fee;

  if (lhs.chp == Wildcard || rhs.chp == Wildcard)
  {
    return false;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  // if(lhs.chp != Wildcard && rhs.chp != Wildcard && lhs.chp != rhs.chp)return lhs.chp < rhs.chp;

  if (lhs.chn == Wildcard || rhs.chn == Wildcard)
  {
    return false;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }
  // if(lhs.chn != Wildcard && rhs.chn != Wildcard && lhs.chn != rhs.chn)return lhs.chn < rhs.chn;

  return false;
}

bool InttMap::OfflineComparator::operator()(
    struct InttMap::Offline_s const& lhs,
    struct InttMap::Offline_s const& rhs) const
{
  if (lhs.layer != rhs.layer)
  {
    return lhs.layer < rhs.layer;
  }
  if (lhs.ladder_phi != rhs.ladder_phi)
  {
    return lhs.ladder_phi < rhs.ladder_phi;
  }
  if (lhs.ladder_z != rhs.ladder_z)
  {
    return lhs.ladder_z < rhs.ladder_z;
  }
  if (lhs.strip_z != rhs.strip_z)
  {
    return lhs.strip_z < rhs.strip_z;
  }
  if (lhs.strip_phi != rhs.strip_phi)
  {
    return lhs.strip_phi < rhs.strip_phi;
  }

  return false;
}

bool InttMap::OfflineWildcardComparator::operator()(
    struct InttMap::Offline_s const& lhs,
    struct InttMap::Offline_s const& rhs) const
{
  if (lhs.layer == Wildcard || rhs.layer == Wildcard)
  {
    return false;
  }
  if (lhs.layer != rhs.layer)
  {
    return lhs.layer < rhs.layer;
  }
  // if(lhs.layer != Wildcard && rhs.layer != Wildcard && lhs.layer != rhs.layer)return lhs.layer < rhs.layer;

  if (lhs.ladder_phi == Wildcard || rhs.ladder_phi == Wildcard)
  {
    return false;
  }
  if (lhs.ladder_phi != rhs.ladder_phi)
  {
    return lhs.ladder_phi < rhs.ladder_phi;
  }
  // if(lhs.ladder_phi != Wildcard && rhs.ladder_phi != Wildcard && lhs.ladder_phi != rhs.ladder_phi)return lhs.ladder_phi < rhs.ladder_phi;

  if (lhs.ladder_z == Wildcard || rhs.ladder_z == Wildcard)
  {
    return false;
  }
  if (lhs.ladder_z != rhs.ladder_z)
  {
    return lhs.ladder_z < rhs.ladder_z;
  }
  // if(lhs.ladder_z != Wildcard && rhs.ladder_z != Wildcard && lhs.ladder_z != rhs.ladder_z)return lhs.ladder_z < rhs.ladder_z;

  if (lhs.strip_z == Wildcard || rhs.strip_z == Wildcard)
  {
    return false;
  }
  if (lhs.strip_z != rhs.strip_z)
  {
    return lhs.strip_z < rhs.strip_z;
  }
  // if(lhs.strip_z != Wildcard && rhs.strip_z != Wildcard && lhs.strip_z != rhs.strip_z)return lhs.strip_z < rhs.strip_z;

  if (lhs.strip_phi == Wildcard || rhs.strip_phi == Wildcard)
  {
    return false;
  }
  if (lhs.strip_phi != rhs.strip_phi)
  {
    return lhs.strip_phi < rhs.strip_phi;
  }
  // if(lhs.strip_phi != Wildcard && rhs.strip_phi != Wildcard && lhs.strip_phi != rhs.strip_phi)return lhs.strip_phi < rhs.strip_phi;

  return false;
}
