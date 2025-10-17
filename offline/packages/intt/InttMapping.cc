#include "InttMapping.h"
#include "InttFelixMap.h"

#include <ffarawobjects/InttRawHit.h>

#include <format>

/// Struct methods
InttNameSpace::RawData_s& InttNameSpace::RawData_s::operator++()
{
  if (++channel < 128)
  {
    return *this;
  }
  channel = 0;

  if (++chip < 26)
  {
    return *this;
  }
  chip = 0;

  if (++felix_channel < 14)
  {
    return *this;
  }
  felix_channel = 0;

  ++felix_server;
  return *this;
}

bool InttNameSpace::operator<(RawData_s const& lhs, RawData_s const& rhs)
{
  if (lhs.felix_server != rhs.felix_server)
  {
    return lhs.felix_server < rhs.felix_server;
  }
  if (lhs.felix_channel != rhs.felix_channel)
  {
    return lhs.felix_channel < rhs.felix_channel;
  }
  if (lhs.chip != rhs.chip)
  {
    return lhs.chip < rhs.chip;
  }
  if (lhs.channel != rhs.channel)
  {
    return lhs.channel < rhs.channel;
  }
  return false;
}

std::ostream& InttNameSpace::operator<<(std::ostream& out, RawData_s const& rawdata)
{
  return out << std::format(
             "RawData_s {{ .felix_server = {}, .felix_channel = {:2}, .chip = {:2}, .channel = {:3} }}",
             rawdata.felix_server, rawdata.felix_channel, rawdata.chip, rawdata.channel);
}

InttNameSpace::Online_s& InttNameSpace::Online_s::operator++()
{
  if (++chn < 128)
  {
    return *this;
  }
  chn = 0;

  if (++chp < 26)
  {
    return *this;
  }
  chp = 0;

  if (++arm < 2)
  {
    return *this;
  }
  arm = 0;

  if (++ldr < (lyr < 5 ? 12 : 16))
  {
    return *this;
  }
  ldr = 0;

  ++lyr;
  return *this;
}

bool InttNameSpace::operator<(Online_s const& lhs, Online_s const& rhs)
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

std::ostream& InttNameSpace::operator<<(std::ostream& out, Online_s const& online)
{
  return out << std::format(
             "Online_s {{ .lyr = {}, .ldr = {:2}, .arm = {:2}, .chp = {:2}, .chn = {:3} }}",
             online.lyr, online.ldr, online.arm, online.chp, online.chn);
}

InttNameSpace::Offline_s& InttNameSpace::Offline_s::operator++()
{
  if (++strip_x < 256)
  {
    return *this;
  }
  strip_x = 0;

  if (++strip_y < (ladder_z % 2 ? 5 : 8))
  {
    return *this;
  }
  strip_y = 0;

  if (++ladder_z < 4)
  {
    return *this;
  }
  ladder_z = 0;

  if (++ladder_phi < (layer < 2 ? 12 : 16))
  {
    return *this;
  }
  ladder_phi = 0;

  ++layer;
  return *this;
}

bool InttNameSpace::operator<(Offline_s const& lhs, Offline_s const& rhs)
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
  if (lhs.strip_x != rhs.strip_x)
  {
    return lhs.strip_x < rhs.strip_x;
  }
  if (lhs.strip_y != rhs.strip_y)
  {
    return lhs.strip_y < rhs.strip_y;
  }
  return false;
}

std::ostream& InttNameSpace::operator<<(std::ostream& out, Offline_s const& offline)
{
  return out << std::format(
             "Offline_s {{ .layer = {}, .ladder_phi = {:2}, .ladder_z = {:2}, .strip_x = {:2}, .strip_y = {:3} }}",
             offline.layer, offline.ladder_phi, offline.ladder_z, offline.strip_x, offline.strip_y);
}

/// Namespace-scope methods
InttNameSpace::RawData_s InttNameSpace::RawFromHit(InttRawHit* hit)
{
  return {
      .felix_server = hit->get_packetid() - 3001,
      .felix_channel = hit->get_fee(),
      .chip = (hit->get_chip_id() + 25) % 26,
      .channel = hit->get_channel_id(),
  };
}

InttNameSpace::Online_s InttNameSpace::ToOnline(Offline_s const& offline)
{
  Online_s online;
  int n_ldr = offline.layer < 5 ? 12 : 16;

  online.lyr = offline.layer - 3;
  // online.ldr = (7 * n_ldr / 4 - offline.ladder_phi + (offline.layer % 2 ? n_ldr - 1 : 0)) % n_ldr;
  online.ldr = (7 * n_ldr / 4 - offline.ladder_phi) % n_ldr;

  online.arm = offline.ladder_z / 2;
  switch (offline.ladder_z)
  {
  case 1:
    online.chp = offline.strip_y + 13 * (offline.strip_x < 128);
    break;

  case 0:
    online.chp = offline.strip_y + 13 * (offline.strip_x < 128) + 5;
    break;

  case 2:
    online.chp = 12 - offline.strip_y + 13 * !(offline.strip_x < 128);
    break;

  case 3:
    online.chp = 4 - offline.strip_y + 13 * !(offline.strip_x < 128);
    break;

  default:
    break;
  }

  online.chn = (offline.strip_x < 128) ? offline.strip_x : 255 - offline.strip_x;

  return online;
}

InttNameSpace::Offline_s InttNameSpace::ToOffline(Online_s const& online)
{
  Offline_s offline;
  int n_ldr = online.lyr < 2 ? 12 : 16;

  offline.layer = online.lyr + 3;
  // offline.ladder_phi = (7 * n_ldr / 4 - online.ldr + (online.lyr % 2 ? 0 : n_ldr - 1)) % n_ldr;
  offline.ladder_phi = (7 * n_ldr / 4 - online.ldr) % n_ldr;

  offline.ladder_z = 2 * online.arm + (online.chp % 13 < 5);
  switch (offline.ladder_z)
  {
  case 1:
    offline.strip_y = online.chp % 13;
    break;

  case 0:
    offline.strip_y = online.chp % 13 - 5;
    break;

  case 2:
    offline.strip_y = 12 - (online.chp % 13);
    break;

  case 3:
    offline.strip_y = 4 - (online.chp % 13);
    break;

  default:
    break;
  }

  offline.strip_x = (online.arm == (online.chp / 13)) ? 255 - online.chn : online.chn;

  return offline;
}

InttNameSpace::RawData_s InttNameSpace::ToRawData(Online_s const& online)
{
  RawData_s rawdata;

  InttFelix::OnlineToRawData(online, rawdata);
  rawdata.chip = online.chp;
  rawdata.channel = online.chn;

  return rawdata;
}

InttNameSpace::Online_s InttNameSpace::ToOnline(RawData_s const& rawdata)
{
  Online_s online;

  InttFelix::RawDataToOnline(rawdata, online);
  online.chp = rawdata.chip;
  online.chn = rawdata.channel;

  return online;
}
