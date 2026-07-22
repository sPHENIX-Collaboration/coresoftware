#include "IdealPadMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

IdealPadMap::IdealPadMap()
{
  m_radius_cm_by_layer.fill(std::numeric_limits<double>::quiet_NaN());
  m_layer_by_key.fill(-1);
}

IdealPadMap::~IdealPadMap() = default;

int IdealPadMap::layer_index(const unsigned int layer) const
{
  if (layer < FIRST_LAYER || layer > LAST_LAYER) { return -1;
}
  return static_cast<int>(layer - FIRST_LAYER);
}

int IdealPadMap::get_region(const unsigned int layer) const
{
  const int ilayer = layer_index(layer);
  if (ilayer < 0) { return -1;
}
  return ilayer / 16;
}

double IdealPadMap::wrap_phi(const double phi) const
{
  double out = phi;
  while (out <= -M_PI) { out += 2.0 * M_PI;
}
  while (out > M_PI) { out -= 2.0 * M_PI;
}
  return out;
}

int IdealPadMap::load_from_cdb(const int /*verbosity*/)
{
  m_is_loaded = false;

  for (auto& v : m_cdb_phi_by_layer) { v.clear();
}
  m_radius_cm_by_layer.fill(std::numeric_limits<double>::quiet_NaN());
  m_layer_by_key.fill(-1);

  CDBInterface* cdb = CDBInterface::instance();
  const std::string calibdir = cdb->getUrl("TPC_FEE_CHANNEL_MAP");

  if (calibdir.empty())
  {
    std::cout << "IdealPadMap::load_from_cdb - no TPC_FEE_CHANNEL_MAP found" << std::endl;
    return -1;
  }

  CDBTTree* cdbttree = new CDBTTree(calibdir);
  cdbttree->LoadCalibrations();

  std::array<std::vector<double>, N_LAYERS> radius_mm_by_layer;

  for (unsigned int fee = 0; fee < N_FEE; ++fee)
  {
    for (unsigned int ch = 0; ch < N_CH; ++ch)
    {
      const unsigned int key = 256U * fee + ch;
      const int layer = cdbttree->GetIntValue(key, "layer");
      m_layer_by_key[key] = layer;

      const int ilayer = layer_index(static_cast<unsigned int>(layer));
      if (ilayer < 0) { continue;
}

      m_cdb_phi_by_layer[ilayer].push_back(cdbttree->GetDoubleValue(key, "phi"));
      radius_mm_by_layer[ilayer].push_back(cdbttree->GetDoubleValue(key, "R"));
    }
  }

  delete cdbttree;
  cdbttree = nullptr;

  for (unsigned int ilayer = 0; ilayer < N_LAYERS; ++ilayer)
  {
    if (m_cdb_phi_by_layer[ilayer].empty() || radius_mm_by_layer[ilayer].empty())
    {
      std::cout << "IdealPadMap::load_from_cdb - missing CDB entries for layer "
                << ilayer + FIRST_LAYER << std::endl;
      return -1;
    }

   
    std::sort(m_cdb_phi_by_layer[ilayer].begin(), m_cdb_phi_by_layer[ilayer].end());

    double radius_mm = 0.0;
    for (const double r : radius_mm_by_layer[ilayer]) { radius_mm += r;
}
    radius_mm /= static_cast<double>(radius_mm_by_layer[ilayer].size());

    // CDB R is in mm.  Most TPC tracking/display code uses cm.
    m_radius_cm_by_layer[ilayer] = radius_mm / 10.0;
  }

  m_is_loaded = true;

 // if (verbosity > 0)
  {
    std::cout << "IdealPadMap::load_from_cdb - loaded " << calibdir << std::endl;
    for (unsigned int region = 0; region < N_REGIONS; ++region)
    {
      const unsigned int first_layer = FIRST_LAYER + 16U * region;
      const unsigned int pads_per_sector = get_pads_per_sector(region);
      std::cout << "  region " << region
                << " first_layer " << first_layer
                << " pads_per_sector " << pads_per_sector
                << " total_phibins " << get_total_phibins(first_layer)
                << std::endl;
      if (pads_per_sector == 0U) { continue;
}

      const unsigned int first_pad = 0U;
      const unsigned int last_pad = pads_per_sector - 1U;
      for (unsigned int side = 0; side < N_SIDES; ++side)
      {
        for (unsigned int sector = 0; sector < N_SECTORS; ++sector)
        {
          std::cout << "    side " << side
                    << " sector " << sector
                    << " first_pad " << first_pad
                    << " first_phi " << get_phi(side, sector, first_layer, first_pad)
                    << " last_pad " << last_pad
                    << " last_phi " << get_phi(side, sector, first_layer, last_pad)
                    << std::endl;
        }
      }
    }
  }

  return 0;
}

unsigned int IdealPadMap::get_pads_per_sector(const unsigned int region) const
{
  if (region >= N_REGIONS) { return 0U;
}

 
  const unsigned int layer = FIRST_LAYER + 16U * region;
  return get_pads_per_sector_for_layer(layer);
}

unsigned int IdealPadMap::get_pads_per_sector_for_layer(const unsigned int layer) const
{
  const int ilayer = layer_index(layer);
  if (ilayer < 0) { return 0U;
}
  return static_cast<unsigned int>(m_cdb_phi_by_layer[ilayer].size());
}

unsigned int IdealPadMap::get_total_phibins(const unsigned int layer) const
{
  return N_SECTORS * get_pads_per_sector_for_layer(layer);
}

double IdealPadMap::get_radius(const unsigned int layer) const
{
  const int ilayer = layer_index(layer);
  if (ilayer < 0) { return std::numeric_limits<double>::quiet_NaN();
}
  return m_radius_cm_by_layer[ilayer];
}

double IdealPadMap::get_layer_thickness(const unsigned int layer) const
{
  const int ilayer = layer_index(layer);
  if (ilayer < 0) { return std::numeric_limits<double>::quiet_NaN();
}

  if (ilayer == 0)
  {
    return m_radius_cm_by_layer[1] - m_radius_cm_by_layer[0];
  }
  if (ilayer == static_cast<int>(N_LAYERS - 1))
  {
    return m_radius_cm_by_layer[N_LAYERS - 1] - m_radius_cm_by_layer[N_LAYERS - 2];
  }

  return 0.5 * (m_radius_cm_by_layer[ilayer + 1] - m_radius_cm_by_layer[ilayer - 1]);
}

double IdealPadMap::get_cdb_local_phi(const unsigned int layer,
                                      const unsigned int local_phibin) const
{
  const int ilayer = layer_index(layer);
  if (ilayer < 0) { return std::numeric_limits<double>::quiet_NaN();
}

  const auto& phi_vec = m_cdb_phi_by_layer[ilayer];
  if (local_phibin >= phi_vec.size()) { return std::numeric_limits<double>::quiet_NaN();
}

  return phi_vec[local_phibin];
}

double IdealPadMap::get_phi(const unsigned int side,
                            const unsigned int layer,
                            const unsigned int phibin) const
{
  const unsigned int pads_per_sector = get_pads_per_sector_for_layer(layer);
  if (pads_per_sector == 0U) { return std::numeric_limits<double>::quiet_NaN();
}

  const unsigned int sector = phibin / pads_per_sector;
  const unsigned int local_phibin = phibin % pads_per_sector;

  return get_phi(side, sector, layer, local_phibin);
}

double IdealPadMap::get_phi(const unsigned int side,
                            const unsigned int sector,
                            const unsigned int layer,
                            const unsigned int local_phibin) const
{
  if (side >= N_SIDES) { return std::numeric_limits<double>::quiet_NaN();
}
  if (sector >= N_SECTORS) { return std::numeric_limits<double>::quiet_NaN();
}

  unsigned int lookup_phibin = local_phibin;
  if (side == 1U)
  {
    const unsigned int pads_per_sector = get_pads_per_sector_for_layer(layer);
    if (pads_per_sector == 0U) { return std::numeric_limits<double>::quiet_NaN();
}
    if (local_phibin >= pads_per_sector) { return std::numeric_limits<double>::quiet_NaN();
}
    lookup_phibin = pads_per_sector - 1U - local_phibin;
  }

  const double cdb_phi = get_cdb_local_phi(layer, lookup_phibin);
  if (!std::isfinite(cdb_phi)) { return std::numeric_limits<double>::quiet_NaN();
}

  const unsigned int mapped_sector = (5U + N_SECTORS - (sector % N_SECTORS)) % N_SECTORS;

  const double phi = ((side == 1U ? 1.0 : -1.0) * (cdb_phi - M_PI / 2.0))
                     + (static_cast<double>(mapped_sector) * M_PI / 6.0);

  return wrap_phi(phi);
}
int IdealPadMap::get_layer_from_fee_channel(const unsigned int fee,
                                            const unsigned int channel) const
{
  if (fee >= N_FEE || channel >= N_CH) { return -1;
}
  const unsigned int key = 256U * fee + channel;
  return m_layer_by_key[key];
}
