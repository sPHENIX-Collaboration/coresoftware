#ifndef IDEALPADMAP_H
#define IDEALPADMAP_H

#include <array>
#include <string>
#include <vector>

class CDBTTree;

class IdealPadMap
{
 public:
  IdealPadMap();
  ~IdealPadMap();

  int load_from_cdb(int verbosity = 0);
  bool is_loaded() const { return m_is_loaded; }

  // TPC layers are real layer numbers: 7..54.
  // region is 0,1,2 for inner/mid/outer.
  int get_region(unsigned int layer) const;
  unsigned int get_pads_per_sector(unsigned int region) const;
  unsigned int get_pads_per_sector_for_layer(unsigned int layer) const;
  unsigned int get_total_phibins(unsigned int layer) const;

  // Radius from CDB, averaged over all pads in that layer. Returned in cm.
  double get_radius(unsigned int layer) const;
  double get_local_radius(unsigned int layer) const { return get_radius(layer); }
  double get_layer_thickness(unsigned int layer) const;

  // Raw local CDB phi for one pad inside one sector.
  // local_phibin is 0..pads_per_sector(layer)-1.
  double get_cdb_local_phi(unsigned int layer, unsigned int local_phibin) const;

  // Convert a full-circle PHG4TpcGeom phibin back to global phi.
  // phibin is the integer returned by layergeom->get_phibin(phi, side).
  double get_phi(unsigned int side, unsigned int layer, unsigned int phibin) const;

  // Same conversion, but with sector already known and local bin inside sector.
  double get_phi(unsigned int side,
                 unsigned int sector,
                 unsigned int layer,
                 unsigned int local_phibin) const;

  // Optional helper if you still need direct FEE/channel lookup.
  int get_layer_from_fee_channel(unsigned int fee, unsigned int channel) const;

 private:
  static constexpr unsigned int N_SIDES = 2;
  static constexpr unsigned int N_SECTORS = 12;
  static constexpr unsigned int N_REGIONS = 3;
  static constexpr unsigned int N_LAYERS = 48;
  static constexpr unsigned int FIRST_LAYER = 7;
  static constexpr unsigned int LAST_LAYER = 54;
  static constexpr unsigned int N_FEE = 26;
  static constexpr unsigned int N_CH = 256;

  int layer_index(unsigned int layer) const;
  double wrap_phi(double phi) const;

  bool m_is_loaded = false;

  // Indexed by layer-7.  Each layer vector is one sector worth of CDB pads.
  std::array<std::vector<double>, N_LAYERS> m_cdb_phi_by_layer;
  std::array<double, N_LAYERS> m_radius_cm_by_layer;

  // Indexed by 256*fee + channel. -1 means no valid TPC layer.
  std::array<int, N_FEE * N_CH> m_layer_by_key;
};

#endif
