#include "PHGarfieldForwardSetup.h"

#include <cstddef>

namespace PHGarfieldForwardSetup
{
  Selection SelectInputs(const Request& request)
  {
    Selection selection;
    const bool explicit_side_coefficients = request.space_charge_scale_override[0] &&
                                            request.space_charge_scale_override[1];
    if (!request.has_node_tree && !explicit_side_coefficients)
    {
      return selection;
    }

    selection.valid = true;
    selection.use_run_dependent_inputs = request.has_node_tree;
    selection.load_tpc_fee_map = request.load_tpc_fee_map;
    selection.magnetic_field_override = request.magnetic_field_override;
    if (request.zero_magnetic_field)
    {
      selection.magnetic_field_source = MagneticFieldSource::Zero;
    }
    else if (!request.magnetic_field_override.empty())
    {
      selection.magnetic_field_source = MagneticFieldSource::ExplicitFile;
    }
    else if (request.has_node_tree)
    {
      selection.magnetic_field_source = MagneticFieldSource::NodeTree;
    }

    // A change to 2D must discard previously resolved automatic 3D paths.
    // Explicit side files retain their existing per-side 3D precedence.
    selection.electric_field_maps_3d = request.electric_field_maps_3d;
    if (request.electric_field_map_override || request.use_2d_electric_field_map)
    {
      for (std::size_t side = 0; side < selection.electric_field_maps_3d.size(); ++side)
      {
        if (!request.electric_field_map_3d_override[side])
        {
          selection.electric_field_maps_3d[side].clear();
        }
      }
    }

    selection.frame_electric_field_maps_3d = request.frame_electric_field_maps_3d;
    const bool suppress_automatic_frames = request.frame_electric_field_map_override ||
                                           (request.electric_field_map_override &&
                                            !request.frame_electric_field_map_3d_override[0] &&
                                            !request.frame_electric_field_map_3d_override[1]);
    if (suppress_automatic_frames)
    {
      for (std::size_t side = 0; side < selection.frame_electric_field_maps_3d.size(); ++side)
      {
        if (!request.frame_electric_field_map_3d_override[side])
        {
          selection.frame_electric_field_maps_3d[side].clear();
        }
      }
    }
    return selection;
  }

  std::string ResolveMagneticField(const Selection& selection, const CdbLookup& lookup)
  {
    switch (selection.magnetic_field_source)
    {
    case MagneticFieldSource::Zero:
      return "ZERO_UNIFORM_FIELD";
    case MagneticFieldSource::ExplicitFile:
      return selection.magnetic_field_override;
    case MagneticFieldSource::NodeTree:
      return "FIELD_MAP node";
    case MagneticFieldSource::CdbFile:
      return lookup("FIELDMAP_TRACKING");
    }
    return {};
  }

  std::string ResolveFeeMap(const Selection& selection, const CdbLookup& lookup)
  {
    return selection.load_tpc_fee_map ? lookup("TPC_FEE_CHANNEL_MAP") : std::string{};
  }

  GasSelection ResolveGas(const GasRequest& request,
                          const CdbLookup& lookup,
                          const FileExists& exists)
  {
    GasSelection selection;
    if (!request.override_filename.empty())
    {
      selection.source = GasSource::ExplicitFile;
      selection.filename = request.override_filename;
      if (!exists(selection.filename))
      {
        selection.valid = false;
        return selection;
      }
    }
    else
    {
      selection.filename = lookup("PHGARFIELD_GAS");
    }

    if (selection.filename.empty() || !exists(selection.filename))
    {
      selection.source = GasSource::DefaultFile;
      selection.unavailable_cdb_filename = selection.filename;
      selection.filename = request.default_filename;
    }
    return selection;
  }
}  // namespace PHGarfieldForwardSetup
