# Generate a forward-drift LUT

After [building the package](README.md), use the standalone generator. No Fun4All
event loop is needed. `PHGarfieldForwardMapGenerator --help` lists all options.

## Example with available SDCC inputs

This recipe uses pinned CVMFS gas, electric/frame and GAP magnetic maps. Use the
same magnetic payload and scale in Geant. It keeps the full gas volume
`r=21.6–76.4 cm`, CM half-width `0.28 cm` and readout at `|z|=102.605 cm`.

```sh
cdb=/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb
PHGarfieldForwardMapGenerator \
  --gas-file "$cdb/PHGARFIELD_GAS/2f/77/2f776187989ee6fbd33e7f85e6cd5099_Ar75_CF20_iso5_v2.gas" \
  --magnetic-field-map "$cdb/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root" \
  --electric-field-map-3d-side0 "$cdb/Tpc_PolySeeding_EField3D_Side0/ba/76/ba7618a6247314a270a9ee4b95530f39_phgarfield_ibf_3d_map_gainOn_alpha_1p8_flatphi_side0.root" \
  --electric-field-map-3d-side1 "$cdb/Tpc_PolySeeding_EField3D_Side1/6b/fd/6bfd3df3e1ef1dc50931e28171240431_phgarfield_ibf_3d_map_gainOn_alpha_1p8_flatphi_side1.root" \
  --frame-field-map-3d-side0 "$cdb/Tpc_PolySeeding_ModuleFrames3D_Side0/05/0c/050c53906b9f27aa2cb69733a3b863df_phgarfield_tpc_moduleframes_3d_side0.root" \
  --frame-field-map-3d-side1 "$cdb/Tpc_PolySeeding_ModuleFrames3D_Side1/8d/c2/8dc290303f3c6caabf4a48fea0c52925_phgarfield_tpc_moduleframes_3d_side1.root" \
  --nominal-e-vcm 375 --scale-side0 0.579394 --scale-side1 0.610215 \
  --frame-charge-scale -180 --magnetic-field-scale 1 \
  --step-ns 1.5625 --interpolation-tolerance-um 10 \
  --output runtime.root
```

Append `--check-inputs` to check the setup without creating a LUT. A difference
between gas-model and reference drift times is reported as a warning; calculated
drift times are unchanged.

Radial space-charge-map gaps of up to **0.3 cm per boundary** use the nearest edge-bin
field value by default, with a startup message and a recorded LUT setting. Gas
boundaries stay unchanged. Set `--radial-field-extension-cm 0` for strict coverage.
This is a constant-field extension. Longitudinal space-charge-map and magnetic-map
coverage remain mandatory.

These are example physics inputs, not universal calibrations. Changing fields,
gas or geometry requires a corresponding LUT and numerical validation.

## Using CDB instead of explicit paths

`--global-tag TAG --timestamp VALUE` resolves omitted gas and magnetic paths from
`PHGARFIELD_GAS` and `FIELDMAP_TRACKING`; override the magnetic path for GAP-based
simulation. Missing electric paths use `Tpc_PolySeeding_EField3D_Side0/Side1`.
A nonzero frame coefficient requests `Tpc_PolySeeding_ModuleFrames3D_Side0/Side1`.
The tag must contain those entries; the default `FieldMapTest` tag does not
guarantee them. Explicit paths take precedence. A gas path may also be a directory
of ordered `MERGED_E###.gas` files.

## Parallel generation interface

The example above generates the whole table in one process. For large tables,
batch wrappers can divide the work and create the ordered shard list automatically.
The underlying interface uses identical input/settings options for each job and
a different contiguous coarse-root interval:

```sh
PHGarfieldForwardMapGenerator [same input/settings options] \
  --adaptive-root-first 0 --adaptive-root-end 100 --output shard0.root
```

The full range is `[0, N)`, where
`N = 2 * phi_nodes * (radial_nodes-1) * (longitudinal_nodes-1)`; the end is exclusive.
The wrapper writes shard paths in ascending range order, one per line in
`shards.txt`. Merge them with:

```sh
PHGarfieldCompactMapMerger --input-list shards.txt --output runtime.root
```

Only the complete merged LUT is suitable for simulation. For intermediate merges,
use `--complete no`. Generation has configurable refinement and allocation limits;
split into smaller root intervals if a job exceeds its allocation limit.
Existing outputs are refused; generator `--force` explicitly allows replacement.

The [format document](PHGarfieldCompactMap.md) describes stored quantities and
interpolation policies. The error guard checks sampled interpolation error;
independent finer-step comparisons are still needed to assess integration error.
