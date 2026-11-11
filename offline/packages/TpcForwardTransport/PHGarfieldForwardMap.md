# Generate an adaptive forward-drift LUT

The standalone executable needs no Fun4All event loop. After setting up sPHENIX
and the installed package, select a CDB tag and timestamp containing the payloads:

```sh
PHGarfieldForwardMapGenerator \
  --global-tag YOUR_TAG --timestamp YOUR_TIMESTAMP \
  --nominal-e-vcm 375 --scale-side0 0.579394 --scale-side1 0.610215 \
  --frame-charge-scale -180 --step-ns 1.5625 \
  --output runtime.root
```

Explicit paths override CDB independently for each side. For example:

```sh
PHGarfieldForwardMapGenerator \
  --electric-field-map-3d-side0 /path/spacecharge_side0.root \
  --electric-field-map-3d-side1 /path/spacecharge_side1.root \
  --gas-file /path/gas.gas \
  --magnetic-field-map /path/fieldmap.root \
  --nominal-e-vcm 375 --scale-side0 0.579394 --scale-side1 0.610215 \
  --step-ns 1.5625 --interpolation-tolerance-um 10 \
  --output runtime.root
```

The fields, scales and step above are examples, not prescribed calibrations. Select
compatible physical inputs and establish numerical convergence for your study.
`--help` lists all options. Optional per-side frame maps use
`--frame-field-map-3d-side0/side1` and `--frame-charge-scale`.

Omitting `--gas-file` or `--magnetic-field-map` resolves `PHGARFIELD_GAS` or
`FIELDMAP_TRACKING` through CDB using `--global-tag` and `--timestamp` (defaults:
`FieldMapTest`, `1`). Missing electric-map paths resolve
`Tpc_PolySeeding_EField3D_Side0/Side1`. A nonzero frame coefficient or any explicit
frame path requests the frame family; missing sides then resolve
`Tpc_PolySeeding_ModuleFrames3D_Side0/Side1`. With coefficient zero and no frame
paths, frame CDB lookup is skipped. Both sides of each requested family are
required; an absent CDB entry fails with its payload name. The selected tag and
timestamp must contain those entries; the defaults do not guarantee this.

An explicit gas path may be a single file or a directory of
ordered `MERGED_E###.gas` files. Resolved paths and payload MD5 values are stored.
The tag does not choose scale coefficients, nominal electric field, geometry or
numerical settings. `--zero-b-field` selects zero magnetic field explicitly.
Inputs must cover the configured gas volume.

## Physical and numerical definitions

Coordinates are in the placed TPC's local frame. The stored displacement is
`delta_r_cm`, `r0 * delta_phi_rad`, plus absolute drift duration in ns. Endpoint
azimuth is `phi0 + stored_angular_value/r0`; it is not divided by endpoint radius.
The gas extends from `central_membrane_half_width_cm` to `readout_abs_z_cm` on
each side. Forward space-charge maps use
`max(0, abs(z)-central_membrane_half_width_cm)` longitudinally. Frame maps retain
their own coordinate convention. This does not change reverse reconstruction.

Direct samples use explicit Euler integration with `--step-ns`. The generator
checks gas/geometry/field consistency and a scale-zero timing reference against
the declared PHG4 calibration velocity. These checks do not establish numerical
convergence; independent finer-step comparisons remain necessary.

A coarse cylindrical grid is refined anisotropically. Each candidate cell uses
27 corner/midpoint stencil positions. Certification requires all samples to
reach the readout, a consistent angular branch, and sampled transverse endpoint
and nominal-velocity-scaled timing errors within the guard. This is a sampled
maximum check, not a mathematical bound between samples. Error discrepancies or
reach-status transitions choose the split axis. Largest unresolved cylindrical
volume is refined first within each root; stable IDs break ties.

`--adaptive-max-depth`, `--adaptive-max-nodes-per-root`, and
`--adaptive-max-queries-per-root` limit refinement. Global point, cell and estimated
allocation caps are also explicit; these do not measure aggregate peak RSS.
Budget-limited cells keep their classification;
they are not silently certified. Physical reach loss is distinct from lack of
interpolation certification.

## Parallel generation

Each job can generate a nonempty coarse-root interval:

```sh
PHGarfieldForwardMapGenerator [identical physical/grid/refinement options] \
  --adaptive-root-first 0 --adaptive-root-end 100 --output shard0.root
```

Repeat for contiguous intervals covering every root. Root count is
`2 * phi_nodes * (radial_nodes-1) * (longitudinal_nodes-1)`. Each output already
contains the runtime format; no conversion is needed. A partial file must not be
used as a complete simulation LUT.

List shard paths in ascending contiguous root-range order, one per line, then:

```sh
PHGarfieldCompactMapMerger --input-list shards.txt --output runtime.root
```

Relative paths are resolved against the input-list directory. Missing,
overlapping, duplicate, out-of-order or incompatible shards are rejected. The
default requires full coverage. For hierarchical merges of smaller groups use
`--complete no`, then merge the partial collections into the complete LUT.
Resource checks reject arrays exceeding supported allocation/index limits.

Outputs are written to temporary files, closed, independently reopened and
checked before publication. Existing outputs are preserved; the generator's
`--force` explicitly permits replacement. The merger never overwrites.
