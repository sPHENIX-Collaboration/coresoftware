# TPC forward-transport LUTs

`PHGarfieldForwardMapGenerator` tabulates deterministic mean electron transport
from an ionization position to a TPC readout plane. It uses the explicit-Euler
`PHGarfield::ForwardDrift` integrator. This does not enable forward transport in
simulation or change the existing reverse-drift reconstruction.

## Generate a regular map

After building and installing PHGarfield, run:

```sh
PHGarfieldForwardMapGenerator --output forward.root \
  --electric-field-map corrections.root \
  --gas-file detector.gas --magnetic-field-map field.root \
  --global-tag YOUR_CDB_TAG --timestamp YOUR_TIMESTAMP \
  --step-ns 12.5 --phi-nodes 36 --radial-nodes 24 --longitudinal-nodes 28
```

The payload paths and numerical/geometry settings must match the intended
simulation. Gas and magnetic-field inputs may instead be resolved from CDB.
Use `--help` for geometry, placement, time calibration and field options.
Side 0 is south (negative z); side 1 is north (positive z).

Forward 2D and 3D space-charge maps use the longitudinal coordinate
`max(0, abs(z) - central_membrane_half_width_cm)`, measured from the membrane's
gas-facing surface. Frame-field maps and upstream reverse reconstruction
retain their existing `abs(z)` lookup convention. LUT source coordinates and
stored displacements remain in the physical local TPC frame.

New 3D LUTs record this mapping in `space_charge_z_coordinate_policy` and in
the canonical configuration fingerprint. Regenerate older 3D tables before
using them with this model; a missing or different policy is rejected.
New 2D tables record the same policy as descriptive metadata, retaining their
established fingerprint and compatibility with older 2D tables.

For per-side 3D space charge, replace `--electric-field-map` with the paired
`--electric-field-map-3d-side0` and `--electric-field-map-3d-side1` options.
Optional frame maps likewise require both sides. Their coefficient is selected
explicitly with `--frame-charge-scale` (default zero). Nominal E/B-only maps
require `--allow-no-electric-corrections --scale-side0 0 --scale-side1 0`.
Use `--zero-b-field` to declare B=0 explicitly.

The generator validates payloads, magnetic-field and space-charge-map coverage, geometry and the
correction-disabled Garfield/PHG4 timing reference. It writes a temporary ROOT
file, checks it, then publishes the output. Existing outputs are preserved
unless `--force` is supplied.
Space-charge maps must cover the configured radii and gas-surface longitudinal
range from zero to `readout_abs_z_cm - central_membrane_half_width_cm` on both
sides. Coverage uses histogram edges; interpolation between an edge and its
nearest bin center retains the boundary value. Axis coverage does not establish
the payload's physical coordinate origin, which must match the declared model.

## Stored quantities and interpolation

The `GarfieldFullTransport` directory stores separate positive/negative-z maps:

| Quantity | Definition | Units |
| --- | --- | --- |
| `hDeltaR_*` | readout radius minus source radius | cm |
| `hDeltaRPhi_*` | source radius times wrapped readout-minus-source phi | cm |
| `hDriftTime_*` | absolute drift duration | ns |
| `hReachesReadoutNode_*` | physical node reachability | Boolean |
| `hReachCellClass_*` | all lost, mixed, all reaching, or numerically unresolved | code |

The reader uses trilinear interpolation within supported cells and wraps phi
periodically. Schema 1 records axisymmetric inputs; schema 2 adds explicit
3D/frame configuration. A positive `--interpolation-tolerance-um` selects schema
3 and masks cells whose sampled transverse or time-equivalent interpolation
error exceeds that value. These finite probes do not prove an error bound over
an entire cell; validate independent points and boundaries for a new grid.

Metadata records geometry, integration settings, nominal timing, field/gas
payload identities and configuration fingerprints. Physical loss, unresolved
interpolation and invalid input are distinct outcomes; invalid transport must
not silently become a legacy displacement correction.

## Simulation integration

Full forward transport requires explicit enablement in the simulation consumer.
It supplies the deterministic endpoint and absolute drift time; simulation
retains its configured stochastic diffusion and electronics response. Arrival
time is ionization time plus LUT drift duration plus stochastic time diffusion.
Legacy displacement maps must not be added to the same full transport.

## Tests

`make check` runs `test_forward_map_formula`, a self-contained check of units,
angle wrapping, configuration validation and canonical metadata. Integration
validation additionally requires gas/field payloads and direct-drift holdouts.

## Adaptive generation

`--adaptive` selects schema 4. It uses the same pinned physics and direct
integration as regular generation, with a binary tree that bisects one cell
axis at a time. The selected axis follows the measured directional transport
curvature or reaching/lost transitions. Reaching leaves interpolate only when
all 27 stencil samples reach and the sampled transverse and time-equivalent
errors pass `--interpolation-tolerance-um` (adaptive default 10 micrometres).

The stored primitive `Points`, `Cells` and `Roots` arrays include classifications,
root ownership, support values, refinement policy and content fingerprints.
Depth, per-root query/node budgets and per-file point/cell/byte limits bound
resource use. Budget exhaustion leaves an explicitly uncertified cell; it does
not turn the cell into an interpolation region. Physical mixed and sampled
all-lost leaves also require a direct query when used by an enabled fallback.
Finite lost samples cannot certify that a whole volume is physically lost.

`--adaptive-root-first` and `--adaptive-root-end` can generate disjoint coarse
root ranges. Partial maps are marked as shards and cannot be used directly in
simulation; combine them into a complete map or compact collection first.
The generator validates and independently reopens its output before publication.

`make check` additionally runs `test_adaptive_map` using synthetic samplers and
small ROOT files. It checks refinement, periodic support, classification,
malformed arrays, complete partitions and ROOT readback without detector inputs.

## Transport implementation

`PHGarfield` owns gas and field resources and performs deterministic forward
integration. Internal `PHGarfieldForwardSetup` helpers select input sources,
resolve standalone/CDB policy, and preserve per-side map precedence; they do
not load transport resources or evaluate fields. `PHGarfieldFieldMapValidation`
inspects ROOT Cartesian grid structure and returns validated domain bounds.
The helpers are not installed public APIs. Ordinary initialization, reverse
drift and electric/magnetic-field evaluation remain in `PHGarfield`.
