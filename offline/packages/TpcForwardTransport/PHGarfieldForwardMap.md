# TPC forward-transport LUTs

`PHGarfieldForwardMapGenerator` tabulates deterministic mean electron transport
from an ionization position to a TPC readout plane. It uses the explicit-Euler
`PHGarfieldForwardTransport::ForwardDrift` integrator and generates adaptive schema-4 tables.
The command can generate one complete table or a coarse-root-range shard.
`--adaptive` remains accepted for existing campaign commands; adaptive
refinement is always used. Regular-grid generation and full adaptive-map merging
are outside this workflow.

## Generate an adaptive map

After building and installing PHGarfield and TpcForwardTransport, run:

```sh
PHGarfieldForwardMapGenerator --output adaptive.root \
  --electric-field-map-3d-side0 south.root \
  --electric-field-map-3d-side1 north.root \
  --gas-file detector.gas --magnetic-field-map field.root \
  --global-tag YOUR_CDB_TAG --timestamp YOUR_TIMESTAMP \
  --step-ns 12.5 --phi-nodes 36 --radial-nodes 24 --longitudinal-nodes 28 \
  --interpolation-tolerance-um 10
```

The payload paths, side coefficients and numerical/geometry settings must match
the intended simulation. Gas and magnetic-field inputs may be resolved from
CDB. Side 0 is south (negative z); side 1 is north (positive z). Both 3D
space-charge maps must be specified. Optional 3D frame maps also require both
sides and an explicit `--frame-charge-scale` (default zero). Use `--zero-b-field`
to declare B=0. `--help` lists all physical inputs and resource limits.

Forward 2D/3D space-charge lookup uses
`max(0, abs(z) - central_membrane_half_width_cm)`, measured from the membrane's
gas-facing surface. Frame-field maps and reverse reconstruction retain their
existing `abs(z)` convention. LUT source coordinates and stored displacements
remain in the physical local TPC frame. The 3D coordinate policy is recorded
in metadata and bound into the canonical configuration fingerprint.

Selected space-charge maps must cover the configured radial range and gas
longitudinal range from zero to
`readout_abs_z_cm - central_membrane_half_width_cm` on both sides. Coverage uses
histogram edges; between an edge and its nearest bin center, interpolation
retains the boundary value. Axis coverage alone does not establish the
payload's physical coordinate origin. Magnetic-field coverage, input payloads,
geometry and the correction-disabled Garfield/PHG4 timing reference are also
checked before generation.

## Refinement and interpolation

The coarse grid is refined by bisecting one cell axis at a time. The selected
axis follows sampled directional transport curvature or reaching/lost
transitions. A reaching leaf permits trilinear interpolation only when all 27
stencil samples reach and the sampled transverse and time-equivalent errors
pass `--interpolation-tolerance-um` (default 10 micrometres). Finite probes do
not prove an error bound throughout a cell; validate independent points and
boundaries for a new table.

The schema-4 `Points`, `Cells` and `Roots` arrays preserve the support values,
classifications, refinement policy and content fingerprints. Each reaching
support value contains delta-r in cm, source-radius-weighted delta-phi in cm,
and absolute drift duration in ns. Lookup wraps phi periodically. Metadata
retains established configuration identifiers and fixed probe/chunk fields
for compatibility with existing adaptive/compact source bindings.

Depth, per-root query/node budgets and per-file point/cell/byte limits bound
resource use. Exhausted budgets leave explicitly uncertified cells. Physical
mixed and sampled all-lost leaves are also noncertified; finite lost samples do
not certify an entire volume as lost. Strict use requires direct drift or
rejection for noncertified leaves. An explicitly approximate consumer may
discard selected leaves or use relaxed interpolation with usable support.
Conversion preserves source classifications and does not upgrade them.

The generator writes a temporary ROOT file, independently reopens and validates
it, then publishes the output. Existing outputs are preserved unless `--force`
is supplied.

## Complete maps and parallel shards

For a small complete adaptive map, use `PHGarfieldCompactMapConverter` to create
a compact schema-7 file. For large campaigns, use `--adaptive-root-first N` and
`--adaptive-root-end N` to generate disjoint coarse-root ranges, convert those
pieces first, and merge compact collections with `PHGarfieldCompactMapMerger`.
This avoids constructing a complete detailed adaptive map in memory. Partial
shards are explicitly marked and cannot be consumed as complete simulation
LUTs. New schema-7/8 compact files retain usable noncertified corner values for
explicit relaxed interpolation, while legacy schema-5/6 output remains
available. See `PHGarfieldCompactMap.md` for commands, version compatibility
and manifest requirements.

## Simulation integration

Simulation consumption requires explicit full-forward enablement in the
separate consumer. The LUT supplies the deterministic endpoint and absolute
drift duration. Simulation retains stochastic diffusion and electronics;
arrival time is ionization time plus LUT duration plus stochastic time
diffusion. Existing reverse reconstruction retains its configuration.

## Tests and implementation

`make check` runs eight self-contained tests covering input policy, field
coordinates/coverage, node encoding, metadata, adaptive refinement, compact
conversion, payload identities and malformed ROOT data. Real gas/field
validation additionally requires independent direct-drift holdouts.

`PHGarfield` supplies the existing gas/field resources. The explicitly enabled
`PHGarfieldForwardTransport` subclass performs forward initialization and direct
transport in this package; ordinary PHGarfield retains upstream behavior.
`PHGarfieldForwardSetup` selects input sources without loading resources or
evaluating fields. `PHGarfieldFieldMapValidation` inspects Cartesian magnetic
map structure. These helpers are internal. `PHGarfieldForwardMap` supplies
shared configuration, physics checks, metadata and encoding to the adaptive
and compact implementation; it does not generate regular-grid tables.
