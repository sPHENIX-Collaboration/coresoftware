# Adaptive LUT runtime storage

Adaptive refers to variable cell sizes. Compact refers to storing only the
routing tree, leaf classifications, usable corner references and unique corner
values needed for lookup. The generator builds samples/refinement diagnostics in
memory and writes this runtime representation directly.

## One output format

Schema 9 is named `PHG4TpcGarfieldAdaptiveTransportLUT`:

```text
GarfieldFullTransport/
  Metadata/   physical configuration, payload identities, generation totals
  Runtime/    Nodes, Corners, Values, Roots and refinement/geometry metadata
```

There is no separate adaptive training file, `SourceMetadata` directory,
conversion executable, external source binding or conversion-proof JSON. Earlier
experimental schema versions are not supported by this package. Existing private
files can still be used with their matching historical software.

`Metadata` records units, coordinate conventions, integration step, gas/fields,
geometry, a canonical physical configuration and its FNV-1a fingerprint. Each
nonzero `max_radial_field_extension_cm` is recorded and fingerprinted; its absence
in earlier schema-9 files means strict radial coverage. Each
shard also records its root interval and execution counters. Merging sums query,
reuse and leaf counters, takes the maximum drift duration and preserves common
physical/refinement metadata. It does not claim to have generated a new set of
direct samples.

`Runtime/Nodes` stores binary split children or a leaf's class/reason and corner
record. `Corners` contains eight indices into `Values`; each value holds radial
and source-radius-weighted angular displacements plus drift duration. Shared
sample values are stored once per generated shard. `Roots` locates the coarse
cells. The exact packing is internal to `PHGarfieldCompactMap`.

A deterministic content fingerprint covers geometry, configuration/refinement
identity and primitive arrays. The reader verifies it together with array types,
counts, references, topology and allocation limits. These are consistency checks,
not proof of physical accuracy. ROOT UUID/timestamps are not physical identities.

## Interpolation and noncertified cells

Strict lookup interpolates certified leaves. The original 27-point
classification is preserved for every leaf; packing never promotes a cell.

For the explicit fast approximation policy, numerical/budget-limited leaves may
also be interpolated when all required corner values reach the readout and have
a usable angular branch. Those values are retained in the file. This deliberately
relaxes the sampled guard; classification/reason remains available to the caller.

Mixed and sampled-all-lost leaves do not provide an interpolated arrival.
Sampled-all-lost is not a proof that every point is physically lost. The consumer
must explicitly select direct fallback, rejection, or counted policy discard;
policy discard must not be reported as a calculated physical loss. Malformed
files and invalid calculations remain errors.

The merger accepts the same runtime format for individual shards and partial
collections. A complete result requires contiguous coverage of the entire grid.
See [generation](PHGarfieldForwardMap.md) for commands. Production simulation
should load only complete tables and verify their physical configuration against
its own detector, drift and magnetic-field settings.
