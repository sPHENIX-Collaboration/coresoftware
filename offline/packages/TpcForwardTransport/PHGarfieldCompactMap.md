# Compact forward-transport LUTs

Compact shard manifests must list shards in ascending `first_root` order,
with adjacent, nonoverlapping ranges. The merger preserves the exact manifest
text as a checksum-bound input identity and rejects out-of-order input.

`PHGarfieldCompactMap` converts adaptive schema-4 maps into schema-7 runtime
arrays. It retains cell routing, classifications and the exact double-precision
corner values needed for interpolation, while omitting generation-only samples
and diagnostics. Conversion changes storage, not the transport model or
interpolation order. It requires OpenSSL development headers and `libcrypto`
for SHA-256 file checks.

## Convert a complete adaptive map

```sh
PHGarfieldCompactMapConverter \
  --input adaptive.root --output compact.root --proof-output conversion.json \
  --expect-source-sha256 SOURCE_SHA256 \
  --expect-source-configuration CONFIGURATION_FNV1A64 \
  --expect-source-content CONTENT_FNV1A64
```

Pin the source file SHA-256 and the configuration/content fingerprints from the
validated adaptive input. The converter refuses changed inputs and existing
outputs. `--help` lists optional point/node/corner/value/byte ceilings.

Conversion validates the source, checks every cell bound, branch route and
retained corner value, writes a temporary ROOT file, then independently
reopens it and verifies the conversion. The JSON record identifies source and
output hashes, conversion software, resource limits and validation counts.
This verifies conversion and file structure; source drift accuracy and physical
coverage require independent validation.

## Runtime format

The file contains:

- `GarfieldFullTransport/Metadata`: compact identity and source bindings.
- `GarfieldFullTransport/SourceMetadata`: original physical configuration.
- `GarfieldFullTransport/Runtime`: primitive node, corner and value arrays.

Node and corner references use checked unsigned 32-bit indices. Corner values
remain three doubles: delta-r in cm, source-radius-weighted delta-phi in cm,
and absolute drift time in ns. Explicit count and byte limits are checked
before allocation; loaded arrays remain usable after closing the ROOT file.

Schema 8 packages compact arrays assembled from multiple adaptive shards.
It uses the same runtime array layout and interpolation as schema 7, with a
source manifest binding all converted pieces. Pieces must have matching
physical/refinement configuration and disjoint, contiguous root coverage; a
partial collection cannot be used as a complete simulation LUT. Generic
conversion/merging tools can combine pieces without building a complete
training map in memory.

A simulation consumer must explicitly select full-forward transport and supply
its expected source bindings and runtime-file hash from the conversion record.
Unsupported or malformed files are errors. Certified interpolation remains
unchanged. Schema 7 additionally retains corner records for eligible
`Uncertified` and `NumericallyUnresolved` leaves when all eight values are
finite, reaching and angularly consistent. Mixed and sampled-all-lost leaves
do not retain interpolation records. Missing or unsafe support stays absent.

Retaining corners does not certify a leaf or promise the interpolation guard.
The strict `Interpolate` API still rejects every noncertified leaf. An explicit
`InterpolateRelaxed` call can use eligible retained records and reports whether
the source was certified. Consumers must choose and report their treatment of
noncertified leaves: direct drift, rejection or an explicitly approximate
discard/interpolation policy. A sampled-all-lost leaf is not a proof of loss
throughout its volume.

Readers also support legacy schemas 5/6, which store corners only for certified
leaves. `--legacy-certified-only yes` on the converter writes schema 5; the
same option in merger `convert` mode writes schema 6. Merger `merge` mode
preserves a homogeneous input version and rejects mixed legacy/new versions.
Legacy layouts and canonical fingerprint/proof bytes remain unchanged.

The new versions keep the 12-byte node, 32-byte corner-record and 24-byte value
sizes. They allow eligible noncertified leaf nodes to reference a corner record
and add a noncertified-leaf conversion-proof count. Layout identifiers and
versioned fingerprints distinguish this extension. Converting an existing
schema-5/6 file cannot recover values it omitted; reconvert its adaptive source
or regenerate a physically matching map.

## Tests

`make check` runs `test_compact_map` with synthetic adaptive maps and small
ROOT files. It verifies exact routing/corner values, periodic and boundary
queries, resource limits, corrupted data rejection and independent readback.
It also covers retained noncertified corners, strict versus relaxed lookup,
legacy output, shard offset handling and mixed-version rejection.

## Convert pieces before merging

`PHGarfieldCompactMapMerger` has two modes. `convert` reads adaptive schema-4
shards and converts them into a compact collection. `merge` combines compact
collections. Each input list has three tab-separated columns per line:
absolute ROOT filename, file SHA-256 and content fingerprint. Supply pieces in
increasing root order, without a header:

```sh
PHGarfieldCompactMapMerger --mode convert --input-list shards.tsv \
  --expect-input-sha256 LIST_SHA256 --output pieces.root \
  --summary pieces.json --complete no
PHGarfieldCompactMapMerger --mode merge --input-list pieces.tsv \
  --expect-input-sha256 LIST_SHA256 --output runtime.root \
  --summary runtime.json --complete yes
```

Use `--complete yes` only when the result covers all coarse roots. The summary
contains the runtime hash and source bindings needed by the simulation reader.
Conversion preserves every retained corner double and checks every remapped
reference. File structure and exact conversion are checked separately from the
accuracy of the original generated transport.
