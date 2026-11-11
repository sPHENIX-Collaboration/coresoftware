# TPC forward transport and lookup tables

This package provides explicitly enabled forward electron transport and adaptive
LUTs. Ordinary PHGarfield initialization and reverse reconstruction retain their
upstream behavior.

Build/install `offline/packages/PHGarfield`, then this package, into the same
prefix. Use an out-of-source `autogen.sh --prefix=...`, `make`, `make install`.
Developers can run `make check` to build and execute the regression tests; those
test files and executables are not needed in an installed workflow. Headers are
installed under `<tpcforwardtransport/...>`.

- `libTpcForwardTransport`: direct transport, LUT reading/interpolation and
  required configuration/integrity checks.
- `libTpcForwardTransportTools`: adaptive construction, packing and shard
  assembly. The runtime library does not link this construction library.
- `PHGarfieldForwardMapGenerator`: generates a ready-to-use adaptive LUT.
- `PHGarfieldCompactMapMerger`: assembles parallel-generated LUT shards.
- `tests/`: synthetic regressions, built only by `make check`.

The generator performs adaptive construction in memory and writes compact runtime
arrays directly. There is one file format and no intermediate adaptive ROOT file,
standalone converter or conversion-proof file. Construction retains the full
27-point stencil and error/reach checks; compact storage preserves the resulting
leaf classifications and usable corner values.

`PHGarfieldForwardTransport` inherits gas/field configuration from PHGarfield and
owns forward-specific initialization/integration. Call
`SetForwardTransportMode(true)` before initialization to enable it. Independent
physical holdouts, convergence studies, benchmarks and plots are separate QA
activities.

See [generation](PHGarfieldForwardMap.md) for commands/physical conventions and
[storage](PHGarfieldCompactMap.md) for the runtime format/interpolation policies.
