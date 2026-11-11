# TPC forward transport and lookup tables

This package contains explicitly enabled forward electron transport and the
adaptive/compact LUT workflow. The upstream PHGarfield initialization and reverse
reconstruction remain the default in the PHGarfield package.

Build and install `offline/packages/PHGarfield` first, then this package, into
the same installation prefix. Build the simulation consumer afterward if needed.
Use the usual out-of-source `autogen.sh --prefix=...`, `make`, `make check` and
`make install` workflow. Installed headers use `<tpcforwardtransport/...>`.

## Components

- `libTpcForwardTransport`: LUT formats, readers, interpolation, required
  integrity checks, and `PHGarfieldForwardTransport` for optional direct drift.
- `libTpcForwardTransportTools`: adaptive construction, generation, conversion
  and shard collection support. The runtime library does not link this library.
- `tools/`: the three command-line programs and construction implementations.
  Existing executable names and command-line options are preserved.
- `tests/`: small synthetic regressions, built and run only by `make check`.
  They are excluded from both production libraries.

`PHGarfieldForwardTransport` inherits PHGarfield's gas/field configuration and
owns forward-specific initialization and integration. Select forward transport
explicitly before initialization. Constructing ordinary PHGarfield retains
upstream behavior. The forward result/status values, integration order and LUT
coordinate conventions are unchanged by this package separation.

The 27-point certification checks remain part of adaptive construction.
Malformed-file and configuration checks remain part of the runtime reader.
Independent physical holdouts, performance campaigns and presentation plots
are separate QA activities and are not linked into these libraries.

See [LUT generation](PHGarfieldForwardMap.md) and
[compact formats and conversion](PHGarfieldCompactMap.md) for commands,
units, version compatibility and interpolation policies.
