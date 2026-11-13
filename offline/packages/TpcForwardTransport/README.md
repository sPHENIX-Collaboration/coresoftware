# TPC forward transport and lookup tables

Generate and read adaptive TPC forward-drift LUTs. Forward transport is opt-in;
ordinary PHGarfield initialization and reverse reconstruction are unchanged.

Build/install `offline/packages/PHGarfield`, then this package, into the same
prefix, using an out-of-source `autogen.sh --prefix=...`, `make`, `make install`.
Headers are installed under `<tpcforwardtransport/...>`.

- [Generate a LUT](PHGarfieldForwardMap.md): commands, inputs and parallel jobs.
- [File format](PHGarfieldCompactMap.md): optional details for developers.

`PHGarfieldForwardMapGenerator` writes the runtime ROOT format directly;
`PHGarfieldCompactMapMerger` combines shards.

The runtime library is `libTpcForwardTransport`; generation code is in the separate
`libTpcForwardTransportTools`. Tests build and run only with `make check` and are
not installed. Independent accuracy studies and benchmarks are separate QA.

For direct transport, enable `PHGarfieldForwardTransport` with
`SetForwardTransportMode(true)` before initialization.
