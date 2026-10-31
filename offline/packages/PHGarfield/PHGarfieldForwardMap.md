# TPC forward transport

PHGarfield exposes explicitly enabled deterministic forward electron transport.
Gas and field inputs, drift geometry and side coefficients must be configured
before initialization. ForwardDrift returns a readout endpoint, absolute drift
duration and an explicit reaching/loss/error status. Transport stops at the
field cages. Ordinary initialization and reverse reconstruction retain upstream
behavior. Self-contained configuration, field and coverage checks run with
`make check`.

## Transport implementation

`PHGarfield` owns gas and field resources and performs deterministic forward
integration. Internal `PHGarfieldForwardSetup` helpers select input sources,
resolve standalone/CDB policy, and preserve per-side map precedence; they do
not load transport resources or evaluate fields. `PHGarfieldFieldMapValidation`
inspects ROOT Cartesian grid structure and returns validated domain bounds.
The helpers are not installed public APIs. Ordinary initialization, reverse
drift and electric/magnetic-field evaluation remain in `PHGarfield`.
