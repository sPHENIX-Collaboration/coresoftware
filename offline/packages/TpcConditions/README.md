# TPC Conditions

The `TpcConditions` package provides event-by-event access to TPC operating
conditions stored in the sPHENIX Conditions Database (CDB).

`TpcConditionsReco` retrieves the `TPC_CONDITIONS` payload for the current
run and uses the event BCO to select the most recent conditions record at or
before the event.  The resulting conditions are made available to other
Fun4All modules through a transient `TpcConditions` object on the node tree.

The `TpcConditions` object is a `PHDataNode` and is intended only for
communication between reconstruction modules.  It is not intended to be
written to an output DST.

## Event Selection

Two flags describe whether the TPC had recovered from an HV trip:

- `FieldOK` indicates that the complete GEM/transfer-gap voltage stack has
  recovered the proper total voltage required for the TPC drift field.
- `GainOK` indicates that the GEMs have recovered their appropriate operating
  voltages.

`TpcConditionsReco` returns `Fun4AllReturnCodes::ABORTEVENT` if either
`FieldOK` or `GainOK` is false.

## Available Conditions

The following quantities are available for the current event:

- Gas temperature
- Gas pressure
- `FieldOK`
- `GainOK`
- Median G4 load current for all undamaged GEMs
- Median G4 load current for the North side
- Median G4 load current for the South side
- Median G4 load current separately for R1, R2, and R3 on each side

The corresponding load-current accessors are:

```cpp
get_LoadCurrent()   // all undamaged GEMs
get_LoadNorth()     // North
get_LoadSouth()     // South

get_LoadNR1()       // North R1
get_LoadNR2()       // North R2
get_LoadNR3()       // North R3

get_LoadSR1()       // South R1
get_LoadSR2()       // South R2
get_LoadSR3()       // South R3
```

The load-current quantities are provided separately so that downstream
calibrations may select the grouping appropriate to their application.
No current-to-calibration conversion is performed by `TpcConditionsReco`.

## Using TpcConditions in a Fun4All Module

A consumer can obtain the current-event conditions using the standard PHOOL
node-tree lookup:

```cpp
#include <tpcconditions/TpcConditions.h>

#include <phool/getClass.h>

TpcConditions *conditions =
    findNode::getClass<TpcConditions>(topNode, "TpcConditions");
```

The conditions can then be accessed directly:

```cpp
float temperature = conditions->get_Temperature();
float pressure    = conditions->get_Pressure();

bool fieldOK = conditions->get_FieldOK();
bool gainOK  = conditions->get_GainOK();

float load      = conditions->get_LoadCurrent();
float loadNorth = conditions->get_LoadNorth();
float loadSouth = conditions->get_LoadSouth();

float loadNR1 = conditions->get_LoadNR1();
float loadNR2 = conditions->get_LoadNR2();
float loadNR3 = conditions->get_LoadNR3();

float loadSR1 = conditions->get_LoadSR1();
float loadSR2 = conditions->get_LoadSR2();
float loadSR3 = conditions->get_LoadSR3();
```

`TpcConditionsReco` must be registered before modules that consume the
`TpcConditions` node so that the conditions have been updated for the
current event before they are used.

## CDB Payload

The conditions are stored in CDB under

```text
TPC_CONDITIONS
```

Each payload contains a time history for one run.  The individual records
are indexed internally by BCO.  For each event, `TpcConditionsReco` selects
the latest conditions sample whose BCO is less than or equal to the event
BCO.

The payload contains the gas temperature, gas pressure, TPC HV status flags,
and individual GEM G4 currents used to calculate the load-current quantities.