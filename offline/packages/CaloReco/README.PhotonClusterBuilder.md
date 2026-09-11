# PhotonClusterBuilder BDT loading

BDT scoring is optional and disabled by default. When enabled, the builder
loads a TMVA `RBDT` named `myBDT` during `InitRun` and stores the score as the
`bdt_score` shower-shape parameter. The caller supplies the ordered features
expected by the model.

## Load a model from CDB

Set `set_bdt_model_cdb_domain(...)` before registering the builder. The lookup
uses the standard `CDBInterface` and the job's existing `CDB_GLOBALTAG` and
`TIMESTAMP`. The builder does not set or change either flag.

For the PPG12 baseV3E photon-ID model, configure the builder as follows in a
Fun4All macro after the input clusters, calibrated towers and geometry have
been configured:

```cpp
#include <TMVA/RBDT.hxx>
#include <caloreco/PhotonClusterBuilder.h>
#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libcalo_reco.so)

auto* photons = new PhotonClusterBuilder();
photons->set_do_bdt(true);
photons->set_bdt_model_cdb_domain("PPG12_PHOTON_ID_BDT");
photons->set_bdt_feature_list({
    "cluster_Et", "cluster_weta_cogx", "cluster_wphi_cogx",
    "vertexz", "cluster_Eta", "e11_over_e33",
    "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4",
    "e32_over_e35"});
photons->Verbosity(1);  // Print the resolved model filename.
Fun4AllServer::instance()->registerSubsystem(photons);
```

The published baseV3E payload is available for these settings:

| Sample | `CDB_GLOBALTAG` | `TIMESTAMP` |
| --- | --- | --- |
| DATA | `newcdbtag` | Actual run number, 47289 through 53999 |
| PPG12 MC | `MDC2` | 28 |

These are job-level reconstruction settings, not an extra model-specific
global tag. In the job's normal CDB setup, for example:

```cpp
#include <phool/recoConsts.h>

auto* rc = recoConsts::instance();
rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
rc->set_uint64Flag("TIMESTAMP", 47289);  // Use the actual DATA run.
// PPG12 MC instead uses "MDC2" and 28.
```

Keep the 11 features in exactly the order shown. Ratios with a non-positive
denominator resolve to zero. The model is the baseV3E classifier used in the
PPG12 `8 <= ET < 35 GeV` range; the separate `base_E` classifier outside that
range is not selected automatically. Working-point cuts, preselection and
photon-selection logic remain the analysis's responsibility.

## Use a local file

Existing file-based configurations continue to work. To override a CDB
selection, call this after configuring the CDB domain, before `InitRun`:

```cpp
photons->set_bdt_model_file("/path/to/model_base_v3E_single_tmva.root");
```

`set_bdt_model_file` clears the CDB selection. A subsequent nonempty
`set_bdt_model_cdb_domain` selects CDB again. Without a CDB domain, the legacy
default filename remains `myBDT_5.root` unless explicitly changed.

An empty feature list, empty resolved path, or model-loading exception aborts
run initialization with a diagnostic. CDB mode never falls back to the local
filename. It follows the standard `CDBInterface` behavior, including any
enabled `_default` domain lookup or calibration-file cache; the builder does
not change these job-wide settings. Nonempty feature lists still must match the model's
input count, order and definitions; this interface does not infer them from
the CDB name or change the existing feature evaluation.
