# rhouecalib: calibration of the rho-based UE subtraction

Makers and payload builder for the three CDB payloads read by the rho-based underlying-event
subtraction of `offline/packages/jetbackground` (`DetermineTowerRho`, `TowerRhoToBackground`,
`MbdChargeKnee`). The Fun4All macros and the condor chain that run them are in the `macros`
repository, `calibrations/jets/rho_ue_calib/`.

## Payloads

| CDB domain | Content | Read by |
|---|---|---|
| `JET_RHO_MBDQ_KNEE` | the MBD charge knee Q_knee of one run (IOV = that run) | `MbdChargeKnee` (Q' = (q_S + q_N) / Q_knee) |
| `JET_RHO_TOWER_SHAPE` | tower shape map s(class, layer, ieta, iphi) for the seed jets | `DetermineTowerRho::set_tower_shape_cdb_tag` |
| `JET_RHO_UE_CALIB` | per layer and eta strip: w, alpha, beta in classes of vertex z and Q' | `TowerRhoToBackground::set_eta_calib_cdb_tag` |

The UE model is UE = alpha x rho x cosh(eta) + beta (affine) or w x rho x cosh(eta)
(multiplicative), interpolated bilinearly in (z, Q'). All CDBHistos; the contents are listed in
`RhoUECalibPayload.h`. No `_default` payload exists for the knee domain by design: a run without a
knee payload uses the per-run knee table carried by the s-map and UE payloads, then the median
of that table (with a warning).

## Steps

The calibration runs three passes over the same minimum-bias events, each on the merged output
of the one before:

1. `MbdChargeHistMaker`: the MBD charge spectrum and t0 per run;
   `RhoUECalibPayload::FitMbdKnees` fits the knee of every run (knee table), and
   `MakeKneePayloads` writes the per-run knee payloads.
2. `TowerShapeMapMaker`: tower E_T sums per (z, Q') class; `MakeTowerShapePayload` writes the s-map.
3. `RhoUECalibMaker`, after the three `DetermineTowerRho` with that s-map: rho vs strip E_T sums;
   `MakeRhoUEPayload` fits w, alpha and beta.

`RhoUECalibBase` holds the event selection common to the makers (minimum bias, vertex class
range).

The calibration must be made with the same `DetermineTowerRho` settings (seed jet acceptance,
tower shape map) as the subtraction it is applied to.
