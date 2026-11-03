# PHJewel — JEWEL 2.6.0 inside Fun4All

`PHJewel` is a SubsysReco that runs JEWEL 2.6.0 in-process, the way `PHPythia8` wraps PYTHIA 8 and
`PHHerwig7` wraps Herwig 7. The seed comes from `PHRandomSeed()`; events go straight to
`PHHepMCGenEventMap`; the generated count, weight sum and cross section go to `RUN/PHGenIntegral`.
Truth triggers run inside the generation loop, so only triggered events reach the node tree. JEWEL
runs in the same job as Geant4 (one pass).

```
user card ──JewelCard (strict)──> canonical card ──> libjewel_<medium>.so (dlopen)
                                                        │ init / initpythia / genevent  (JEWEL's own code)
                                                        │ HepMC2 lines in memory (JEWEL's own writer logic)
                       PHJewel <──HepMC::IO_GenEvent────┘
                          │ b [fm], PdfInfo, process id added
                          │ PHJewelGenTrigger(s) (regenerate until pass)
                          └──> PHHepMCGenEventMap, RUN/PHGenIntegral ──> DST
```

## Usage

```cpp
#include <phjewel/PHJewel.h>
#include <phjewel/PHJewelJetTrigger.h>
R__LOAD_LIBRARY(libPHJewel.so)

PHJewel *jw = new PHJewel();
jw->set_card("auau200_0010_recoil.dat");  // main card; it names the medium card
jw->set_medium("simple");                 // simple | brick | vac
// optional
jw->set_seed(0);              // 0 (default): PHRandomSeed();  k*1000 reproduces a standalone NJOB k run
jw->set_output_stem("jw");    // jw.log, jw.card.dat, jw.medium.dat  (default PHJewel_<seed>)

PHJewelJetTrigger *trig = new PHJewelJetTrigger();
trig->SetJetR(0.4);
trig->SetMinJetPt(10.);
trig->SetEtaMax(1.0);
jw->register_trigger(trig);   // several: set_trigger_OR() (default) or set_trigger_AND()

Input::ApplysPHENIXBeamParameter(jw);   // vertex and crossing angle, as for the other generators
se->registerSubsystem(jw);
```

`se->run(n)` counts **events written** (passing the trigger), as for PHPythia8. There is no limit
on trigger attempts unless `set_max_trigger_attempts(n)` is called.

Only one PHJewel is allowed per process, because JEWEL keeps its state in Fortran common blocks. Do not
add PHPythia6 to the same process: the PYTHIA 6 symbols are isolated, but the Fortran I/O unit numbers
are shared by the runtime.

Other settings: `set_fortran_hepmc_file()` (JEWEL's own HepMC writer, every JEWEL event, for
validation), `set_evinfo_file()` (EVINFOFILE), `set_isospin_block()` (see Isospin),
`allow_table_generation()`, `set_core_library()`, `save_integrated_luminosity(false)`.

### Truth triggers

| class | passes if |
|---|---|
| `PHJewelJetTrigger` | the leading anti-k_T jet (status-1 particles, no neutrinos) with \|η\| < `SetEtaMax` has p_T ≥ `SetMinJetPt`, and lies inside `SetLeadJetWindow(lo, hi)` if set (p̂_T-bin stitching) |
| `PHJewelPhotonTrigger` | a status-1 photon has p_T ≥ `SetMinPt` and \|η\| ≤ `SetEtaMax`; optionally isolated (`SetIsolation(on, R, frac)`) and/or prompt (`SetPrompt()`: no hadron parent) |

With `KEEPRECOILS T` the jet trigger sees **unsubtracted** jets: the status-3 thermal momenta are
excluded but the recoils are in. Those jets are harder than subtracted ones, so a trigger at X keeps
every event whose subtracted jet is above X.

## Input cards

PHJewel reads **JEWEL's own card format**, so one file drives both PHJewel and standalone JEWEL. It reads
it **strictly**, because JEWEL's reader never fails. An unknown key, an unreadable value, a missing
medium file or a forgotten key all leave a built-in default in place, and those defaults are LHC
Pb+Pb (SQRTS 2760, MASS 208, TI 0.59, SIGMANN 6.2). A RHIC card with one typo runs to completion as a
different sample, and nothing in the output says so.

A run uses a **pair** of cards: the main card, plus the medium card named by `MEDIUMPARAMS` (not needed
for `vac`). The medium type is not in either card. `set_medium()` gives it, and the medium card must
contain exactly that medium's keys. A brick card given to `simple` fails on `TBRICK`.

### Grammar (the subset standalone JEWEL reads the same way)

| rule | why |
|---|---|
| `KEY value`: key in column 1, separated by **one space** | JEWEL splits at the first space; a key is everything before it |
| `#` comments **in column 1 only** | an indented `#` gives JEWEL an empty key and the line is ignored |
| no tabs, no DOS line endings | a tab ends up inside the key; `\r` ends up inside the value |
| key lines ≤ 80 characters (main), ≤ 100 (medium) | standalone JEWEL silently truncates longer lines; use relative paths |
| every key at most once | JEWEL takes the last one |
| integers without a decimal point (`NF 3`, not `3.`) | Fortran list-directed integer read fails on `3.` and keeps the default |
| reals as Fortran literals (`200.`, `1.5`, `1.d-3`) | `1,5` or `20GeV` fail the same way |
| logicals `T` / `F` | anything else fails the same way |
| paths: one token, nothing after it | JEWEL reads a path as the whole rest of the line, comment included |
| a trailing `# comment` is allowed after numbers and logicals | list-directed reads stop at the blank |
| ≤ 1000 lines | JEWEL reads the first 1000 |

**Relative paths are relative to the card file's directory.** Standalone JEWEL resolves them relative
to the working directory, so to run the same card standalone, run from the card's directory.

### Main card keys

**Required.** These have no default in PHJewel, because JEWEL's defaults are for the LHC:

| key | type | meaning |
|---|---|---|
| `PROCESS` | PPJJ PPYJ PPYQ PPYG PPZJ PPZQ PPZG PPWJ PPWQ PPWG PPDY | hard process (EEJJ refused) |
| `SQRTS` | real | √s_NN [GeV] |
| `PTMIN` | real ≥ 3 | p̂_T minimum [GeV]; JEWEL would raise anything lower to 3 without saying so |
| `PTMAX` | real | p̂_T maximum [GeV], `-1.` = none |
| `ETAMAX` | real > 0 | max \|η\| of ME partons and of the medium |
| `PDFSET` | int | LHAPDF ID of the proton PDF (13100 = CT14nlo) |
| `MASS`, `NPROTON` | int | nucleus A and Z: **1 1 for p+p** (JEWEL's default is 208 82) |
| `ISOCHANNEL` | XX PP PN NP NN | XX mixes nucleon pairs by Z/A; otherwise one pair only |
| `WEIGHTED` | T/F | p̂_T-weighted generation; with T, `WEXPO` is also required |
| `KEEPRECOILS` | T/F | keep medium recoils; with T, `WRITESCATCEN T` is required |
| `SPLITINTFILE`, `PDFFILE`, `XSECFILE` | path | JEWEL's tables; must exist (see Tables) |
| `MEDIUMPARAMS` | path | medium card (required for simple and brick, ignored for vac) |

**Optional.** JEWEL's built-in value applies if the key is absent:

| key | default | | key | default |
|---|---|---|---|---|
| `WEXPO` | 5. | | `HADROTYPE` | 0 (0 1) |
| `NF` | 3 (must equal medium NF) | | `PI0DEC` | T |
| `LAMBDAQCD` | 0.4 | | `ISRSCAT` | T |
| `LAMBDAQCDPS` | 0.4 | | `SHORTHEPMC` | T |
| `Q0` | 1.5 | | `COMPRESS` | T |
| `CHANNEL` | MUON (ELEC) | | `WRITESCATCEN` | F |
| `PDFALPHAS` | 0.118 | | `RECSOFTCUT` | 0. |
| `XSECSCALEFAC` | 1. | | `RECHARDCUT` | 5. |
| `MPI` | F | | `KINMODE` | 1 (0-3) |
| `MPIFSR` | T | | `RECMODE` | 0 (0-3; 2 turns recoils off) |
| `ANGORD` | 1 (0-5) | | `COHLENGTHFAC` | 1. |
| `AOFIRSTSPLIT` | F | | `PCOHREJ` | 1. |
| `SCATRECOIL` | F | | `COHSCAT` | T |
| `HADRO` | T (must be T) | | `RSUB`, `PUTONSHELL` | 1., 0 (no effect without DOSUBTRACTION) |

**Owned by PHJewel.** These are accepted but ignored, with a note in the log:

| key | replaced by |
|---|---|
| `NEVENT` | `se->run(n)` |
| `NJOB` | `set_seed()` / `PHRandomSeed()`; `MRPY(1) = seed`, so NJOB k ≡ seed 1000k |
| `LOGFILE` | `<stem>.log` |
| `HEPMCFILE` | the node tree (`set_fortran_hepmc_file` for validation) |
| `EVINFOFILE` | `set_evinfo_file` |
| `READRAN` | never |

**Refused values.** These are legal for JEWEL but wrong in front of a detector simulation:

| setting | why |
|---|---|
| `WRITEDUMMIES T` | adds one **status-1** π⁰ placeholder per thermal parton (E = 10⁻⁶ GeV along its direction, a ghost for constituent subtraction). There are hundreds per event; they carry no visible energy, but Geant and every truth-particle count (multiplicities, constituents, fragmentation) treat them as real particles. Recoil subtraction needs only `WRITESCATCEN T`. |
| `DOSUBTRACTION T` | rewrites the final state inside JEWEL; a DST needs the unsubtracted event |
| `KEEPRECOILS T` with `WRITESCATCEN F` | the recoil contribution could never be subtracted |
| `HADRO F` | partons in the final state |
| `PROCESS EEJJ` | not a hadron collider |
| `RECMODE 2` with `KEEPRECOILS T` | JEWEL silently turns recoils off |

### Medium cards

`simple` (Bjorken-expanding Glauber medium). Units: TAUI fm/c, temperatures GeV, **SIGMANN fm²
(42 mb = 4.2)**, centrality %.

| key | required | default | check |
|---|---|---|---|
| `TAUI`, `TI`, `TC` | yes | | TI > TC > 0; TI > 2 refused as "MeV" |
| `CENTRMIN`, `CENTRMAX` | yes | | 0 ≤ min < max ≤ 100 |
| `SIGMANN` | yes | | > 15 refused as "mb" |
| `NF` | yes | | = main-card NF |
| `MDFACTOR`, `MDSCALEFAC` | yes | | > 0 |
| `WOODSSAXON` | no | T | |
| `N0`, `D` | no | 0.17, 0.54 | > 0 |

Also refused: MASS < 2 with the simple medium.

`brick` (static uniform brick): `TBRICK` (GeV), `LBRICK` (fm), `NF`, `MDFACTOR`, `MDSCALEFAC`, all required.

`vac`: no medium card.

### What JEWEL is actually given

PHJewel writes `<stem>.card.dat` and `<stem>.medium.dat`, then hands those to JEWEL:
- keys that PHJewel owns are dropped, and paths are made absolute;
- every key you did not write appears as a `#default KEY value` comment, so the file is the complete record of the run;
- defaults are **not** written as active keys. JEWEL's built-in defaults are single-precision literals (`lqcd = 0.4` is 0.4000000059), so `LAMBDAQCD 0.4` would change the run in the 8th digit and break event identity with standalone JEWEL.

A 64-bit hash of the canonical pair goes into the log and into `PHGenIntegral::get_Description()`,
together with the medium, card path, seed and embedding id. After `init` PHJewel reads JEWEL's own
values back (A, Z, NF, √s, p̂_T range, η range, Q0, Λ, the recoil, subtraction and hadronisation
switches, process and isospin). It stops if any differ from its own reading of the card.

### Tables

`SPLITINTFILE`, `PDFFILE` and `XSECFILE` must exist. If one is missing, JEWEL integrates it, which is
slow, and writes it to the path, which is a race when parallel jobs share a directory.
`allow_table_generation()` permits that for a dedicated one-off run.

JEWEL checks only that a table exists, not what it was built with. The splitting integrals depend on
Q0 and the energy reach. The PDF and cross-section tables also depend on the medium's Debye-mass
range. Lookups outside a table extrapolate, and JEWEL counts them. **PHJewel prints the fractions at
End.** With `pdfs_ct14nlo.dat` and `xsecs_ct14nlo.dat`, the RHIC 0-10 % simple medium extrapolates in
about 11-12 % of PDF and cross-section lookups; the brick and vac media extrapolate in none.

---

## The event record

| | |
|---|---|
| status 1 | final state (what Geant transports and the trigger sees) |
| status 3 | thermal momenta of the scattering centres (`WRITESCATCEN T`): subtract them from jets at analysis level |
| status 4 | beam particles; also non-final partons |
| event number | JEWEL's NGOOD, which matches the EVINFOFILE line and the standalone HepMC file; after a trigger the numbers have gaps |
| weight | `weights()[0]` = PYTHIA `PARI(10)` |
| `GenCrossSection` | `PARI(1)` [pb] at that event, a running estimate for the current isospin channel |
| `HeavyIon::impact_parameter` | **b in fm** (simple), LBRICK (brick), -1 (vac); event-plane angle 0 because JEWEL lays b along +x. JEWEL itself writes the centrality in this slot, and PHJewel replaces it |
| `PdfInfo` | id1, id2 (`MSTI(15,16)`), x1, x2 (`PARI(33,34)`), Q = √`PARI(22)`; JEWEL leaves it empty |
| `signal_process_id` | PYTHIA `MSTI(1)` |
| vertices | all at the origin, as JEWEL writes them; the collision vertex comes from PHHepMCGenHelper's settings |

### Normalisation

`RUN/PHGenIntegral` holds `N_Generator_Accepted_Event` (events JEWEL wrote, before the trigger),
`N_Processed_Event` (events in the DST) and `Sum_Of_Weight` (Σw over the JEWEL events).
`Integrated_Lumi` = N/σ [pb⁻¹] follows the PHPythia8 convention. With `WEIGHTED T` a histogram filled
with w is normalised by **σ / Sum_Of_Weight**, where σ = N/L, as in `hwfilter`.

### Isospin (ISOCHANNEL XX)

Standalone JEWEL draws all four channel counts up front from NEVENT. A trigger-driven run does not
know how many events it will generate. PHJewel draws the channel per block of `set_isospin_block(n)`
JEWEL events (default 1000) with probabilities Z², Z(A−Z), (A−Z)Z, (A−Z)² over A². The mix is exact in
expectation, and PYTHIA is re-initialised only when the channel changes. σ is PARI(1) at the end of
each PYTHIA initialisation, weighted by the events made under it. For a fixed ISOCHANNEL the call
sequence is the standalone one, which is why seed 1000k reproduces NJOB k event for event.
### Recoils and 4MomSub

With `KEEPRECOILS T` the medium partons a jet scatters off are kept and hadronised with it (status 1).
Their momentum includes the thermal momentum they had before the scattering, which in data belongs to
the underlying event. `WRITESCATCEN T` writes those thermal momenta as status-3 entries. 4MomSub
(arXiv:1707.04142): cluster the status-3 entries as ghosts together with the final state, then subtract
the full four-momenta of the ghosts in each jet from that jet. Geant never sees the status-3 entries.

## Files

| file | role |
|---|---|
| `PHJewel.{h,cc}` | the SubsysReco |
| `PHJewelGenTrigger.h`, `PHJewelJetTrigger.{h,cc}`, `PHJewelPhotonTrigger.{h,cc}` | truth triggers |
| `JewelCard.{h,cc}` | strict card reader, canonical cards, hash |
| `JewelCore.{h,cc}` | dlopen loader, event loop, isospin channels, bookkeeping (no Fun4All, no HepMC) |
| `JewelRun.{h,cc}` | validate → canonical card → init → read-back check |
| `jewel/jewel_capi.cc`, `jewel/jewel-api.f`, `jewel/exports.map` | the core's C API, its Fortran side, and the list of exported symbols (`jw_*` only) |
| `jewel/jewel-lib.f`, `jewel/medium-{simple,brick,vac}.f`, `pythia6425mod-lhapdf6.f`, `meix.f`, `subtraction.C` | JEWEL 2.6.0 and its PYTHIA 6.4.25 |

Only `PHJewel.h` and the trigger headers are installed. The build makes `libPHJewel.so` and one core
per medium, `libjewel_simple.so`, `libjewel_brick.so` and `libjewel_vac.so`. Each core carries its own
private JEWEL and PYTHIA 6 and exports only `jw_*`; PHJewel `dlopen`s the one named by `set_medium()`
with `RTLD_DEEPBIND`, so it binds to LHAPDF 6 even when the Geant4 decayer has already loaded
PYTHIA 6 and LHAPDF 5.

`jewel/jewel-lib.f` is `jewel-2.6.0.f` without `PROGRAM JEWEL`, plus the API hooks. They are inert
unless `jw_init` sets the API mode, so the standalone program built from the same file is unchanged:
- the parameter file is read from a common block instead of `getarg`;
- PHJewel owns the output files, NJOB and the seed;
- `out.ran` is not written;
- each `WRITE` in `CONVERTTOHEPMC` goes through a helper that writes the same line and also passes the
  numbers to the event buffer;
- file-name buffers are 512 characters instead of 80.
