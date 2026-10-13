# Simulated truth: IO yellowfin, replicate 4

Operating-model output for replicate 4 of the Indian Ocean yellowfin
operating model of the Spatial Stock Assessment Simulation Experiment, run
for the Spatial Stock Assessment Methods Workshop
(<https://aaronmberger-nwfsc.github.io/Spatial-Assessment-Modeling-Workshop/>;
repository <https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop>).
The SS3 model in `model-files/ss3/xx-io-yft-sim-iattc-sample-tree-4` is fitted
to data simulated from this replicate.

The operating model is a spatial population model (SPM) with 4 regions, ages
1–28 (plus group) and 256 years of 3 time steps each, i.e. one model year per
quarter, aligned one-to-one with SS3 years 1–256.

## Files

Copied unchanged from the replicate's `4/summary_output/`, except
`SSB_all.txt`, which is the SPM report `4/SSB_all` renamed so the repository
tracks it.

| File | Contents | Source md5 |
|---|---|---|
| `ssb_yr_reg_ESS5_run4.csv` | Spawning biomass by year (rows) and region (V1–V4) | `1e68004df42114571eeeb4d348010011` |
| `biomass_yr_reg_ESS5_run4.csv` | Total biomass by year and region | `7f0f911f9171e522b1ebbda9e235d5f7` |
| `vuln_bio_reg_ESS5_run4.csv` | Vulnerable biomass by year and region | `35ff09475fbd4ec4fa37ab38d0fc6622` |
| `recruitment_yr_reg_ESS5_run4.csv` | Recruitment by year and region (V5 is empty) | `f717e750e07a72ddd5231052544313a9` |
| `SSB_all.txt` | SPM `derived_quantity` report for SSB, including both initialisation phases | `72be4fdfecd3235707b6976c88f49c51` |

## Definitions (from the replicate's `config.spm`, md5 `d1efcd78fc880ded3d781caa1424a82c`)

* **SSB** is the biomass of the `mature` category at time step `one`, which
  is the start of the year before fishing (time step `two`). Maturity is a
  process that moves fish between categories, not an ogive. Weight is
  2.459e-8 × L^2.9667 in tonnes, so SSB is in tonnes, the same units as SS3's
  spawning biomass for this model.
* **B0** is the SSB at the end of initialisation phase `phase2`
  (`@process BHrecruitment: b0 phase2`). The initialisation time steps
  contain recruitment, movement, natural mortality, maturation and ageing,
  with no fishing. In `SSB_all.txt` this is the last value of the second
  `initialisation_phase` line: 3.92888e+06 t. Year-1 SSB is 3,928,884 t.
* Total SSB summed over regions equals the `SSB_all` year values to the
  report's six significant figures.

`code/rtmb/scripts/02-extract-io-yft-sim-truth.r` reads these files and
writes the cache the vignette uses.
