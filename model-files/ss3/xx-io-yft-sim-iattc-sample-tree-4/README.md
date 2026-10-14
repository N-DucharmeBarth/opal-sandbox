# SS3 model: IO yellowfin, simulated data (replicate 4)

The Stock Synthesis (v3.24AB) estimation model used as the baseline for the
SS3-to-opal conversion in `code/rtmb/`.

## Credit

This model was developed by the IATTC analyst team for the Spatial Stock
Assessment Simulation Experiment (Spatial Stock Assessment Methods
Workshop): **Mark N. Maunder, Haikun Xu and Carolina Minte-Vera**
(Inter-American Tropical Tuna Commission). The workshop's team listing gives
the team's platform as SS3*VAST:
<https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop/blob/main/vignettes/Team_info.Rmd>.

The yellowfin assessments from the experiment, including the IATTC team's,
are evaluated in:

> Goethel, D.R., et al. (2024). 'Drivin' with your eyes closed': Results from
> an international, blinded simulation experiment to evaluate spatial stock
> assessments. *Fish and Fisheries* 25(3): 471–490.
> <https://doi.org/10.1111/faf.12819>

Please cite that paper and acknowledge the IATTC team when using this model.

## What it is

- Fitted to data simulated from replicate 4 of the experiment's Indian Ocean
  yellowfin operating model
  (<https://aaronmberger-nwfsc.github.io/Spatial-Assessment-Modeling-Workshop/articles/OM_description_YFT.html>).
  The simulated truth for that replicate is in
  `model-files/sim-truth/io-yft-iattc-sample-tree-4/`.
- Single region, 256 quarterly "years" (`nseas = 1`, `months_per_seas = 3`),
  ages 0–28, 16 fishing fleets and one longline CPUE survey.
- `GrowthModel 3` (age-specific K), age-based selectivity (double-normal for
  the fisheries, logistic for the survey), Beverton-Holt recruitment.
- Run with `SS-V3.24AB-safe; 12/20/2016` (Stock Synthesis by Richard Methot,
  NOAA; ADMB 11.2).

## Files

| File | md5 |
|---|---|
| `control.ss_new` | `fb58229ad06bc48df036042601ac16d1` |
| `data.ss_new` | `5574da47c311fc151d4107b3b2bc4e15` |
| `ss3.par` | `83e50a0065e9dc07c6063bb95fd730a1` |
| `Report.sso` | `50771e4d5d474f3f919e55af1887dea4` |

The original input files (`data.ss`, `control.ss`, named in the `ss_new`
headers as `test_data.ss` and `YFT_IO.ctl`) are not in this folder; the
`ss_new` files echo them as read by SS3.
