# opal-sandbox (development)

- 2026/11/03: South-west Pacific albacore case study (`case-studies/swpo-alb/`).
  - An MFCL bridge from the 2024 WCPFC diagnostic model (two regions, MFCL
    2.2.7.0, now in `executables/`) to model region 1: M0 reproduces the
    supplied objective; M1 sets movement to zero; M2 keeps region 1 only (par
    built from M1's estimates); M3 fixes growth and drops the conditional
    age-at-length data. Region 1's spawning biomass roughly doubles at each
    step (1954: 0.87, 1.7, 3.4, 6.6 Mt; 2022 depletion 0.50, 0.57, 0.62, 0.75).
  - New shared code `shared/mfcl-to-opal/`: base-R readers for MFCL par and
    plot-report files and `mfcl_to_opal()`, which maps MFCL's
    October-September model year (`.frq` month_1), catch, CPUE (time-varying
    CVs), length compositions and biology onto a quarterly opal model.
    `parse_frq()` now reads `.frq` files without weight frequencies.
  - At MFCL's values opal (O0) reproduces MFCL's initial numbers and
    recruitment exactly and its spawning biomass within 2.1%; replays put
    almost all of the gap on catch timing (start-of-quarter removal against
    MFCL's Baranov at each incident), little on seasonal selectivity. Fitted
    with MFCL's selectivity (O1), opal's scale is unbounded: CPUE and the
    recruitment penalty favour MFCL's scale but the annual length
    compositions favour an unexploited stock.
  - The cause is length-at-age timing: MFCL predicts each length sample at
    its date with lengths that grow through the year, opal from start-of-year
    lengths. With opal's growth at MFCL's mid-year lengths (O1-mid) the fit
    is bounded, with B0 at MFCL's Beverton-Holt B0 (log 13.96 against 14.04);
    its stock is still about half MFCL's in 2022, because MFCL scales its
    initial population apart from its stock-recruitment relationship;
    estimating opal's initial population (O1i) brings it to 709 against
    860 kt (depletion 0.72 against 0.75). Fits with opal's own double-normal
    selectivity (O2, O2q, O2-split; three parameters free) collapse to
    knife-edge curves: the longline compositions are narrower than the
    model's spread of length within an age. Estimating the SD of length at the
    oldest age (4-5 cm against MFCL's 7.6), with tail pooling and a
    Dirichlet-multinomial, gives the best fits and a stock near MFCL's.
    Age-based double normals fit the longline compositions better than
    MFCL's selectivity but miss the troll fisheries; a hybrid (age-based
    longline, length-based troll and driftnet; O4) fits best and gives 886 kt
    and depletion 0.76 in 2022 (MFCL 860, 0.75). Seasonal fishery
    splits (`mfcl_split_seasonal_fisheries()`) change little; a quarterly-age
    set-up (`mfcl_quarterly_ages()`: quarterly ages and steps, recruitment
    kept annual) also profiles to MFCL's B0 and fits the indices best.
    The M1 -> M2 transplant difference (up to 1.3% in adult biomass with
    identical inputs) is MFCL's period structure: movement weeks add periods
    only in multi-region models and M is shared equally among periods;
    aligning the movement weeks with fishing dates removes nearly all of it.
    The vignette `vignettes/swpo-alb-mfcl-to-opal.qmd` (HTML and PDF)
    documents the case study and the opal features the stock would need.

- 2026/10/28: North Pacific blue shark case study (`case-studies/npo-bsh/`).
  - An SS3 bridge from the 2022 ISC assessment (S6, two sexes) to a
    single-sex SS3 model, one change per step: single sex, F18 spline to
    double normal, F8's selectivity blocks as two fleets, no F19 selectivity
    deviations, no initial regime offset, no F-ballpark penalty. Every step
    converges with a positive-definite Hessian. Likelihood profiles over
    ln(R0) show that the single-sex models' scale is weakly determined.
  - The last step (B2e) is converted to opal. At SS3's values opal differs
    from SS3 through catch timing (start-of-year removal). Fitted with SS3's
    selectivity (O1), opal converges but estimates terminal depletion of
    0.25 against SS3's 0.54; diagnostic refits trace this to opal predicting
    length compositions from start-of-year lengths (with SS3's mid-year
    length-at-age, or without the compositions, opal matches SS3). Fitting
    opal's own double normals (O2) does not converge. The vignette
    `vignettes/npo-bsh-ss3-to-opal.qmd` documents these and the other places
    where opal falls short of SS3 or differs from it.
  - The shared converter reads SS3 3.30 models: 3.30 par files and variance
    adjustments (including SS3's additive index CV), fishing fleets numbered
    after the surveys, initial F, early recruitment deviations as the initial
    age structure, the female fraction of `Nsexes = -1`, and recruitment
    deviations aligned to SS3's recruitment year when SS3 age 0 is kept.
    `ss3_seed_par_sel_size_dblnorm()` carries SS3's size double normals into
    opal exactly. The IO yellowfin conversion is unchanged.
  - The vignette summarises the assessment report (ISC 2022, Annex 12): the
    stock, fisheries, data, data weighting, biology, the three-model ensemble
    (S6_base has half its weight) and S6_base's diagnostics. It checks the S6
    files against the report, notes where they differ from the report's text
    (recruitment settles in month 7, not at the start of the year; sigma_R
    0.4; the updated growth parameters; constant fecundity, so spawning output
    is in pups; the F-ballpark penalty and initial regime offset), and
    compares the S6 re-run with the ensemble's reference points. The bridge
    table's F2020 column is now labelled F2020/F_MSY, which is what SS3
    reports with this starter file, and the bridge cache adds F_MSY.
  - New tests: `test-npo-bsh-ss3-bridge.r` (each bridge step's convergence,
    Hessian, structure and data) and `test-ss3-to-opal-npo-bsh.r` (the
    conversion, recruitment alignment, selectivity and opal's dynamics at
    SS3's values).

- 2026/10/28 ([#8](https://github.com/N-DucharmeBarth/opal-sandbox/pull/8)):
  Continuous integration.
  - `.github/workflows/tests.yaml` runs `tests/testthat` on every pull
    request and on every push to `main`. It installs opal from `dev`,
    caches the R package library between runs and uploads the tests'
    coverage of `shared/` to Codecov (`tests/coverage.r`, `codecov.yml`).
  - The new test `test-ss3-to-opal-io-yft.r` checks the SS3-to-opal
    converter on the IO yellowfin case study. The structural checks and the
    reconciliation with SS3 pass, the opal object builds, and replaying the
    converted inputs with SS3's catch timing reproduces SS3's spawning
    biomass to within 1%.
  - The SS3 length- and weight-composition tests read `CompReport.sso`,
    which is not in the repository, so they errored in a fresh clone. They
    are now skipped until the file is added.
  - `extract_mfcl_cpue()` now returns `Time` as a double, like
    `extract_ss3_cpue()`. Its rows are now sorted by fishery, year and month
    before `Time` is assigned; the sort had been computed but not kept.
  - The CPUE extractor test no longer overwrites the tracked
    `mfcl/v11/cpue.csv`.
- 2026/10/28 ([#7](https://github.com/N-DucharmeBarth/opal-sandbox/pull/7)):
  The repository now holds opal case studies and tests.
  - Each stock has its own folder under `case-studies/` (`wcpo-bet`,
    `io-yft`). Code used by more than one case study is in `shared/`: the
    SS3-to-opal converter (`ss3-to-opal/`), the SS3 and MFCL extractors
    (`extractors/`) and the cloud install script (`setup/`). The opal
    simulation studies are in `simulations/`.
  - The southwest Pacific striped marlin SS3 base file moved to
    `case-studies/wcpo-bet/`, the case study it was used to build.
  - The README was rewritten for the new purpose. The WCPO bigeye
    background moved to `case-studies/wcpo-bet/README.md`.
- 2026/10/13–2026/10/27: IO yellowfin case study.
  - A YAML-driven SS3-to-opal converter, with validation, reconciliation and
    bridging checks.
  - Vignettes comparing opal with SS3 and with the simulated truth,
    including the age-specific K and age-based selectivity fits, and a
    step-by-step workflow guide.
  - Catch-timing and selectivity-basis simulations.
  - Exports of the data, fits and articles to the opal package
    (`opal-export/`, `opal-patches/`).
- Earlier: the WCPO bigeye single-region MFCL model was rebuilt in SS3, and
  the SS3 and MFCL output extractors were added with unit tests.
