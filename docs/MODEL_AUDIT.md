# Model audit: enzyme model vs. the PNAS paper, and the dimensionless mode

This note records (a) what the code's model actually is, (b) how it differs
from the PNAS paper the researcher was working from, (c) the specific artifacts
found in the original results, and (d) what the new
`paper_inspired_dimensionless` mode changes to address them.

## 1. What the code models

The code is an **enzyme reaction-cycle mass-action model** integrated as ODEs
(`system_ode.m`) and simulated stochastically with the Gillespie SSA
(`run_gillespie_ssa.m`). Species are E, ET, ETS, S, P. The observable is the
steady-state output ratio `P / (P + S)`.

The experiment:

1. Sweep total enzyme E and find `E*`, the value where the deterministic output
   ratio is at the midpoint of its range.
2. Define two nearby input states `E0 = E0_frac·E*` and `E1 = E1_frac·E*`.
3. Run SSA at E0 and E1, build histograms of the output, and measure how
   distinguishable the two output distributions are (KL / symmetric KL / JS).
4. Plot distinguishability against an energy-dissipation quantity (`phi`) and a
   chemical-driving quantity (`delta_mu`).

## 2. How this differs from the PNAS paper

- **Different model.** The paper analyzes a **gene-locus / transcriptional
  Markov chain** and derives an **analytic information rate / decision metric**.
  This code is a continuous enzyme kinetics model and measures the
  **distinguishability of simulated output distributions**. We deliberately do
  **not** reimplement the paper's Markov model or its analytic derivation.
- **Different metric.** Paper: information rate / optimal-decision bounds. Code:
  histogram KL / JS between SSA outputs. These are related in spirit (both are
  "how well can you tell two conditions apart") but are not the same quantity.
- **Different input.** Paper's "input" is a signal driving the locus. Here the
  input is **total enzyme concentration** (E0/E1). If the researcher decides the
  intended biological input is substrate S or influx `alpha`, the E0/E1
  definition must change before final interpretation.

The `paper_inspired_dimensionless` mode keeps the enzyme model but borrows the
paper's *conventions* (dimensionless units, symmetric parameter ranges, a small
±5% input separation, and bounded metrics) so the output is scientifically
defensible and comparable in spirit.

## 3. Artifacts found in the original results

1. **KL saturation at ~34.** Raw histogram KL adds a smoothing floor
   `epsilon = 1e-15` to every bin. When the E0 and E1 histograms do not overlap,
   KL approaches `-log(epsilon) ≈ 34.5`. The ubiquitous "~34" plateau was this
   numerical ceiling, **not** maximal biological distinguishability.
2. **Asymmetric, high-driving parameter ranges.** The original log-uniform
   ranges (e.g. forward rates 10^3..10^7 against reverse 10^-1..10^3) pushed
   almost every sample into a strongly driven, far-from-equilibrium regime, so
   the near-equilibrium part of the distinguishability-vs-energy curve was never
   sampled.
3. **Dimensionally inconsistent `delta_mu`.** The driving was not consistently
   reduced to a dimensionless argument before taking logs.
4. **R·T-scaled `phi` not comparable to kBT.** `phi` was scaled by
   `R·T ≈ 2494 J/mol`, so values were not directly interpretable as kBT and not
   comparable to the paper's dimensionless quantities.
5. **Complex / negative `phi` rows.** Small negative ODE concentrations produced
   complex logs or negative entropy-production terms.
6. **Cost.** The 10k × SSA brute force with `t_end_stoch = 1500` and `nBins = 40`
   was extremely expensive (≈48 h) and still produced mostly saturated KL.

## 4. What `paper_inspired_dimensionless` changes

| Aspect | Legacy | Dimensionless |
|--------|--------|---------------|
| ATP | 1e-3 | **1** |
| R, T (so R·T) | 8.314, 300 (≈2494) | **1, 1 (R·T = 1)** |
| `phi`, `delta_mu` units | J/mol-scaled | **kBT-like** (`phi_kBT_*`, `delta_mu_kBT_*`) |
| Rate ranges (log10) | asymmetric per rate | **symmetric 10^-2..10^2 for all six** |
| Enzyme sweep range | 1e-8..5e-6 M | **1..500** |
| Vol_factor (SSA counts) | 1e9 | **1** |
| E0 / E1 | 0.9·E* / 1.1·E* | **0.95·E* / 1.05·E*** |
| Headline metric | symmetric KL (saturates) | **JS divergence (bounded 0..ln2)** |

What it does **not** change: the enzyme ODE/SSA model itself, the validity
flagging, or the staged pipeline. Raw KL columns are retained for continuity,
but the headline figure is JS divergence and a `log1p(symmetric_KL)` plot is
provided to compress the saturation tail.

### Why these defaults, and what may still need tuning

- The symmetric rate range and ±5% input separation are the main levers that
  remove the saturation artifact: with overlapping E0/E1 outputs, KL and JS no
  longer pin to their ceilings, and JS becomes an honest 0..ln2 measure.
- `alpha = beta = 1`, `ATP = 1`, and the enzyme sweep `1..500` are reasonable
  dimensionless defaults, **not** a calibrated cell model. They make the
  deterministic ratio sweep span 0→1 so that `E*` is well defined, and keep SSA
  copy numbers small and fast. If the researcher wants specific biological
  numbers, these are the knobs to revisit (in `model_constants.m`).
- Validity filtering still excludes ODE failures, complex/negative `phi`,
  `max_steps_hit`, sparse bins, degenerate histograms, and NaN/Inf rows.

## 5. Sanity-audit findings (dimensionless first-result, n=73 valid)

A focused audit of `results/dimensionless_first_result/` (auto-reported in
`summary.md`) found the following. The KL≈34 epsilon artifact is gone, but three
real issues remain.

### 5.1 `phi` and `delta_mu` are NOT a copy bug (resolved)

They look almost identical per row, but that is the correct thermodynamic
identity, not duplication:

```
phi = sum_j (Jf_j - Jb_j)*RT*log(Jf_j/Jb_j) = Jnet * RT*log(forward/backward) = Jnet * delta_mu
```

Verified numerically: max relative error of `phi` vs `Jnet*delta_mu` = 2.9e-3,
**median 7.8e-13**. At steady state the net cycle flux equals the substrate
influx, `Jnet ~ alpha = 1` (95% of valid rows have |Jnet-1|<0.05), so
`phi ~ delta_mu` numerically. **Consequence:** the JS-vs-phi and JS-vs-delta_mu
plots are nearly the same plot, and `phi`'s spread across parameter sets is
driven almost entirely by the affinity `delta_mu`, not by flux. This is fine but
should be stated; do not present them as two independent axes.

### 5.2 Stochastic output is boundary-dominated (caveat)

Median `output_mean` ≈ 0.22–0.24, but median `fraction_at_0` ≈ 0.43 and
**55/73 rows are boundary-dominated** (`fraction_at_0 + fraction_at_1 > 0.5`).
The SSA output `P/(P+S)` piles up at 0 and 1 rather than forming a smooth
interior distribution. Cause: `Vol_factor = 1` gives small molecule counts
(E* often only a few to a few-hundred molecules) and `t_end_stoch = 25` is
short, so individual trajectories are near-deterministic and land on a boundary.
This makes histogram KL fragile (hence the 17 sparse / 13 degenerate exclusions)
and keeps JS small. Levers: raise `Vol_factor` (more molecules, smoother
interior), raise `N_stoch`, or lengthen `t_end_stoch` — after a convergence
check, not blindly.

### 5.3 No strong energy-discrimination trend (key result)

Spearman rank correlations on valid rows:

- `average_phi` vs `JS_divergence`: **rho = -0.07**
- `average_phi` vs `log1p_symmetric_KL`: **rho = -0.19**

**Under the current settings the corrected enzyme model does not show a strong
monotonic "more dissipation -> more distinguishable" trend.** This is an honest
negative result, not a plotting problem. The dominant cause is 5.4.

### 5.4 The E* operating point is poorly resolved (action recommended)

`fig_example_response_curves.png` shows that `P/(P+S)` vs total enzyme is **not
a broad sigmoid across the swept range**. It rises almost vertically over the
first few enzyme units and then saturates flat from roughly E≈20 out to E=500.
Because `E_min=1, E_max=500` is sampled with a *linear* 100-point grid, only a
handful of points land on the steep rising region, and `E*` (the ratio midpoint)
sits at very small enzyme. For low-dynamic-range parameter sets the curve is
flat at ≈0 and `E*` is effectively arbitrary. With `E*` tiny, `E0/E1 = 0.95/1.05*E*`
are a few molecules apart, so the SSA outputs barely differ -> small JS.

**Recommended fix (not yet applied; outside this audit's "no big runs" scope):**

- Sample the enzyme sweep on a **log grid** (e.g. `logspace`) so the steep rising
  region E∈[~0.5, ~30] is well resolved, and/or narrow `E_max`.
- Optionally widen the input separation (e.g. E0/E1 = 0.8/1.2*E*) so the two
  conditions are actually distinguishable.
- Re-check `fig_example_response_curves.png` after either change to confirm `E*`
  lands mid-way on a genuinely monotonic stretch.

### 5.5 Input-signal semantics (caveat, by design)

`E0/E1` vary **total enzyme**, not an upstream signal concentration `c` as in the
paper. A runtime warning is now emitted and a `substrate_or_alpha` preset is
flagged as a TODO in `model_constants.m`. Decide the intended input before any
final biological interpretation.

### 5.6 Linear vs log enzyme sweep (comparison result)

The audit recommendation in 5.4 was tested directly with a log-spaced sweep
(`scripts/run_dimensionless_logspacing_first_result_pipeline.m`,
`enzyme_spacing='log'`). Same bounded settings, separate output dirs.

| Quantity | Linear sweep | Log sweep |
|----------|--------------|-----------|
| `E_star_at_sweep_boundary` rows | (several) | **0** |
| boundary-dominated SSA rows | 55/73 | 44/67 |
| Spearman average_phi vs JS | -0.07 | -0.17 |
| Spearman average_phi vs log1p(symKL) | -0.19 | -0.42 |

**What improved:** log spacing fixes the E* resolution problem. The response
curves (`figures/dimensionless_logspacing_first_result/fig_example_response_curves.png`,
semilogx) show E* sitting mid-transition on genuinely sigmoidal parameter sets
(e.g. runs 205, 361, 406), and `E_star_at_sweep_boundary` drops to 0. Boundary
pile-up also eases somewhat.

**What did NOT improve:** the JS-vs-phi trend is still weak/flat (and what
correlation exists is mildly *negative*, not the expected "more dissipation ->
more distinguishable"). So the flat trend is **not primarily a sweep-resolution
artifact**. The dominant limitation is the experiment design itself:

- the input is a ±5% change in *total enzyme*, not an upstream signal `c`
  (see 5.5), so the two conditions are intrinsically hard to separate; and
- `phi ~ Jnet * delta_mu` with `Jnet` pinned at alpha (5.1), so "energy
  dissipation" mostly tracks the affinity rather than an independent control.

**Recommendation:** adopt the **log-spaced** pipeline as the researcher-facing
dimensionless output (cleaner E* definition, parameter-set responses are
genuine sigmoids). But report honestly that, even with correct E* resolution,
the corrected enzyme model does not reproduce a strong energy-discrimination
trend. The next scientific lever is the **input definition** (substrate/alpha
sweep, 5.5) and/or a larger input separation — not further sweep tuning.

### 5.7 Focused diagnostic deep-dive (decoupling, raw P/S, Vol_factor)

A no-large-run diagnostic (`src/diagnostics/diagnose_dimensionless.m` +
`scripts/run_dimensionless_volfactor_sensitivity.m`, outputs in
`results/diagnostics/` and `figures/diagnostics/`) pinned down *which* of the
candidate causes actually drives the flat trend.

1. **JS is decoupled from the deterministic response geometry**
   (`fig_response_curves_by_JS_extremes.png`). Among valid rows, some of the
   *highest*-JS rows have **flat** P/(P+S) responses, while the *lowest*-JS rows
   include a **textbook sharp sigmoid** with E0/E1 straddling the transition.
   So stochastic distinguishability is **not** driven by input-output response
   geometry or E* placement — this rules out "poorly resolved E*" (5.4) as the
   dominant cause and confirms why fixing the sweep (5.6) did not help.

2. **The output is genuinely low-copy** (`fig_raw_P_S_distributions.png`). The
   *raw* P and S counts are single digits (P ∈ {0..4}, S ∈ {0..1} for the
   low-JS example). So `P/(P+S)` is not manufacturing the boundary pile-up out of
   otherwise-smooth P/S — the counts themselves are too discrete to form an
   interior distribution.

3. **Vol_factor confirms boundary pile-up is a discreteness artifact.** Bounded
   sweep (`results/dimensionless_volfactor_smoke/volfactor_summary.md`):

   | Vol_factor | valid_for_kl | median boundary frac | median output_var |
   |-----------:|-------------:|---------------------:|------------------:|
   | 1   | 6 | **0.74** | 0.073 |
   | 10  | 7 | **0.14** | 0.0005 |
   | 100 | 7 | **0.00** | 0.0001 |

   Raising molecule counts collapses the 0/1 pile-up, so the boundary-heaviness
   at `Vol_factor=1` (5.2) is **low-copy-number discreteness, not true model
   behavior**. (Output variance shrinks because the ratio concentrates on its
   deterministic value as counts grow.) The per-Vol Spearman phi-vs-JS is on only
   ~7 rows and is **not statistically reliable**; a real trend must be confirmed
   on a full run at higher `Vol_factor` (with `max_steps` raised accordingly,
   since SSA cost grows with molecule count).

**Verdict on the four candidate causes:** the flat trend is *partly* a
low-copy-number/`Vol_factor` artifact (boundary pile-up, fixable) and *primarily*
an experiment-design issue — the input (±5% total enzyme with substrate influx
`alpha` fixed) does not make energy dissipation an independent control of
discrimination, and JS is decoupled from the response geometry. It is **not** a
bad-ratio-metric problem per se and **not** an E*-resolution problem.

## 6. Reproducing

```
matlab -batch "run('scripts/run_dimensionless_smoke_pipeline.m')"                      % minutes
matlab -batch "run('scripts/run_dimensionless_first_result_pipeline.m')"               % ~1-2 min
matlab -batch "run('scripts/run_dimensionless_logspacing_first_result_pipeline.m')"    % ~1-2 min
```

Focused diagnostics (no large runs; depend on the log-spaced result above):

```
matlab -batch "run('scripts/run_dimensionless_volfactor_sensitivity.m')"               % <1 min
matlab -batch "addpath('scripts'); setup_repo_paths(); diagnose_dimensionless('ssaFile','results/dimensionless_logspacing_first_result/results_ssa_kl.csv','enzyme_spacing','log')"  % <1 min
```

The diagnostics write `results/diagnostics/{six_run_table.csv,DIAGNOSTIC_REPORT.md}`
and `figures/diagnostics/{fig_response_curves_by_JS_extremes.png,fig_raw_P_S_distributions.png}`.

Outputs land under `results/dimensionless_*` and `figures/dimensionless_*`.
The legacy mode remains available via `scripts/run_smoke_pipeline.m` and
`scripts/run_first_result_pipeline.m`.

## 7. Current conclusion

- The original "KL saturates at ~34" finding was a **numerical artifact** caused
  by histogram KL with epsilon = 1e-15 smoothing (the ceiling is `-log(epsilon)`),
  compounded by asymmetric high-driving rate ranges and RT-scaled units.
- The corrected dimensionless runs **remove that artifact**: KL no longer pins at
  34, JS divergence is bounded and honest, phi spans near-equilibrium to strongly
  driven, units are kBT-like, and 0 trajectories hit the step cap.
- **However, the enzyme model does not show a strong monotonic, paper-like
  energy/discrimination trend** under the current setup. Spearman correlations of
  `average_phi` vs JS / log1p(symmetric_KL) are weak and mildly *negative*
  (|rho| <= ~0.4), in both linear and log sweeps.
- **Testing a log-spaced enzyme sweep did not resolve this.** Log spacing fixed
  the E* resolution problem (E* now sits mid-sigmoid, `E_star_at_sweep_boundary`
  drops to 0) but the JS-vs-phi trend stayed weak/flat.
- **Focused diagnostics (5.7) localized the cause.** JS is decoupled from the
  deterministic response geometry (so E* placement is not the issue), the raw P/S
  counts are single-digit, and a Vol_factor sweep shows the 0/1 boundary pile-up
  is a low-copy-number artifact (median boundary fraction 0.74 -> 0.00 as
  Vol_factor 1 -> 100). The flat trend is therefore an experiment-design issue
  (input definition), not a metric or sweep-resolution bug.
- **Therefore the remaining mismatch is most likely conceptual, not numerical:**
  the choice of model (enzyme mass-action vs the paper's transcriptional Markov
  chain), the definition of the input signal (total enzyme vs substrate/alpha vs
  an upstream concentration `c`), and/or the choice of output metric
  (SSA output-distribution KL/JS vs the paper's information-rate/decision metric).
  These are decisions for the researcher/supervisor, not further parameter tuning.

## 8. What needs researcher/supervisor confirmation

Before any further modeling work, please confirm:

1. **Input signal.** Should the input signal be total enzyme concentration,
   substrate concentration, substrate influx alpha, or something else?
2. **Model intent.** Is the enzyme mass-action model intended as an analogy, or
   should the actual PNAS transcriptional Markov model be implemented?
3. **Output metric.** Is the desired output metric KL/JS between stochastic
   product distributions, or the paper's information-rate/decision metric?
4. **Input separation.** Should E0/E1 be +/-5% around the midpoint, matching
   delta_c = 0.1*c*, or something else?

The current code answers (1) "total enzyme", (2) "analogy", (3) "KL/JS of SSA
output distributions", (4) "+/-5%". If any of these differ from the intended
study, that — not the numerics — is the reason the paper-like trend is absent.
