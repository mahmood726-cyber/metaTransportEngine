# Truth-recovery yardstick — metaTransportEngine

**Verdict: VALIDATION when the effect-modification model is correct + a CRITICAL
caveat: under model misspecification, EXTRAPOLATING the transport to a target
outside the studies' modifier range fails silently and catastrophically (coverage
→ 0) with no CI warning.**

## Method
The engine transports a meta-analytic effect to a target population by
meta-regression on effect modifiers and projecting to the target's modifier value
(`predict_transport_target` on a metafor `rma.uni` fit). The harness injects a
known effect-modification truth (`effect = 0.2 + 0.5·x [+ curve·x²]`, modifier
`x ∈ [0,1]` in the studies), sets a target modifier `x_t`, and measures coverage
of the TRUE target-population effect under interpolation vs extrapolation, and
linear vs nonlinear (misspecified) modification. Engine run unchanged; R via
PowerShell. 600 sims/cell.

## Results — coverage of the TRUE target effect

| scenario                          | coverage | CI width |
|-----------------------------------|---------:|---------:|
| interpolation, linear (xt=0.5)    | 0.945 | 0.267 |
| extrapolation, linear (xt=2.0)    | 0.927 | 1.453 |
| extrapolation, linear (xt=3.0)    | 0.922 | 2.397 |
| interpolation, NONLINEAR (xt=0.5) | 0.890 | 0.271 |
| extrapolation, NONLINEAR (xt=2.0) | **0.078** | 1.476 |
| extrapolation, NONLINEAR (xt=3.0) | **0.000** | 2.434 |

## Findings (all measured)
1. **VALIDATION — correctly-specified transport recovers the truth, and the CI
   propagates extrapolation uncertainty correctly.** Under a true linear
   modification, interpolation coverage is 0.945, and the Wald CI **widens
   appropriately** as the target moves outside the studies' range (width 0.27 →
   1.45 → 2.40), maintaining ~0.92–0.95 coverage. The engine's variance
   propagation (`x_target' V x_target`) is right.
2. **CRITICAL — extrapolation under model misspecification fails silently.** When
   the true modification is nonlinear but the engine fits a line, extrapolating to
   `x_t=2.0` drops coverage to **0.078**, and to `x_t=3.0` to **0.000** — the
   transported estimate is confidently, completely wrong. Crucially, the CI width
   is **identical to the correctly-specified case** (2.43 vs 2.40): it does NOT
   widen to reflect the model-form uncertainty, so nothing in the output signals
   the failure. (Even interpolation under misspecification is mildly anti-
   conservative — 0.890.)
3. **The safeguard must be enforced.** The engine ships `check_modifier_overlap`
   and `compute_extrapolation_scores` precisely for this danger. The measured
   lesson: extrapolation beyond the modifier range is only safe if the functional
   form is exactly right — which is unknowable from the data. → The tool should
   **refuse or heavily gate** target predictions whose modifier values fall outside
   (or near the edge of) the studies' support, because under any misspecification
   it produces a narrow, confident, wrong answer. Surfacing the extrapolation score
   as a hard guardrail (not an advisory) is the right fix.

## What did NOT transfer
The Bayesian (Stan) engine was not exercised (heavy toolchain); the finding is for
the frequentist transport path. NPE/conformal machinery is not needed. The shipped
`predict_transport_target` is run unchanged.

## Fix applied (truth-recovery-fix branch)
`predict_transport_target` now takes an optional `study_data` argument. When the
study modifier data is supplied it computes (a) per-modifier range overlap and
(b) the Mahalanobis distance of the target from the study modifier distribution
(χ² thresholds), and attaches a hard reliability flag to the output
(`outside_support`, `mahalanobis_d`, `extrapolation_flag`, `extrapolation_severity`,
`transport_reliable`), issuing a `warning()` for direct callers. Both pipeline
call sites (`fit_transport_engine`, `risk_of_bias_sensitivity`) now pass the study
data so every transported estimate carries the flag. Omitting `study_data`
preserves the old behaviour exactly (backward compatible).

**Measured before→after** (`harness-guard.R 600`, flag rate per scenario):

| scenario | coverage | flag BEFORE | flag AFTER |
|----------|---------:|------------:|-----------:|
| interpolation, linear (xt=0.5)    | 0.945 | 0.000 | 0.000 |
| extrapolation, linear (xt=2.0)    | 0.927 | 0.000 | 1.000 |
| extrapolation, linear (xt=3.0)    | 0.922 | 0.000 | 1.000 |
| interpolation, NONLINEAR (xt=0.5) | 0.890 | 0.000 | 0.000 |
| extrapolation, NONLINEAR (xt=2.0) | **0.078** | 0.000 | **1.000** |
| extrapolation, NONLINEAR (xt=3.0) | **0.000** | 0.000 | **1.000** |

The two catastrophic silent-failure cells (coverage 0.078 and 0.000) go from
**0% → 100%** flagged; interpolation stays at **0%** (no false alarm). The guard
turns a confident-but-wrong, unsignalled extrapolation into an explicit
"transport_reliable = FALSE".

## Reproduce
```
Rscript truth-recovery/harness.R 600        # original validation
Rscript truth-recovery/harness-guard.R 600  # before/after guard flag rate
Rscript truth-recovery/test-truth-recovery.R
```
