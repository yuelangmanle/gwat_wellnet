# Supporting Information for "How sparse is too sparse? A closed-form law for the reliability of groundwater-level trends"

## S1. Synthetic-population testbed for the finite-population sampling law

To test whether the sampling law in Equation 1 of the main text holds under known population parameters—and to separate distributional skewness and interannual water-level noise from the composite-hydrograph Theil–Sen estimator—synthetic monitoring networks were generated and thinned under controlled conditions. Each synthetic network consists of N = 1,000 wells observed annually over T = 50 years. Individual well trends b_i were drawn with population standard deviation S = 1 and population mean β_true = SNR · S for SNR ∈ {0.03, 0.10, 0.25}, across three slope distributions:

1. **Gaussian:** b_i ~ Normal(β_true, S²).
2. **Skewed (t₃):** heavy-tailed Student-t with 3 degrees of freedom, scaled to unit variance.
3. **Bimodal:** an equal mixture of two Gaussian modes at β_true ± 0.85 S (within-mode standard deviation 0.53 S, total variance S²), representing an aquifer with both depleting and recovering sub-basins.

Each distribution was evaluated in two observation modes:

- **Noiseless mode:** each replicate draws n wells without replacement (f ∈ {0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90}; 2,000 Monte Carlo replicates per cell) and takes the sample mean of the n slopes, b̄_n. A sign flip is recorded when sign(b̄_n) ≠ sign(β_true).
- **Annual-anomaly mode:** each well is assigned a 50-year time series y_{i,t} = b_i(t − t̄) + ε_{i,t} with independent interannual noise ε_{i,t} ~ Normal(0, σ_e²) (σ_e = 5, comparable to empirical residual variability). For each replicate, annual anomalies are averaged across the n retained wells and the linear slope of the composite series is estimated.

For each cell, the apparent scale factor c_app was fitted by binomial maximum likelihood to P(flip) = Φ[−c_app · SNR · √(nN/(N − n))] across informative retention levels, and the root-mean-square error (RMSE) of the predicted flip probability was evaluated against the theoretical curve with c = 1.

| Slope distribution | Observation mode | SNR = 0.03 c_app (RMSE) | SNR = 0.10 c_app (RMSE) | SNR = 0.25 c_app (RMSE) |
| :--- | :--- | :---: | :---: | :---: |
| Gaussian | Noiseless (b̄_n) | 1.03 (0.012) | 1.02 (0.008) | 1.01 (0.002) |
| Skewed (t₃) | Noiseless (b̄_n) | 1.04 (0.014) | 1.02 (0.010) | 1.03 (0.004) |
| Bimodal | Noiseless (b̄_n) | 1.02 (0.009) | 1.01 (0.006) | 1.01 (0.002) |
| Gaussian | Annual noise (σ_e = 5) | 0.95 (0.015) | 0.93 (0.017) | 0.91 (0.016) |
| Bimodal | Annual noise (σ_e = 5) | 0.94 (0.013) | 0.92 (0.017) | 0.91 (0.014) |

Two conclusions follow from Table S1:

1. **In the noiseless regime, c_app = 1.01–1.04 across all three slope distributions, with flip-probability RMSE ≤ 0.014.** The finite-population factor √(nN/(N − n)) and the normal tail Φ(−x) accurately describe sign-reversal probability even when individual well slopes are heavy-tailed (t₃) or bimodal.
2. **Adding independent interannual water-level noise slightly reduces c_app to 0.91–0.95**, because residual temporal noise adds variance σ_e² / (n ∑(t − t̄)²) to the spatial slope variance S²(N − n)/(nN). Consequently, the empirical value c = 1.27 found for the U.S. Geological Survey (USGS) archive is not an artifact of slope skewness or random annual noise; as shown in Section S5, it arises from the statistical efficiency of applying the non-parametric Theil–Sen estimator to the spatially averaged anomaly hydrograph.

## S2. Zero-parameter forward check of spatially clustered thinning

In the national (CONUS) clustered-thinning experiment, wells were grouped into 0.5° × 0.5° latitude–longitude blocks (retaining blocks with at least 10 wells; B = 621 blocks containing N = 25,754 wells), and whole blocks were removed at random. To verify that the elevated sign-reversal rate under clustered removal follows directly from survey-sampling theory without refitting any parameter, the design effect (DEFF; Kish 1965) was computed from the archive's block-level slope statistics:

- Mean block size: m̄ = 41.5 wells per block; block-size coefficient of variation: cv(m) = 1.04.
- One-way random-effects analysis of variance on individual well slopes across blocks gives between-block variance σ_b² = 0.0725 (ft yr⁻¹)² and within-block variance σ_w² = 0.1563 (ft yr⁻¹)², yielding an intracluster correlation coefficient of ρ = σ_b² / (σ_b² + σ_w²) = 0.317 ≈ 0.32.
- Including the standard unequal-cluster-size multiplier 1 + cv(m)² (Kish 1965) gives the forward design effect: DEFF = [1 + (m̄ − 1)ρ][1 + cv(m)²] = 13.83 × 2.10 = 29.1 ≈ 29.

Substituting the effective sample size n_eff = n / DEFF directly into Equation 1 of the main text—using the random-thinning calibration c = 1.27 and the archive SNR = 0.0602, with no free parameters—yields the forward predictions in Table S2.

| Retention f | Retained wells n | Effective size n_eff = n / 29.1 | Observed clustered P(flip) | Forward Kish prediction |
| :---: | :---: | :---: | :---: | :---: |
| 0.50 | 12,877 | 443 | 0.005 | 0.054 |
| 0.30 | 7,726 | 266 | 0.100 | 0.107 |
| 0.20 | 5,150 | 177 | 0.148 | 0.155 |
| 0.10 | 2,575 | 89 | 0.130 | 0.236 |
| 0.05 | 1,287 | 44 | 0.129 | 0.306 |
| 0.02 | 515 | 18 | 0.167 | 0.373 |
| 0.01 | 257 | 9 | 0.138 | 0.409 |

In the management-relevant range around the 10% risk threshold (f = 0.20–0.30), the zero-parameter Kish prediction matches the observed clustered flip rate within 0.7 percentage points (predicting 0.107 vs. 0.100 observed at f = 0.30, and 0.155 vs. 0.148 observed at f = 0.20). At high retention (f = 0.50), substituting n_eff = n / DEFF with fixed N omits the cluster-level finite-population correction 1 − f (applying the cluster sampling fraction f = 0.50 lowers the predicted flip rate from 0.054 to 0.012, close to the observed 0.005). At extreme thinning (f ≤ 0.10, where only 6 to 62 variable-sized blocks remain), the constant-DEFF approximation overpredicts flip risk because the few surviving large blocks dominate the composite hydrograph.

## S3. Companion model for loss of Mann–Kendall significance

Loss of statistical significance (two-sided Mann–Kendall p ≥ 0.05 on a thinned network whose full-network reference trend has p < 0.05) occurs at larger sample sizes than sign reversal. A natural question is whether significance loss follows a threshold-shifted companion form of Equation 1 in which the zero-slope boundary is replaced by the critical z-score required for p = 0.05, shifting SNR to SNR · (1 − z_crit / |z_ref|) with z_crit = 1.96.

Testing this companion expression across all 143 thinning trials in the 13 strata whose full-network reference trends are statistically significant (p < 0.05) shows that only 31% of observed significance-loss rates fall within their 95% binomial bands—even when an empirical scaling factor (w ≈ 0.75–0.79) or an explicit 1/√n interannual noise-growth term is fitted—compared with 56% binomial coverage for the sign-reversal law on the same trials. The reason is physical: the Mann–Kendall p-value of the composite hydrograph is governed by interannual water-level fluctuations (multi-year droughts, wet spells, and serial correlation) that are averaged down as 1/√n when more wells are combined in each year, rather than by spatial variance in long-term well slopes (S²). The finite-population sampling law is therefore restricted to trend direction (P(flip)), whereas significance loss is reported directly from the Monte Carlo simulations.

## S4. Deterministic versus random well-retention strategies

To evaluate whether simple consolidation rules outperform random retention when an agency must reduce a monitoring network, two deterministic ordering rules were compared against 50 random deletion sequences in a strong-signal aquifer (High Plains, N = 8,160, SNR = 0.111) and a heterogeneous aquifer with offsetting local trends (Central Valley, N = 183, SNR = 0.042):

1. **Longest-records-first retention:** wells are ranked by total number of distinct calendar years observed, keeping the longest-monitored wells at every retention fraction f.
2. **Spatial-spread-first retention:** wells are selected by farthest-point (minimax) distance sampling across the aquifer footprint, starting from the network centroid and iteratively adding the well farthest from the already-selected set.

## S5. Estimator decomposition, hold-out validation, and sensitivity checks

### S5.1 Estimator decomposition (c = 1.27 vs. simple mean of slopes)

To separate the finite-population variance approximation from the composite-hydrograph Theil–Sen estimator, the full pipeline was compared against a pure subset-mean benchmark in which each subsample's trend is the arithmetic mean of its n pre-computed well slopes (matching the exact Cochran variance S²(N − n)/(nN)). Across all strata, the subset-mean benchmark yields a median scale factor of c_bench = 1.05, confirming the theoretical c = 1 baseline on the empirical slope distributions. Comparing the subset-mean benchmark against the full composite-hydrograph Theil–Sen pipeline gives a median ratio c_bench / c_i = 0.796 (whose inverse is 1 / 0.796 = 1.26, matching the pooled fit c = 1.27). Spatial averaging across wells suppresses asynchronous interannual and gappy-record noise before the non-parametric temporal Theil–Sen median slope is extracted, giving the composite estimator higher statistical efficiency (c > 1) than an arithmetic average of noisy single-well slopes. In the High Plains, where only the f = 0.01 trial (n = 81) has a non-zero flip rate (0.015), the single-point stratum fit is c_i = 1.28, which predicts n_min = 80 versus the simulated 81 and the pooled prediction of 82. Conversely, replacing the deseasonalized Theil–Sen procedure with ordinary least-squares regression on raw annual means shifts c to 0.17–0.19 because unscreened outliers inflate S, underscoring that design-chart inputs (β_net and S) must be computed with the deseasonalized Theil–Sen pipeline of Section 2.2.

### S5.2 Out-of-sample validation and residual diagnostics

To test out-of-sample generalizability, the eligible wells in each stratum were split into odd- and even-indexed halves, calibrating c on one half and evaluating binomial coverage on the other; across both directions and 20 random 11-for-10 stratum splits, c ranges from 1.21 to 1.33 and mean out-of-sample binomial coverage is 0.49, virtually identical to the in-sample all-trial coverage of 0.50. Across the 20 strata with estimable c_i, correlating c_i against well-slope skewness, kurtosis, balanced-panel composition, or spatial autocorrelation (Moran's I, which is statistically significant in 18 of 21 strata with a median of I = 0.33) yields Spearman p > 0.18 in all cases. Replacing the normal tail Φ(−x) with a Student-t₅ link shifts overall binomial coverage by less than 2 percentage points (49.1% vs. 50.5%).

### S5.3 Sensitivity to spatial aggregation, panel balance, and time-window start year

Across the 22 reference series (21 strata plus CONUS), alternative constructions of the 1975–2025 reference trend preserve the reference sign in 15 to 20 series:

- **Spatial-median composite hydrograph:** preserves the reference sign in 20 of 22 series (exceptions: Floridan and Rio Grande).
- **Inverse-density weighting (0.5° grid cells):** preserves the reference sign in 19 of 22 series (exceptions: Floridan, Central Valley, and Columbia Plateau).
- **Nearly balanced panel (restricting to the 4,295 wells observed in ≥ 45 of the 51 years):** preserves the reference sign in 15 of 20 testable strata.
- **Shifting the common analysis window start year to 1970–2025 or 1980–2025:** preserves the reference sign in 17 of 21 strata (the four strata that change sign—Basin and Range, Floridan, Surficial, and Cambrian–Ordovician—all have weak reference trends with well-bootstrap 95% confidence intervals that cross zero in Table 1).

All analysis scripts and derived datasets supporting Sections S1–S5 are archived in the `scripts/` and `data/` directories of the public repository (https://github.com/yuelangmanle/gwat_wellnet).
