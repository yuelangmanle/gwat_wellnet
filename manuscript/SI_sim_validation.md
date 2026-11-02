# Supporting Information for "Monitoring Network Density and Groundwater-Level Trend Reliability: A Closed-Form Law"

## S1. Synthetic-population testbed for the finite-population sampling law

To test whether the sampling law in Equation 1 of the main text holds under known population parameters—and to separate distributional skewness and interannual water-level noise from the composite-hydrograph Theil–Sen estimator—synthetic monitoring networks were generated and thinned under controlled conditions. Each synthetic network consists of N = 400 wells observed annually over T = 51 years (1975–2025). Individual well trends b_i were drawn with population standard deviation S = 1 and population mean β_true = SNR · S = 0.15 (SNR = 0.15) across three standardized slope distributions:

1. **Normal (Gaussian):** b_i ~ Normal(β_true, S²).
2. **Skewed (lognormal):** standardized lognormal slopes with skewness ≈ 2, scaled to mean β_true and variance S².
3. **Offsetting mixture:** a bimodal mixture of 60% depleting wells centered at +1.0 S and 40% recovering wells centered at −0.9 S (within-mode standard deviation 0.365 S, total variance S²), representing an aquifer with both depleting and recovering sub-basins.

Each distribution was evaluated in two observation modes across six retention fractions (f ∈ {0.05, 0.08, 0.12, 0.18, 0.28, 0.45}; 500 Monte Carlo replicates per cell) using the composite-hydrograph Theil–Sen pipeline of the main text:

- **Noiseless mode:** each well follows y_{i,t} = b_i(t − t̄) without interannual noise, isolating pure spatial sampling of the N well slopes.
- **Annual-noise mode:** each well is assigned independent Gaussian annual anomalies ε_{i,t} chosen so the per-well slope estimation standard error equals 0.3 S, comparable to empirical residual variability.

For each cell, the apparent scale factor c_app was fitted by binomial maximum likelihood to P(flip) = Φ[−c_app · SNR · √(nN/(N − n))] across informative retention levels (Table S1).

**Table S1.** Apparent scale factor c_app fitted across synthetic slope distributions and observation modes (N = 400, T = 51 yr, SNR = 0.15).

| Slope distribution | Observation mode | Informative trials | Fitted c_app |
| :--- | :--- | :---: | :---: |
| Normal (Gaussian) | Noiseless | 6 | 1.01 |
| Skewed (lognormal, skew ≈ 2) | Noiseless | 5 | 1.04 |
| Offsetting mixture (60% / 40%) | Noiseless | 6 | 1.03 |
| Normal (Gaussian) | Annual noise (slope SE = 0.3 S) | 6 | 0.93 |
| Skewed (lognormal, skew ≈ 2) | Annual noise (slope SE = 0.3 S) | 6 | 0.95 |
| Offsetting mixture (60% / 40%) | Annual noise (slope SE = 0.3 S) | 6 | 0.91 |

Two conclusions follow from Table S1:

1. **In the noiseless regime, c_app = 1.01–1.04 across all three slope distributions.** The finite-population factor √(nN/(N − n)) and the normal tail Φ(−x) accurately describe sign-reversal probability even when individual well slopes are strongly skewed or bimodal.
2. **Adding independent interannual water-level noise slightly reduces c_app to 0.91–0.95**, because residual temporal noise adds variance σ_e² / (n ∑(t − t̄)²) to the spatial slope variance S²(N − n)/(nN). Consequently, the empirical value c = 1.27 found for the U.S. Geological Survey (USGS) archive is not an artifact of slope skewness or random annual noise; as shown in Section S5, it arises from the statistical efficiency of applying the non-parametric Theil–Sen estimator to the spatially averaged anomaly hydrograph across gappy, asynchronous well records.

## S2. Zero-parameter forward check of spatially clustered thinning

In the national (CONUS) clustered-thinning experiment, wells were grouped into 0.5° × 0.5° latitude–longitude blocks (retaining blocks with at least 10 wells; B = 621 blocks containing N = 25,754 wells), and whole blocks were removed at random. To verify that the elevated sign-reversal rate under clustered removal follows directly from survey-sampling theory without refitting any parameter, the design effect (DEFF; Kish 1965) was computed from the archive's block-level slope statistics:

- Mean block size: m̄ = 41.5 wells per block; block-size coefficient of variation: cv(m) = 1.04.
- One-way random-effects analysis of variance on individual well slopes across blocks gives between-block variance σ_b² = 0.00674 (m yr⁻¹)² and within-block variance σ_w² = 0.01452 (m yr⁻¹)², yielding an intracluster correlation coefficient of ρ = σ_b² / (σ_b² + σ_w²) = 0.317 ≈ 0.32.
- Including the standard unequal-cluster-size multiplier 1 + cv(m)² (Kish 1965) gives the forward design effect: DEFF = [1 + (m̄ − 1)ρ][1 + cv(m)²] = 13.96 × 2.08 = 29.06 ≈ 29.

Substituting the effective sample size n_eff = n / DEFF directly into Equation 1 of the main text—using the random-thinning calibration c = 1.27 and the archive SNR = 0.060 (0.05955), with no free parameters—yields the forward predictions in Table S2.

**Table S2.** Observed sign-reversal probabilities under 0.5° spatially clustered removal in CONUS versus zero-parameter forward predictions using the Kish (1965) design effect (DEFF = 29.06).

| Retention f | Retained wells n | Effective size n_eff = n / 29.06 | Observed clustered P(flip) | Forward Kish prediction |
| :---: | :---: | :---: | :---: | :---: |
| 0.50 | 12,877 | 443 | 0.005 | 0.054 |
| 0.30 | 7,726 | 266 | 0.101 | 0.108 |
| 0.20 | 5,150 | 177 | 0.151 | 0.156 |
| 0.10 | 2,575 | 89 | 0.135 | 0.238 |
| 0.05 | 1,287 | 44 | 0.148 | 0.307 |
| 0.03 | 772 | 27 | 0.166 | 0.348 |
| 0.02 | 515 | 18 | 0.180 | 0.375 |
| 0.01 | 257 | 9 | 0.224 | 0.411 |

In the management-relevant range around the 10% risk threshold (f = 0.20–0.30), the zero-parameter Kish prediction matches the observed clustered flip rate within 0.7 percentage points (predicting 0.108 vs. 0.101 observed at f = 0.30, and 0.156 vs. 0.151 observed at f = 0.20). At high retention (f = 0.50), substituting n_eff = n / DEFF with fixed N omits the cluster-level finite-population correction 1 − f (applying the cluster sampling fraction f = 0.50 lowers the predicted flip rate from 0.054 to 0.012, close to the observed 0.005). At extreme thinning (f ≤ 0.10, where only 6 to 62 variable-sized blocks remain), the constant-DEFF approximation overpredicts flip risk because the few surviving large blocks dominate the composite hydrograph.

## S3. Companion model for loss of Mann–Kendall significance

Loss of statistical significance (two-sided Mann–Kendall p ≥ 0.05 on a thinned network whose full-network reference trend has p < 0.05) occurs at larger sample sizes than sign reversal. A natural question is whether significance loss follows a threshold-shifted companion form of Equation 1 in which the zero-slope boundary is replaced by the critical z-score required for p = 0.05, shifting SNR to SNR · (1 − z_crit / |z_ref|) with z_crit = 1.96.

Testing this companion expression across all 143 thinning trials in the 13 strata whose full-network reference trends are statistically significant (p < 0.05) shows that only 31% of observed significance-loss rates fall within their 95% binomial bands—even when an empirical scaling factor (w ≈ 0.75–0.79) or an explicit 1/√n interannual noise-growth term is fitted—compared with 56% binomial coverage for the sign-reversal law on the same trials. The reason is physical: the Mann–Kendall p-value of the composite hydrograph is governed by interannual water-level fluctuations (multi-year droughts, wet spells, and serial correlation) that are averaged down as 1/√n when more wells are combined in each year, rather than by spatial variance in long-term well slopes (S²). The finite-population sampling law is therefore restricted to trend direction (P(flip)), whereas significance loss is reported directly from the Monte Carlo simulations.

## S4. Deterministic versus random well-retention strategies

To evaluate whether simple consolidation rules outperform random retention when an agency must reduce a monitoring network, two deterministic ordering rules were compared against 50 random deletion sequences in a strong-signal aquifer (High Plains, N = 8,160, SNR = 0.111) and a heterogeneous aquifer with offsetting local trends (Central Valley, N = 183, SNR = 0.042):

1. **Longest-records-first retention:** wells are ranked by total number of distinct calendar years observed, keeping the longest-monitored wells at every retention fraction f.
2. **Spatial-spread-first retention:** wells are selected by farthest-point (minimax) distance sampling across the aquifer footprint, starting from the network centroid and iteratively adding the well farthest from the already-selected set.

At each 2% retention step, the composite anomaly hydrograph over 1975–2025 and its Theil–Sen slope β̂ were recomputed (Figure S1). In the High Plains aquifer (Figure S1a), none of the three ordering rules reverses the reference depletion sign at any tested sample size; however, at 10% retention (816 wells), keeping the longest records first produces 5.4 times the absolute slope bias |β̂ − β_net| of random retention (0.017 vs. 0.003 m yr⁻¹) because the oldest monitoring wells are concentrated in early, heavily pumped irrigation districts. In the Central Valley (Figure S1b), random and longest-records-first deletions reverse the reference trend sign below ~150 wells, whereas spatial-spread-first retention preserves the reference trend sign down to 14 wells (8% retention) at the cost of inflating the composite slope to 3.1 times the full-network reference value.

## S5. Estimator decomposition, hold-out validation, and sensitivity checks

### S5.1 Estimator decomposition (c = 1.27 vs. simple mean of slopes)

To separate the finite-population variance approximation from the composite-hydrograph Theil–Sen estimator, the full pipeline was compared against a pure subset-mean benchmark in which each subsample's trend is the arithmetic mean of its n pre-computed well slopes (matching the exact Cochran variance S²(N − n)/(nN)). Across all strata, the subset-mean benchmark yields a median scale factor of c_bench = 1.05, confirming the theoretical c = 1 baseline on the empirical slope distributions. Comparing the subset-mean benchmark against the full composite-hydrograph Theil–Sen pipeline gives a median ratio c_bench / c_i = 0.796 (whose inverse is 1 / 0.796 = 1.26, matching the pooled fit c = 1.27). Spatial averaging across wells suppresses asynchronous interannual and gappy-record noise before the non-parametric temporal Theil–Sen median slope is extracted, giving the composite estimator higher statistical efficiency (c > 1) than an arithmetic average of noisy single-well slopes. In the High Plains, where only the f = 0.01 trial (n = 81) has a non-zero flip rate (0.015), the single-point stratum fit is c_i = 1.28, which predicts n_min = 80 versus the simulated 81 and the pooled prediction of 82. Conversely, replacing the deseasonalized Theil–Sen procedure with ordinary least-squares regression on raw annual means shifts c to 0.17–0.19 because unscreened outliers inflate S, underscoring that design-chart inputs (β_net and S) must be computed with the deseasonalized Theil–Sen pipeline of Section 2.2.

### S5.2 Out-of-sample validation and residual diagnostics

Across 200 random half-splits of the 22 networks (calibrating c on 11 networks and evaluating binomial coverage on the held-out 11 networks), the median split-sample fit is c = 1.24 (range: 1.11–1.53) and mean out-of-sample binomial coverage is 0.49, virtually identical to the in-sample all-trial coverage of 0.50. Across the 20 strata with estimable c_i, correlating c_i against well-slope skewness, kurtosis, balanced-panel composition, or spatial autocorrelation (Moran's I, which is statistically significant in 18 of 21 strata with a median of I = 0.33) yields Spearman p > 0.18 in all cases. Replacing the normal tail Φ(−x) with a heavy-tailed Student-t₃ link (ν = 3) reduces transition-zone binomial coverage from 23% to 16% and overall coverage from 50.5% to 36.6%.

### S5.3 Sensitivity to spatial aggregation, panel balance, and time-window start year

Across the 22 reference series (21 strata plus CONUS), alternative constructions of the 1975–2025 reference trend preserve the reference sign in 15 to 20 series:

- **Spatial-median composite hydrograph:** preserves the reference sign in 20 of 22 series (exceptions: Floridan and Rio Grande).
- **Inverse-density weighting (0.5° grid cells):** preserves the reference sign in 19 of 22 series (exceptions: Floridan, Central Valley, and Columbia Plateau).
- **Nearly balanced panel (restricting to the 5,492 wells observed in ≥ 45 of the 51 years):** preserves the reference sign in 15 of 21 aquifer strata.
- **Shifting the common analysis window start year to 1970–2025 or 1980–2025:** preserves the reference sign in 17 of 21 strata (the four strata that change sign—Basin and Range, Floridan, Surficial, and Mississippi embayment—all have weak reference trends with well-bootstrap 95% confidence intervals that include zero in Table 1).

### S5.4 USGS national aquifer codes and official stratum names

Table S3 lists the U.S. Geological Survey National Water Information System (NWIS) national aquifer codes (nat_aqfr_cd) and official USGS principal aquifer names corresponding to the 21 monitoring-well strata reported in Table 1 of the main text.

Table S3. Mapping between Table 1 stratum labels, USGS NWIS national aquifer codes (nat_aqfr_cd), official USGS principal aquifer names, and eligible long-term monitoring well counts (N).

| Table 1 stratum label | USGS code (nat_aqfr_cd) | Official USGS national aquifer name | Eligible wells (N) |
| :--- | :--- | :--- | ---: |
| High Plains | N100HGHPLN | High Plains aquifer | 8,160 |
| Unassigned wells | (unassigned) | Wells without an assigned USGS national aquifer code | 6,932 |
| Northern Atlantic Coastal Plain | S100NATLCP | Northern Atlantic Coastal Plain aquifer system | 1,959 |
| Basin and Range basin-fill | N100BSNRGB | Basin and Range basin-fill aquifers | 1,557 |
| Sand and gravel (glaciated) | N100GLCIAL | Sand and gravel aquifers (glaciated regions) | 1,422 |
| Mississippi River Valley alluvial | N100MSRVVL | Mississippi River Valley alluvial aquifer | 1,237 |
| Coastal lowlands | S100CSLLWD | Coastal lowlands aquifer system | 1,145 |
| Alluvial aquifers | N100ALLUVL | Alluvial aquifers | 1,101 |
| California Coastal Basins | N100CACSTL | California Coastal Basin aquifers | 863 |
| Floridan | S400FLORDN | Floridan aquifer system | 819 |
| Other aquifers | N9999OTHER | Other aquifers (designated minor or local aquifers) | 765 |
| Rio Grande | S100RIOGRD | Rio Grande aquifer system | 742 |
| Mississippi embayment | S100MSEMBM | Mississippi embayment aquifer system | 477 |
| Biscayne | N400BISCYN | Biscayne aquifer | 232 |
| Surficial aquifer system | S100SURFCL | Surficial aquifer system | 195 |
| Intermediate aquifer system | S500INTRMD | Intermediate aquifer system | 186 |
| Central Valley | S100CNRLVL | Central Valley aquifer system | 183 |
| Roswell Basin | S400RSWLBS | Roswell Basin aquifer system | 155 |
| Colorado Plateaus | N300COPLTS | Colorado Plateaus aquifers | 111 |
| Cambrian-Ordovician | S300CAMORD | Cambrian-Ordovician aquifer system | 102 |
| Columbia Plateau basaltic-rock | N600CMBPLV | Columbia Plateau basaltic-rock aquifers | 100 |

All analysis scripts and derived datasets supporting Sections S1–S5 are archived in the public repository (https://github.com/yuelangmanle/gwat_wellnet).
