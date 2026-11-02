# Monitoring Network Density and Groundwater-Level Trend Reliability: A Closed-Form Law

## Authors

Zijie Deng

Jiangxi Provincial Key Laboratory of Water Resources Allocation and Efficient Utilization, Jiangxi University of Water Resources and Electric Power, No. 289 Tianxiang Avenue, Nanchang 330099, China

Corresponding author: Zijie Deng, 3335196397vip@gmail.com (ORCID 0009-0002-4824-9764).

**Running head:** Network density and trend reliability

---

## Abstract

Groundwater-level records provide primary evidence for assessing aquifer depletion, yet monitoring networks face persistent fiscal contraction. This study evaluates how reliably downsized networks preserve the regional trend direction of full networks. From 8.8 million quality-controlled U.S. Geological Survey depth-to-water measurements (1950–2026), 29,717 long-term monitoring wells were identified across the conterminous United States and grouped into 21 strata. Regional trends and network thinning were evaluated over a common 51-year window (1975–2025). Virtual thinning across 1,000 Monte Carlo replicates per scenario (216 distinct scenarios across strata from 90% down to 1% retention) demonstrates that trend sign-reversal probability follows a finite-population sampling law. This relationship is governed by network retention and the signal-to-noise ratio (SNR, reference trend magnitude divided by individual well slope dispersion). Inverting this relationship at a 10% reversal risk yields a closed-form minimum network size, n_min = Nz²/(z²+(c·SNR)²N) with c = 1.27 (95% CI: 1.04–1.55), reproducing simulated minima within a factor of two for 17 of 19 resolvable networks. Strong-signal aquifers preserve regional trend sign with as few as 11 wells (Coastal Lowlands) or 81 of 8,160 wells in the High Plains (representing regional sign thresholds rather than local management targets), whereas systems with offsetting trends (Basin and Range, surficial, Floridan) require 88–94% of their wells. Spatially clustered cuts are substantially more damaging: the national trend tolerates random thinning down to 297 wells (1%), whereas retaining 7,726 wells (30%) in spatial clusters incurs an observed 10.1% reversal rate, equivalent to only ~180–280 randomly distributed wells.

## Impact Statement

A closed-form sampling law predicts the minimum network size needed to preserve regional groundwater trend direction under well loss.

## Keywords

groundwater monitoring; water-level trends; network design; trend detection; depletion; Theil–Sen

---

## 1. Introduction

Observation-well records provide primary evidence for groundwater depletion (Konikow 2015; Jasechko et al. 2024), satellite gravimetry validation (Rodell et al. 2018), and national assessments (Stets et al. 2025). Yet monitoring networks face persistent funding cuts and attrition. Network designs optimize spatial interpolation (Ohmer et al. 2019) or seasonal dynamics (Stahl et al. 2024), but rarely examine how contraction alters the *trend conclusion* itself: if downsized, how often does the retained subset reverse regional trend direction or lose significance? Addressing this is urgent because sparse networks detect trends in only a minority of areas (Nolte et al. 2026), and mapping can outrun underlying field observations (Bradbury 2025).

This study investigates how monitoring density controls regional trend reliability across CONUS using 8.8 million measurements from 29,717 long-term wells (1950–2026), evaluating composite trends and thinning over 1975–2025 across three questions: (1) How often does contraction reverse regional trend sign or lose significance under random versus clustered removal? (2) Can reversal risk be predicted by a closed-form sampling law from archive statistics? (3) Does retaining longest records or spatially dispersed wells outperform random thinning?

## 2. Data and Methods

### 2.1 Groundwater-level data and analysis time windows

Discrete depth-to-water measurements (USGS parameter code 72019, 1950–2026) were downloaded from the U.S. Geological Survey Water Data API (U.S. Geological Survey 2016), yielding 9,456,561 raw measurements from 843,302 sites. Quality assurance retained approved records within 0–1,524 m (0–5,000 ft), removed pumping and dry-well flags (234,752 excluded), aggregated same-day measurements to daily medians (8,814,831 well-days), and excluded sites outside CONUS (103,049 excluded). Requiring at least 30 measurements across at least 20 calendar years yielded 29,717 eligible wells (3.5% of all sites) across 48 conterminous states and Washington, D.C.; 3,064 wells (10.3%) share coordinates with another well (nested piezometers or replacements) and were retained as separate series. Grouping wells by USGS principal aquifer code identifies 21 strata with at least 100 eligible wells (Table 1): 20 principal aquifers plus 6,932 unassigned wells.

Because record spans vary, two complementary windows are analyzed: (1) full archive (1950–2026, ≥ 20 annual values) screens eligible wells and calculates single-well slope dispersion S; (2) a common 51-year window (1975–2025) constructs regional anomaly hydrographs and executes thinning experiments under dense nationwide monitoring (29,332 wells, 98.7%, have active observations in 1975–2025; median coverage is 29 of 51 years, and 3.5% cover all 51 years).

### 2.2 Well-scale trend estimation

To reduce seasonal sampling bias (Hirsch et al. 1982), each well's record was deseasonalized by subtracting month-of-year median depth to water for months with at least two observations (sparser months were centered on the overall median), aggregating residuals to annual medians. Trends were estimated with the non-parametric Theil–Sen slope (Theil 1950; Sen 1968). Significance (two-sided p < 0.05) was tested via the Mann–Kendall test (Mann 1945; Kendall 1975) using Yue–Pilon trend-free pre-whitening (Yue et al. 2002) where lag-1 serial correlation was significant at the 95% level (60.2% of wells). Positive depth-to-water slopes indicate declining water levels; 62.3% of the 29,717 wells have positive slopes and 70.1% are statistically significant (p < 0.05).

### 2.3 Network-thinning experiments

For each of the 21 strata and CONUS (22 series; CONUS serves as a national statistical benchmark), each well's annual series was centered on its 1975–2025 median and averaged across wells in each calendar year to form the regional anomaly hydrograph; the reference trend (β_net) is the 1975–2025 Theil–Sen slope of this composite series. Alternative aggregation schemes—spatial medians, inverse-density weighting, balanced panels, or shifting window starts between 1970 and 1980—preserve the reference sign in 15 to 20 of 22 series, with changes limited to weak-signal strata whose bootstrap intervals include zero (Table 1; Supporting Information Section S5).

Network contraction was simulated across 11 retention fractions (f = 0.90 to 0.01). At each fraction, n = max(10, ⌊fN⌋) wells were drawn without replacement 1,000 times to recompute composite slopes and Mann–Kendall p-values over 1975–2025. A replicate was defined as a *sign flip* if its slope had the opposite sign to β_net, and as a *loss of significance* if the reference trend had p < 0.05 but the subsample had p ≥ 0.05 or flipped sign. In *uniform random thinning*, wells were selected with equal probability. In *spatially clustered thinning*, wells were grouped into 0.5° × 0.5° blocks (excluding blocks with < 10 wells, leaving N = 25,754 in CONUS) and whole blocks were deleted at random. The minimum reliable network size (n_min) is the smallest n with sign-reversal probability ≤ τ (τ = 0.10). Because subsampling ceases at 10 wells, floor values represent upper bounds; pooling duplicate floor sizes in small strata leaves 216 distinct random-thinning trials. In the High Plains and Central Valley, 50 random deletion orders were compared against *longest-records* and *spatial-spread* retention.

### 2.4 Finite-population sampling law for sign-reversal risk

Under uniform random selection of n from N wells without replacement, classical finite-population sampling theory applies (Cochran 1977; Supporting Information Section S1). Subset sampling variance depends on the finite-population variance S² of the N individual well slopes. To first order, the composite anomaly slope behaves like the sample mean of the n retained slopes, with sampling variance S²(N − n)/(nN). Approximating the subsample slope distribution around β_net as normal yields a one-parameter law for sign-reversal probability:

  P(flip) = Φ[−c · SNR · √(nN/(N − n))],   SNR = |β_net| / S,   (1)

where Φ is the standard normal cumulative distribution function and c is an empirical efficiency factor (c = 1.27, 95% cluster-bootstrap confidence interval: 1.04–1.55; leave-one-out range: 1.24–1.42; stratum fits c_i span 0.55 to 4.05). The signal-to-noise ratio SNR = |β_net| / S divides 1975–2025 reference trend magnitude by the standard deviation S of 1950–2026 individual well slopes (≥ 20 annual values). Estimator comparisons indicate that c > 1 reflects the higher precision of composite-hydrograph Theil–Sen slopes relative to individual well slope dispersion (Supporting Information Section S5), where spatial averaging attenuates asynchronous interannual noise before trend fitting. Design-chart inputs must be evaluated using the deseasonalized Theil–Sen pipeline of Sections 2.2–2.3.

Solving for n at acceptable reversal risk τ (z = Φ⁻¹(1 − τ); τ = 0.10, z = 1.282) provides the closed-form minimum network size:

  n_min = N z² / (z² + (c·SNR)² N).   (2)

Where regional signals are strong ((c·SNR)²N ≫ z²), N cancels and n_min ≈ (z/(c·SNR))² depends solely on SNR. Where signals are weak ((c·SNR)²N ≲ z²), n_min approaches N because localized opposing trends cancel. For clustered removal, inverting Equation 1 yields effective random sample size n_eff.

### 2.5 Software

Analyses were performed in Python 3.9 (pandas, NumPy, SciPy, Matplotlib).

## 3. Results

### 3.1 National monitoring network and reference trends

Monitoring wells concentrate in the High Plains, California, Florida, the Mississippi embayment, and the glaciated north (Figure 1a). Water levels deepen at 18,513 wells (62.3%), rise at 11,134 (37.5%), and remain stable at 70 (0.2%); median slope is +0.010 m yr⁻¹ (interquartile range: −0.010 to +0.083; extremes: −9.8 to +5.4 m yr⁻¹). Across 21 strata, 1975–2025 reference slopes range from −0.170 to +0.074 m yr⁻¹ (Table 1); in Coastal Lowlands, regional composite water levels indicate recovery (−0.170 m yr⁻¹) following post-1970s reductions in coastal pumping, contrasting with historical declines across 53.1% of wells over 1950–2026.

### 3.2 Sign-reversal probability under random and clustered thinning

Under random thinning, sign-reversal curves separate distinctly by aquifer SNR (Figure 2a). High-SNR systems maintain reference trend signs down to small fractions of wells (Coastal Lowlands, n_min = 11; Mississippi River Valley alluvium, 12), whereas the Basin and Range aquifer (SNR = 0.006, reference slope +0.002 m yr⁻¹) flips sign in up to one-third of replicates across all tested sizes.

On CONUS, random and clustered thinning diverge sharply (Figure 2b). Random removal produces zero sign flips down to 5% retention (1,485 wells) and a 0.075 flip rate at the 1% floor (297 of 29,717 wells). Under 0.5° clustered deletion, the reversal rate reaches an observed 0.101 (10.1%) at 30% retention of block-eligible wells (7,726 of 25,754 wells; 95% binomial CI: 8.3–12.1%), approximately the 10% risk threshold, equivalent to only n_eff ≈ 280 random wells (95% CI: 236–331); across 10–30% clustered retention, effective sample size remains constrained at ~180–280 random wells. A cluster sampling model (Kish 1965) based on block statistics (ρ = 0.32, mean block size 41.5, cv = 1.04) yields design effect DEFF ≈ 29; setting n_eff = n / 29 predicts observed flip rates within ~1 percentage point across 20–30% retention (Supporting Information Section S2). Using 0.25° or 1.0° blocks yields CONUS flip rates of 0.11–0.15 at 30% retention. At the aquifer scale, 30% clustered retention preserves stability in the High Plains (0.002 flip rate) but elevates flip rates to 0.31 in Coastal Lowlands and 0.17 in Central Valley.

Statistical significance (p < 0.05) erodes at larger sample sizes than trend direction (at 1% in CONUS, 18.5% of replicates lose significance while 7.5% flip sign; at 30% in glaciated aquifers, 36% lose significance and none flips). Significance loss is governed by interannual noise growth in the composite hydrograph rather than spatial slope variance; a threshold-shifted companion approximation of Equation 1 captures only 31% of trials within 95% binomial bands (Supporting Information Section S3).

### 3.3 Minimum reliable network sizes across aquifers

Simulated minimum reliable sizes (n_min for ≤ 10% flip risk; Fig. 1b and Table 1) categorize networks into three regimes. In the strong-signal regime, six strata and CONUS remain reliable down to the grid floor of 1% of N or 10 wells (Coastal Lowlands at 11 wells, Mississippi alluvium at 12, High Plains at 81, and CONUS at 297, all upper bounds; Table 1), while three coastal and glaciated strata stabilize between 2% and 5% retention (25–71 wells). In the intermediate regime, nine heterogeneous aquifers require 20–80% retention: Biscayne requires ~20%, whereas Rio Grande (519 of 742 wells), Mississippi embayment (333 of 477), Cambrian–Ordovician (71 of 102), and Central Valley (146 of 183) require 70–80%, where local declines up to 1.4 m yr⁻¹ are offset by recovery areas. Rio Grande and Roswell Basin composite series show significant temporal trends (p < 0.05) despite bootstrap intervals spanning zero, as localized declines and recoveries offset across wells (Table 1). In the weak-signal regime, three aquifers with offsetting local trends (Basin and Range, surficial, and Floridan) never reach ≤ 10% reversal risk on the standard grid; Equation 2 predicts required retentions of 94%, 91%, and 88%; dense simulations at 92–98% retention confirm the first two within ~2 percentage points while placing the Floridan threshold near 93% (versus 88% predicted).

Across the 19 resolvable networks, Equation 2 reproduces simulated n_min within a factor of two for 17 networks (89%; median log₁₀ error 0.15 dex, comparable to the 1.5× median grid step); glaciated aquifers (213 predicted vs. 71 simulated) and Roswell Basin (17 vs. 77) exhibit the two larger deviations, reflecting localized efficiency factors (c_i = 2.75 and 0.55) departing from pooled c = 1.27.

### 3.4 Collapse across networks and the design chart

Plotting all 216 random-thinning trials against normalized sample size x = c · SNR · √(nN/(N − n)) collapses observations toward the standard normal curve Φ(−x) with no fitted shape parameter (Figure 3a). In the low-risk regime (< 5% flip rate), 74% of trials fall within 95% binomial Monte Carlo bands. In the transition zone (23% coverage), slope skewness introduces departures between composite Theil–Sen slopes and simple means (Supporting Information Section S5). In strongly skewed Rio Grande and Floridan systems, seven deep-thinning trials reach flip rates of 0.51–0.59 outside theoretical bands, all falling within the low-SNR zone identified as unreliable. Contours of required retention n_min / N across (SNR, N) provide a direct design chart (Figure 3b).

### 3.5 Deterministic retention strategies: slope bias versus sign risk

In the High Plains (N = 8,160) and Central Valley (N = 183), deterministic retention alters slope bias rather than sign-reversal risk (Figure S1; Supporting Information Section S4). In the High Plains, no strategy flips the depletion sign, but at 10% retention (816 wells) retaining longest records produces 5.4× the slope bias of random retention (0.017 vs. 0.003 m yr⁻¹) due to clustering in early pumping districts. In the Central Valley, spatial-spread retention preserves the reference sign down to 8% retention (14 wells) while inflating slope magnitude to 3.1× the reference value.

## 4. Discussion

**Designing and consolidating monitoring networks.** For *network size*, agencies can estimate β_net and S via the deseasonalized Theil–Sen pipeline, evaluate SNR = |β_net| / S, and determine n_min from Figure 3b. Unlike geostatistical designs minimizing kriging variance (Ohmer et al. 2019), this formulation identifies the threshold required to preserve composite trend direction; localized tasks such as cone-of-depression tracking require denser monitoring. In vast systems such as the High Plains (~450,000 km²), n_min = 81 represents an aquifer-scale sign threshold rather than an operational target for Groundwater Management Districts managing localized sub-basin depletions (e.g., severe depletion in western Kansas and Texas versus relative stability in parts of Nebraska). Where (c·SNR)²N ≫ z², modest well counts preserve regional trend sign; where (c·SNR)²N ≲ z², localized trends cancel and sub-basins should be analyzed separately. In direct application without local simulation, a 1.5–2× safety factor accommodates regional noise synchrony and slope skewness (covering stratum efficiencies down to c_i ≈ 0.7). For *spatial layout*, well-slope autocorrelation within 0.5° blocks (ρ = 0.32) inflates sampling variance by Kish design effect DEFF ≈ 29 under district-scale cuts, requiring reductions to be distributed across the domain. For *measurement frequency*, annual subsampling maintains median annual-to-full-archive SNR ratio near unity (1.06, interquartile range 0.82–1.14), but degrades SNR in seven strata into high-retention regimes ((c·SNR)²N ≲ z²), indicating that stratum-specific SNR should be checked before reducing site visit frequency.

**Implications for regional syntheses.** Regional compilations (Konikow 2015; Rodell et al. 2018; Jasechko et al. 2024; Stets et al. 2025) inherit underlying network vulnerabilities: in low-SNR aquifers, retaining half the domain carries 10–35% sign-reversal risk. Analogous to time-window sensitivity (Toran et al. 2026), regional trends depend intrinsically on spatial density, for which archive SNR provides a transparent metric.

**Scope and limitations.** Depth-to-water records (parameter 72019) maximize national coverage; co-located nested piezometers (10.3% of wells) were retained as separate series. Stratifying unstable systems by USGS aquifer type confirms that low SNR reflects opposing geographical trends rather than confined-unconfined mixing, though local investigations should incorporate head. Depth to water reflects levels relative to land surface; in compacting unconsolidated basins (e.g., Central Valley and Coastal Lowlands), subsidence introduces datum drift unless tied to benchmarks. In the Central Valley (N = 183), this count reflects federal NWIS records meeting our filtering criteria, whereas state portals (e.g., California SGMA) compile thousands of wells for compliance. Where regional trends are near zero, reference signs can vary with spatial aggregation rules, flagged whenever bootstrap intervals include zero (Table 1). Preferential well attrition (e.g., shallow wells going dry) introduces bias analogous to longest-record retention. Within its calibrated range (N ≥ 100, n ≥ 10–15, τ = 0.10, Section 2.2 pipeline), Equation 2 serves as a reliable first-order benchmark; Monte Carlo simulations remain advisable in transition aquifers.

## 5. Conclusions

1. Across 29,717 long-term U.S. monitoring wells (1950–2026) evaluated over 1975–2025, regional trend sign-reversal probability follows a finite-population sampling law, P(flip) = Φ[−c · SNR · √(nN/(N − n))] with c = 1.27 (95% CI: 1.04–1.55). Its closed-form inverse, n_min = Nz²/(z²+(c·SNR)²N), predicts minimum reliable network sizes within a factor of two for 17 of 19 resolvable networks. Significance erodes at larger sizes than trend direction due to interannual noise.

2. Minimum reliable network sizes for ≤ 10% reversal risk exhibit three regimes. Six strong-signal strata and CONUS remain reliable down to the 1% or 10-well floor (including 11 wells in Coastal Lowlands and 81 of 8,160 in High Plains, both upper bounds). Nine heterogeneous strata require 20–80% retention (including 70–80% in Rio Grande and Central Valley). Three weak-signal systems (Basin and Range, surficial, Floridan) require 88–94% retention.

3. Clustered well losses are far more damaging than random thinning. The national trend tolerates uniform random deletion down to 297 of 29,717 wells (1%), whereas 0.5° clustered removal reaches an observed 10.1% reversal risk (approximately the 10% threshold) at 7,726 wells (30%), equivalent to only ~180–280 random wells (Kish design effect ≈ 29). Network consolidation must be distributed geographically rather than concentrated by district; agencies should avoid retaining longest-running records alone, which concentrates monitoring in early pumping districts and inflates slope bias.

## Acknowledgments

The author thanks the U.S. Geological Survey for maintaining the National Water Information System (NWIS).

## Author's Note

The author has no conflicts of interest or financial disclosures to report.

## Data Availability Statement

Groundwater-level measurements were retrieved from the U.S. Geological Survey Water Data API (parameter code 72019), accessed October 2026. Analysis scripts and derived data are available at https://github.com/yuelangmanle/gwat_wellnet.

## References

Bradbury, K.R. 2025. Groundwater potential mapping: A misused and dubious concept. *Groundwater* 63, no. 5: 664–665. https://doi.org/10.1111/gwat.70013

Cochran, W.G. 1977. *Sampling Techniques*. 3rd ed. New York: John Wiley & Sons.

Hirsch, R.M., J.R. Slack, and R.A. Smith. 1982. Techniques of trend analysis for monthly water quality data. *Water Resources Research* 18, no. 1: 107–121. https://doi.org/10.1029/WR018i001p00107

Jasechko, S., H. Seybold, D. Perrone, et al. 2024. Rapid groundwater decline and some cases of recovery in aquifers globally. *Nature* 625, no. 7996: 715–721. https://doi.org/10.1038/s41586-023-06879-8

Kendall, M.G. 1975. *Rank Correlation Methods*. 4th ed. London: Charles Griffin.

Kish, L. 1965. *Survey Sampling*. New York: John Wiley & Sons.

Konikow, L.F. 2015. Long-term groundwater depletion in the United States. *Groundwater* 53, no. 1: 2–9. https://doi.org/10.1111/gwat.12306

Mann, H.B. 1945. Nonparametric tests against trend. *Econometrica* 13, no. 3: 245–259. https://doi.org/10.2307/1907187

Nolte, A., S. Bender, J. Hartmann, S. Baltruschat, and N. Moosdorf. 2026. Coastal groundwater-level trends reveal global susceptibility to seawater intrusion. *Nature Water* 4, no. 5: 673–682. https://doi.org/10.1038/s44221-026-00619-8

Ohmer, M., T. Liesch, and N. Goldscheider. 2019. On the optimal spatial design for groundwater level monitoring networks. *Water Resources Research* 55, no. 11: 9454–9473. https://doi.org/10.1029/2019WR025728

Rodell, M., J.S. Famiglietti, D.N. Wiese, et al. 2018. Emerging trends in global freshwater availability. *Nature* 557, no. 7707: 651–659. https://doi.org/10.1038/s41586-018-0123-1

Sen, P.K. 1968. Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association* 63, no. 324: 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

Stahl, M.O., T. Mar, and Y. Jameel. 2024. Seasonal groundwater level dynamics in unconfined aquifers across the United States. *Groundwater* 62, no. 6: 876–888. https://doi.org/10.1111/gwat.13422

Stets, E.G., A.A. Archer, J.R. Degnan, et al. 2025. The National integrated water availability assessment, water years 2010–20. U.S. Geological Survey Professional Paper 1894-A. https://doi.org/10.3133/pp1894A

Theil, H. 1950. A rank-invariant method of linear and polynomial regression analysis. *Proceedings of the Koninklijke Nederlandse Akademie van Wetenschappen* 53: 386–492.

Toran, L., T. Baraza, K. Çallı, et al. 2026. Long-term discharge records from United States and European karst springs: Data limitations in trend analysis. *Groundwater* 64, no. 5: 634–649. https://doi.org/10.1111/gwat.70093

U.S. Geological Survey. 2016. National Water Information System—web services. https://doi.org/10.5066/F7P55KJN (accessed October 2026)

Yue, S., P. Pilon, B. Phinney, and G. Cavadias. 2002. The influence of autocorrelation on the ability to detect trend in hydrological series. *Hydrological Processes* 16, no. 9: 1807–1829. https://doi.org/10.1002/hyp.1095
