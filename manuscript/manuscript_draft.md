# How sparse is too sparse? A closed-form law for the reliability of groundwater-level trends

## Authors

Zijie Deng

Jiangxi Provincial Key Laboratory of Water Resources Allocation and Efficient Utilization, Jiangxi University of Water Resources and Electric Power, No. 289 Tianxiang Avenue, Nanchang 330099, China

Corresponding author: Zijie Deng, 3335196397vip@gmail.com (ORCID 0009-0002-4824-9764).

**Running head:** Network density and trend reliability

---

## Abstract

Long-term water-level trends from monitoring networks are widely used to assess aquifer depletion, allocate pumping, and evaluate regulatory compliance, yet many networks are losing wells to budget cuts and attrition. This study examines how reliably a reduced monitoring network reproduces the regional trend direction of the full network. Using 8.8 million quality-controlled U.S. Geological Survey depth-to-water measurements collected between 1950 and 2026, 29,717 monitoring wells with at least 20 years of record were identified across the conterminous United States. Regional anomaly hydrographs were constructed over a common 51-year period (1975–2025) for 21 aquifer strata and the national network. Each network was virtually thinned 1,000 times at retention fractions from 90% down to 1%, removing wells either uniformly at random or in 0.5° spatial clusters. Across all 216 thinning trials, the probability that a thinned network reverses the full-network trend sign follows a finite-population sampling law controlled by the retained fraction and the signal-to-noise ratio (SNR, the magnitude of the regional trend divided by the standard deviation of individual well trends). Inverting this relationship at a 10% sign-reversal risk gives a closed-form minimum network size, n_min = Nz²/(z²+(c·SNR)²N), with c = 1.27 (95% confidence interval: 1.04–1.55), which reproduces simulated minimum network sizes within a factor of two for 17 of 19 resolvable networks. Strong-signal aquifers preserve their regional trend sign with as few as 11 wells (Coastal Lowlands) or 81 of 8,160 wells (High Plains)—a mathematical threshold for the aquifer-wide sign rather than a target for local management—whereas three systems with weak, spatially offsetting trends (Basin and Range, surficial, Floridan) require 88–94% of their wells. Spatially clustered well losses are much more damaging than random losses: the national trend tolerates random removal down to 297 wells (1%), but 7,726 wells (30%) retained in spatial clusters carry the same reversal risk as ~180–280 random wells.

## Impact Statement

A closed-form sampling law predicts the minimum number of monitoring wells needed to detect reliable aquifer-scale water-level trends.

## Keywords

groundwater monitoring; water-level trends; network design; trend detection; depletion; Theil–Sen

---

## 1. Introduction

Water-level records from observation wells are the primary evidence used to estimate groundwater depletion (Konikow 2015; Jasechko et al. 2024), validate satellite gravimetry (Rodell et al. 2018), and support national water-availability assessments (Stets et al. 2025). Regulatory frameworks such as California's Sustainable Groundwater Management Act also rely on monitored water-level trends to evaluate whether aquifers are being pumped sustainably.

In practice, long-term monitoring networks are costly to maintain and rarely match management needs (Taylor and Alley 2001). The U.S. National Ground-Water Monitoring Network still has substantial geographic gaps (Wunsch and Schreiber 2024), and state or local networks periodically lose wells when funding ends or wells fail. Most prior network-design studies have focused on choosing well locations to minimize spatial interpolation error (Loaiciga et al. 1992; Ohmer et al. 2019) or on classifying seasonal water-level cycles (Stahl et al. 2024). Less attention has been paid to how network size affects the *trend conclusion* itself: if a monitoring network is cut, how often does the remaining subset reverse the regional trend direction or lose statistical significance? This question is practical because regional storage estimates often disagree across methods (Richey et al. 2015), sparse coastal networks detect significant trends in only about one-fifth of grid cells (Nolte et al. 2026), and data-driven mapping can easily outrun the underlying field observations (Bradbury 2025).

This study evaluates how monitoring density controls regional trend reliability across the conterminous United States (CONUS) using 8.8 million discrete water-level measurements from 29,717 long-term wells (1950–2026), analyzed over a common 51-year period (1975–2025). Three questions are addressed: (1) How often does a reduced network reverse the full-network trend sign or lose statistical significance under random versus spatially clustered well removal? (2) Can sign-reversal probability be described by a closed-form finite-population sampling law that predicts the minimum reliable network size (n_min) from existing records? (3) When wells must be cut, does keeping the longest records or spatially spread wells outperform random retention?

## 2. Data and Methods

### 2.1 Groundwater-level data and analysis time windows

Discrete depth-to-water measurements (USGS parameter code 72019, depth to water level below land surface) collected between January 1950 and October 2026 were downloaded from the U.S. Geological Survey Water Data API (U.S. Geological Survey 2016), yielding 9,456,561 raw measurements from 843,302 sites. Quality-assurance filtering retained approved records within 0–1,524 m (0–5,000 ft) below land surface, removed measurements flagged for pumping or dry-well conditions (234,752 excluded), reduced same-day repeat measurements at a well to their daily median (8,814,831 well-days), and excluded sites outside CONUS (103,049 excluded). Each well was required to have at least 30 ambient measurements spanning at least 20 distinct calendar years. This screen retains 29,717 long-term wells (3.5% of all sites) across all 48 conterminous states and Washington, D.C.; 3,064 wells (10.3%) share coordinates with another well (nested piezometers or replacement wells) and were kept as separate series. Grouping wells by USGS principal aquifer code gives 21 strata with at least 100 eligible wells (Table 1): 20 principal aquifer systems plus 6,932 unassigned wells.

Because individual wells have different record spans within 1950–2026, two time windows are used. Individual well trends were calculated over each well's full available record within 1950–2026 (≥ 20 annual values). Regional composite anomaly hydrographs and network-thinning experiments were restricted to the common 51-year period from 1975 to 2025, when coverage is dense across all 21 strata and complete calendar years are available. Within 1975–2025, the median eligible well has observations in 29 of the 51 years, and 3.5% cover all 51 years.

### 2.2 Well-scale trend estimation

To remove seasonal sampling bias (cf. Hirsch et al. 1982), each well's record was deseasonalized by subtracting its long-term month-of-year median depth to water for months with at least two observations (sparser months were centered on the well's overall median), and the residuals were aggregated to annual medians. Long-term trends were estimated with the non-parametric Theil–Sen slope (Theil 1950; Sen 1968). Statistical significance (two-sided p < 0.05) was tested with the Mann–Kendall test (Mann 1945; Kendall 1975) using Yue–Pilon trend-free pre-whitening (Yue et al. 2002) where lag-1 serial correlation was significant at the 95% level (60.2% of wells). Positive depth-to-water slopes indicate falling water levels; 62.3% of the 29,717 wells have positive slopes and 70.1% are statistically significant (p < 0.05).

### 2.3 Network-thinning experiments

For each of the 21 strata and the CONUS network (22 series), each well's annual series was centered on its 1975–2025 median and averaged across wells in each calendar year to form the regional anomaly hydrograph; the full-network reference trend (β_net) is the 1975–2025 Theil–Sen slope of this composite series. This definition targets the network-average trend direction rather than rank concordance alone (cf. Fang et al. 2019). Alternative reference definitions—spatial-median aggregation, inverse-density weighting, balanced panels, or shifting the start year between 1970 and 1980—preserve the reference sign in 15 to 20 of the 22 series, with sign changes limited to weak-signal strata whose bootstrap 95% confidence intervals include or touch zero (Table 1; Supporting Information Section S5).

Network contraction was simulated at 11 retention fractions (f = 0.90 to 0.01). At each fraction, n = max(10, ⌊fN⌋) wells were drawn without replacement 1,000 times and the composite slope and Mann–Kendall p-value were recomputed over 1975–2025. A replicate was counted as a *sign flip* if its slope had the opposite sign to β_net, and as a *loss of significance* if the reference trend had p < 0.05 but the subsample had p ≥ 0.05 or flipped sign. In *uniform random thinning*, wells were drawn with equal probability. In *spatially clustered thinning*, wells were grouped into 0.5° × 0.5° blocks (excluding blocks with < 10 wells, leaving N = 25,754 in CONUS) and whole blocks were removed at random to represent district-scale cuts. The minimum reliable network size (n_min) is the smallest n with sign-reversal probability ≤ τ (τ = 0.10). Because subsampling stops at 10 wells, floor values are upper bounds; merging duplicate floor sizes in small strata leaves 216 distinct random-thinning trials. Finally, in the High Plains and Central Valley aquifers, 50 random deletion orders were compared against *longest-records-first* and *spatial-spread-first* retention.

### 2.4 Finite-population sampling law for sign-reversal risk

Physically, a thinned monitoring network is much less likely to reverse its regional trend direction when the aquifer-wide trend is large relative to well-to-well variability in slopes. Because random thinning draws n of N wells without replacement from a fixed network, this intuition can be formalized with classical finite-population sampling theory (Cochran 1977). Spatial correlation among wells is fixed in the population variance S² of the N individual well slopes and modifies sampling variance only under clustered removal. To first order, the slope of the averaged anomaly hydrograph behaves like the mean of the n retained well slopes, which has sampling variance S²(N − n)/(nN). Treating the subsample slope distribution around β_net as approximately normal yields a one-parameter sampling law for sign-reversal probability:

  P(flip) = Φ[−c · SNR · √(nN/(N − n))],   SNR = |β_net| / S,   (1)

where Φ is the standard normal cumulative distribution function, SNR = |β_net| / S is the network signal-to-noise ratio (the 1975–2025 reference trend magnitude divided by the standard deviation S of the individual well slopes in that stratum), and c is an empirical constant capturing the efficiency of the composite Theil–Sen slope relative to a simple mean of well slopes. Binomial maximum-likelihood fitting across the informative trials (those with observed flip rates strictly between 0 and 1) gives c = 1.27, with a 95% cluster-bootstrap confidence interval of 1.04–1.55 (leave-one-stratum-out range: 1.24–1.42; individual stratum fits c_i range from 0.55 to 4.05). Synthetic-population tests, estimator decompositions, and half-split validations confirm that c = 1.27 reflects the composite Theil–Sen estimator rather than slope skewness or temporal noise, and that out-of-sample coverage matches in-sample coverage (Supporting Information Sections S1 and S5). Because c = 1.27 is calibrated to the composite Theil–Sen procedure, SNR inputs must be computed with the deseasonalized Theil–Sen pipeline of Sections 2.2–2.3.

Solving for n at a maximum acceptable flip probability τ (with normal quantile z = Φ⁻¹(1 − τ); τ = 0.10, z = 1.282) gives the minimum reliable network size:

  n_min = N z² / (z² + (c·SNR)² N).   (2)

When the regional signal is strong ((c·SNR)²N ≫ z²), N cancels out and n_min ≈ (z/(c·SNR))² is independent of total network size. When the signal is weak ((c·SNR)²N ≲ z²), n_min approaches N because opposing local trends nearly cancel. Varying c across 1.04–1.55 changes predicted n_min by ~2.2× in strong-signal networks (e.g., 55–122 wells around the predicted 82 wells in the High Plains). For clustered removal, inverting Equation 1 maps the observed flip rate to its effective random sample size (n_eff).

### 2.5 Software

All analyses were conducted in Python 3.9 (pandas, NumPy, SciPy, Matplotlib).

## 3. Results

### 3.1 National monitoring network and reference trends

Eligible wells are concentrated in the High Plains, California, Florida, the Mississippi embayment, and the glaciated north (Fig. 1a). Water levels fall at 18,513 wells (62.3%), rise at 11,134 (37.5%), and are unchanged at 70 (0.2%); the median well slope is +0.010 m yr⁻¹ (interquartile range: −0.010 to +0.083; extremes: −9.8 to +5.4 m yr⁻¹). Across the 21 strata, 1975–2025 reference slopes span −0.170 to +0.074 m yr⁻¹, and well-bootstrap 95% confidence intervals exclude zero only where regional signals are strong (Table 1).

### 3.2 Sign-reversal probability under random and clustered thinning

Under random thinning, sign-reversal curves differ sharply across aquifers (Fig. 2a). High-SNR systems keep their reference trend sign down to a few percent of their wells (Coastal Lowlands, n_min = 11; Mississippi River Valley alluvium, 12), whereas the Basin and Range aquifer (SNR = 0.006, reference slope +0.002 m yr⁻¹) reverses sign in up to one-third of replicates at every size tested.

On the CONUS network, random and spatially clustered removals behave very differently (Fig. 2b). Random thinning produces zero sign flips down to 5% retention (1,485 wells) and a 0.075 flip rate at the 1% floor (297 of 29,717 wells). Under 0.5° clustered removal, the flip rate reaches 0.10 at 30% retention (7,726 wells) and drops below 10% only at 50% retention (12,877 wells). Inverting the sampling law shows that between 10% and 30% clustered retention, 7,726 clustered wells carry the flip risk of only n_eff ≈ 180–280 random wells, matching the number of surviving 0.5° blocks. A two-stage cluster sampling calculation (Kish 1965) from archive block statistics (intracluster correlation ρ = 0.32, mean block size 41.5, coefficient of variation 1.04) gives a design effect of DEFF ≈ 29; setting n_eff = n / 29 with no free parameters predicts observed clustered flip rates within ~1 percentage point across 20–30% retention (Supporting Information Section S2). Using 0.25° or 1.0° blocks also yields CONUS flip rates above 10% at 30% retention (0.11–0.15). At the aquifer scale, 30% clustered retention keeps the High Plains stable (0.002 flip rate) but raises flip rates to 0.31 in the Coastal Lowlands and 0.17 in the Central Valley.

Across all networks, statistical significance (p < 0.05) is lost at larger sample sizes than trend direction (at 1% in CONUS, 18.5% of replicates lose significance while 7.5% flip sign; at 30% in glaciated aquifers, 36% lose significance and none flips). Losing significance is driven by increased interannual noise in the averaged hydrograph as fewer wells are averaged each year, rather than spatial slope variance (a threshold-shifted companion form of Equation 1 captures only 31% of trials within 95% binomial bands; Supporting Information Section S3).

### 3.3 Minimum reliable network sizes across aquifers

Simulated minimum reliable sizes (n_min for ≤ 10% flip risk; Fig. 1b and Table 1) divide the 21 strata and the national network into three regimes. In the strong-signal regime, six strata and the CONUS network stay below 10% flip risk all the way down to the grid floor of 1% of N or 10 wells (such as Coastal Lowlands at 11 wells, Mississippi alluvium at 12, High Plains at 81, and CONUS at 297, all upper bounds; Table 1), and three additional coastal and glaciated strata resolve between 2% and 5% retention (25–71 wells). In the intermediate regime, nine heterogeneous aquifers require 20–80% of their wells, ranging from ~20% in the Biscayne aquifer to 70–80% in the Rio Grande (519 of 742 wells), Mississippi embayment (333 of 477), Cambrian–Ordovician (71 of 102), and Central Valley (146 of 183, where severe local declines up to 1.4 m yr⁻¹ are partly offset by managed-recovery areas). Notably, both the Rio Grande and Roswell Basin reference series are statistically significant in time (Mann–Kendall p < 0.05), yet their spatial bootstrap confidence intervals across wells include zero (Table 1), showing that temporal significance does not guarantee spatial sampling stability. Finally, in the weak-signal regime, three aquifers with spatially offsetting trends (Basin and Range, surficial, and Floridan) never reach ≤ 10% flip risk on the standard grid; Equation 2 predicts required retentions of 94%, 91%, and 88%, and dense simulations at 92–98% retention confirm the first two within ~2 percentage points while placing the Floridan threshold near 93%.

Across the 19 uncensored networks, Equation 2 reproduces simulated n_min within a factor of two for 17 networks (89%; median log₁₀ error 0.15 dex, comparable to the 1.5× median grid step). The two larger misses occur in transition aquifers in opposite directions (glaciated aquifers: 213 predicted vs. 71 simulated; Roswell Basin: 17 predicted vs. 77 simulated).

### 3.4 Collapse across networks and the design chart

Plotting all 216 random-thinning trials against the normalized sample size x = c · SNR · √(nN/(N − n)) collapses the observations onto the standard normal curve Φ(−x) with no fitted shape parameter (Fig. 3a). In the low-risk regime (< 5% flip rate), 74% of trials fall within their 95% binomial Monte Carlo intervals. Agreement is looser in the transition zone (23% coverage), where slope skewness causes the composite Theil–Sen estimator to depart from a simple mean (Supporting Information Section S5). In the strongly skewed Rio Grande and Floridan systems, seven deep-thinning trials reach flip rates of 0.51–0.59; because Φ(−x) cannot exceed 0.50, none of these seven trials falls within the binomial bands, though all lie inside the low-SNR zone already flagged as unreliable. Plotting required retention contours n_min / N across (SNR, N) gives a direct design chart (Fig. 3b).

### 3.5 Deterministic retention strategies: slope bias versus sign risk

In the High Plains (N = 8,160) and Central Valley (N = 183), deterministic retention rules mainly alter slope bias rather than sign-reversal risk (Figure S1; Supporting Information Section S4). In the High Plains, no rule flips the depletion sign, but at 10% retention (816 wells) keeping the longest records first produces 5.4× the slope bias of random retention (0.017 vs. 0.003 m yr⁻¹) because the oldest wells cluster in early pumping districts. In the Central Valley, keeping spatially spread wells preserves the reference sign down to 8% retention (14 wells) but inflates the slope to 3.1× the reference value.

## 4. Discussion

**Designing and consolidating monitoring networks.** These findings translate directly into practical guidance for sizing, laying out, and scheduling observation networks. To evaluate *network size*, an agency can compute β_net and S from its existing archive using the deseasonalized Theil–Sen procedure, form SNR = |β_net| / S, and read n_min from Figure 3b. Unlike geostatistical designs that minimize kriging variance across a water-level map (Loaiciga et al. 1992; Ohmer et al. 2019), this calculation gives the mathematical threshold needed to preserve the sign of the aquifer-wide composite trend; local management tasks such as tracking drawdown cones, enforcing water-rights compliance, or mapping hydraulic gradients require much denser monitoring. Where (c·SNR)²N ≫ z², a modest fixed number of wells suffices for the regional sign, whereas where (c·SNR)²N ≲ z², opposing local trends nearly cancel and sub-basin trends should be reported separately. For *spatial layout*, because well trends are correlated within 0.5° blocks (ρ = 0.32), retiring entire districts inflates sampling variance by a Kish design effect of ~29 nationally, meaning that unavoidable cuts must be spread across the network footprint. For *measurement frequency*, subsampling each well to one measurement per year leaves the median ratio of annual-visit SNR to full-archive SNR near unity (1.06, interquartile range 0.82–1.14), but reduces SNR enough in seven of the 21 strata to shift them into the high-retention regime ((c·SNR)²N ≲ z²), so agencies should check stratum-specific SNR before reducing field visits.

**Implications for regional syntheses.** Regional and global compilations (Konikow 2015; Rodell et al. 2018; Jasechko et al. 2024; Stets et al. 2025) inherit the sampling risk of their underlying networks: in moderate- to low-SNR aquifers, a network covering half the domain carries a 10–35% chance of reversing the regional trend sign. Just as trend estimates depend on the time window analyzed (Toran et al. 2026), they depend equally on which wells are sampled. Reporting archive SNR and implied sign-reversal risk alongside regional trend estimates and model calibrations (Wöhling et al. 2025) provides a simple check on network sparsity.

**Scope and limitations.** Several practical boundaries apply to this framework. Depth-to-water records (parameter 72019) were used to maximize national coverage, and co-located nested piezometers (10.3% of wells) were retained as separate series; splitting the three unstable strata by USGS aquifer-type code confirms that their low SNR reflects spatially opposing trends rather than mixing confined and unconfined wells, though local confined-aquifer studies should screen by depth and analyze hydraulic head. In addition, where regional trends are near zero, the reference sign can vary with the spatial aggregation rule, a condition flagged whenever the bootstrap confidence interval in Table 1 crosses zero. Preferential attrition—such as shallow wells going dry and dropping out of the network first—also biases the surviving sample in the same way as longest-record-first retention. Within its calibrated domain (N ≥ 100, n ≥ 10–15, τ = 0.10, using the procedure of Section 2.2), Equation 2 provides a reliable first-order benchmark, whereas direct Monte Carlo thinning is recommended in transition aquifers.

## 5. Conclusions

1. Across 29,717 long-term U.S. monitoring wells (1950–2026) evaluated over 1975–2025, the probability that a thinned network reverses the full-network trend direction follows a finite-population sampling law, P(flip) = Φ[−c · SNR · √(nN/(N − n))] with c = 1.27 (95% CI: 1.04–1.55). Its closed-form inverse, n_min = Nz²/(z²+(c·SNR)²N), reproduces simulated minimum network sizes within a factor of two for 17 of 19 uncensored networks. Statistical significance (p < 0.05) is lost at larger sample sizes than trend direction because thinning increases interannual hydrograph noise.

2. Minimum reliable network sizes for ≤ 10% sign-reversal risk fall into three regimes: six strong-signal strata and the national network remain stable down to the 1% (or 10-well) grid floor (such as 11 wells in the Coastal Lowlands and 81 of 8,160 in the High Plains, all upper bounds), nine heterogeneous strata require 20–80% of their wells (including 70–80% in the Rio Grande and Central Valley), and three weak-signal systems with spatially offsetting trends (Basin and Range, surficial, and Floridan) require 88–94% of their wells.

3. Spatially clustered well losses are far more damaging than random losses. The national trend tolerates uniform random thinning down to 297 of 29,717 wells (1%), whereas removing 0.5° spatial blocks exceeds 10% reversal risk at 7,726 wells (30% retention), equivalent to only ~180–280 random wells (Kish design effect ≈ 29). Well reductions should therefore be spread across the network rather than concentrated by district, and consolidation should avoid keeping only the longest records, which inflates bias.

## Acknowledgments

The author thanks the U.S. Geological Survey for maintaining the National Water Information System.

### Author's Note

The author has no conflicts of interest or financial disclosures to report.

## Data Availability Statement

Groundwater-level measurements were retrieved from the U.S. Geological Survey Water Data API (parameter code 72019), accessed October 2026. The complete analysis code and derived datasets are available at https://github.com/yuelangmanle/gwat_wellnet, and the synthetic validation testbed is described in the Supporting Information.

## References
Bradbury, K. 2025. Groundwater potential mapping: A misused and dubious concept. *Groundwater* 63, no. 5: 664–665. https://doi.org/10.1111/gwat.70013

Cochran, W.G. 1977. *Sampling Techniques*. 3rd ed. New York: John Wiley & Sons.

Fang, C., S. Sun, S. Jia, and Y. Li. 2019. Groundwater level analysis using regional Kendall test for trend with spatial autocorrelation. *Groundwater* 57, no. 2: 320–328. https://doi.org/10.1111/gwat.12800

Hirsch, R.M., J.R. Slack, and R.A. Smith. 1982. Techniques of trend analysis for monthly water quality data. *Water Resources Research* 18, no. 1: 107–121. https://doi.org/10.1029/WR018i001p00107

Jasechko, S., H. Seybold, D. Perrone, et al. 2024. Rapid groundwater decline and some cases of recovery in aquifers globally. *Nature* 625, no. 7996: 715–721. https://doi.org/10.1038/s41586-023-06879-8

Kendall, M.G. 1975. *Rank Correlation Methods*. 4th ed. London: Charles Griffin.

Kish, L. 1965. *Survey Sampling*. New York: John Wiley & Sons.

Konikow, L.F. 2015. Long-term groundwater depletion in the United States. *Groundwater* 53, no. 1: 2–9. https://doi.org/10.1111/gwat.12306

Loaiciga, H.A., R.J. Charbeneau, L.G. Everett, G.E. Fogg, B.F. Hobbs, and S. Rouhani. 1992. Review of ground-water quality monitoring network design. *Journal of Hydraulic Engineering* 118, no. 1: 11–37. https://doi.org/10.1061/(ASCE)0733-9429(1992)118:1(11)

Mann, H.B. 1945. Nonparametric tests against trend. *Econometrica* 13, no. 3: 245–259. https://doi.org/10.2307/1907187

Nolte, A., S. Bender, J. Hartmann, S. Baltruschat, and N. Moosdorf. 2026. Coastal groundwater-level trends reveal global susceptibility to seawater intrusion. *Nature Water* 4, no. 5: 673–682. https://doi.org/10.1038/s44221-026-00619-8

Ohmer, M., T. Liesch, and N. Goldscheider. 2019. On the optimal spatial design for groundwater level monitoring networks. *Water Resources Research* 55, no. 11: 9454–9473. https://doi.org/10.1029/2019WR025728

Richey, A.S., B.F. Thomas, M.-H. Lo, et al. 2015. Quantifying renewable groundwater stress with GRACE. *Water Resources Research* 51, no. 7: 5217–5238. https://doi.org/10.1002/2015WR017349

Rodell, M., J.S. Famiglietti, D.N. Wiese, et al. 2018. Emerging trends in global freshwater availability. *Nature* 557, no. 7707: 651–659. https://doi.org/10.1038/s41586-018-0123-1

Sen, P.K. 1968. Estimates of the regression coefficient based on Kendall's tau. *JASA* 63, no. 324: 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

Stahl, M.O., T. Mar, and Y. Jameel. 2024. Seasonal groundwater level dynamics in unconfined aquifers across the United States. *Groundwater* 62, no. 6: 876–888. https://doi.org/10.1111/gwat.13422

Stets, E.G., A.A. Archer, J.R. Degnan, et al. 2025. The National integrated water availability assessment, water years 2010–20. U.S. Geological Survey Professional Paper 1894-A. https://doi.org/10.3133/pp1894A

Taylor, C.J., and W.M. Alley. 2001. Ground-water-level monitoring and the importance of long-term water-level data. U.S. Geological Survey Circular 1217.

Theil, H. 1950. A rank-invariant method of linear and polynomial regression analysis. *Proc. K. Ned. Akad. Wet.* 53: 386–492.

Toran, L., T. Baraza, K. Çallı, et al. 2026. Long-term discharge records from United States and European karst springs: Data limitations in trend analysis. *Groundwater* 64, no. 5: 634–649. https://doi.org/10.1111/gwat.70093

U.S. Geological Survey. 2016. National Water Information System—web services. https://doi.org/10.5066/F7P55KJN (accessed October 2026)

Wunsch, D.R., and R.P. Schreiber. 2024. Calling all groundwater professionals: Support the National Groundwater Monitoring Network. *Groundwater* 62, no. 3: 328–329. https://doi.org/10.1111/gwat.13398

Wöhling, T., A. Delgadillo, M. Kraft, and A. Guthke. 2025. Comparing physics-based, conceptual and machine-learning models to predict groundwater levels by BMA. *Groundwater* 63, no. 4: 484–505. https://doi.org/10.1111/gwat.13487

Yue, S., P. Pilon, B. Phinney, and G. Cavadias. 2002. The influence of autocorrelation on the ability to detect trend in hydrological series. *Hydrological Processes* 16, no. 9: 1807–1829. https://doi.org/10.1002/hyp.1095
