# A Scalable Risk-Scoring System Based on Consumer-Grade Wearables for Inpatients With COVID-19: Statistical Analysis and Model Development

Simon Föll (1*); Adrian Lison (1*); Martin Maritsch (1*); Karsten Klingberg (2*); Vera Lehmann (3);
	Thomas Züger (1,3,4); David Srivastava (2); Sabrina Jegerlehner (2); Stefan Feuerriegel (1,5);
	Elgar Fleisch (1,6); Aristomenis Exadaktylos (2); Felix Wortmann (1,6)

(1) Department of Management, Technology, and Economics, ETH Zürich, Zürich, CH\
(2) Department of Emergency Medicine, Inselspital, Bern, University Hospital, University of Bern, Bern, CH\
(3) Department of Diabetes, Endocrinology, Nutritional Medicine and Metabolism, Inselspital, Bern, University Hospital, University of Bern, Bern, CH\
(4) Department of Endocrinology, Diabetes and Metabolic Diseases, Kantonsspital Olten, Olten, CH\
(5) Institute of AI in Management, LMU Munich, Munich, DE\
(6) Institute of Technology Management, University of St. Gallen, St. Gallen, CH\
(*) These authors contributed equally to this work.

## Abstract
**Background**: To provide effective care for inpatients with COVID-19, clinical practitioners need systems that monitor patient health and subsequently allow for risk scoring. Existing approaches for risk scoring in patients with COVID-19 focus primarily on intensive care units (ICUs) with specialized medical measurement devices but not on hospital general wards.

**Objective**: In this paper, we aim to develop a risk score for inpatients with COVID-19 in general wards based on consumer-grade wearables (smartwatches).

**Methods**: Patients wore consumer-grade wearables to record physiological measurements, such as the heart rate (HR), heart rate variability (HRV), and respiration frequency (RF). Based on Bayesian survival analysis, we validated the association between these measurements and patient outcomes (ie, discharge or ICU admission). To build our risk score, we generated a low-dimensional representation of the physiological features. Subsequently, a pooled ordinal regression with time-dependent covariates inferred the probability of either hospital discharge or ICU admission. We evaluated the predictive performance of our developed system for risk scoring in a single-center, prospective study based on 40 inpatients with COVID-19 in a general ward of a tertiary referral center in Switzerland.

**Results**: First, Bayesian survival analysis showed that physiological measurements from consumer-grade wearables are significantly associated with patient outcomes (ie, discharge or ICU admission). Second, our risk score achieved a time-dependent area under the receiver operating characteristic curve (AUROC) of 0.73-0.90 based on leave-one-subject-out cross-validation.

**Conclusions**: Our results demonstrate the effectiveness of consumer-grade wearables for risk scoring in inpatients with COVID-19. Due to their low cost and ease of use, consumer-grade wearables could enable a scalable monitoring system.

**Trial Registration**: Clinicaltrials.gov NCT04357834; https://www.clinicaltrials.gov/ct2/show/NCT04357834

## Contents of this Repository
This repository contains the code used for the statistical analysis and modeling in this study:
- [Main analysis](analysis.Rmd)
- [Seasonality analysis](seasonality.Rmd)
- [MCMC diagnostics](diagnostics.Rmd)
- [Model robustness checks](robustness.Rmd)
- [Preprocessing](preprocessing)
- [Utility functions](utils)

The Bayesian survival models were implemented using the package brms, which interfaces to stan, see [here](https://github.com/paul-buerkner/brms) and [here](https://mc-stan.org/). The utility functions for performance evaluation using the time-dependent AUROC are based on code by Bansal and Heagerty, see [here](https://github.com/aasthaa/meanrankROC_package).

## Citation
If using parts of this code in your work, please consider citing the corresponding [paper](https://doi.org/10.2196/35717):

	Simon Föll; Adrian Lison; Martin Maritsch; Karsten Klingberg; Vera Lehmann;
	Thomas Züger; David Srivastava; Sabrina Jegerlehner; Stefan Feuerriegel;
	Elgar Fleisch; Aristomenis Exadaktylos; Felix Wortmann.
	A Scalable Risk-Scoring System Based on Consumer-Grade Wearables for Inpatients
	With COVID-19: Statistical Analysis and Model Development.
	JMIR Form Res 2022; 6(6):e35717
	DOI: 10.2196/35717
