# DT-RLR: Drift-Triggered Robust Local Reference

Reproducibility package for the manuscript:

**Drift-Triggered Adaptive Reference Calibration for Open-Set Electrical Fault Diagnosis under Temporal Drift**

Author: Alexandre Giacomelli Leal  
Affiliation: Independent Researcher, Brazil

## Overview

This repository contains the experimental scripts, derived result tables, and principal figures used to evaluate Drift-Triggered Robust Local Reference (DT-RLR). DT-RLR is an adaptive local-reference strategy for open-set electrical fault diagnosis under temporal drift.

The method maintains:
- a slow robust diagnostic reference;
- a fast robust recent-context reference;
- commissioning-derived drift thresholds;
- a persistent ADAPT/FREEZE trigger;
- a fast-context safe-update gate to reduce fault contamination.

The package also includes:
- Static RLR and continuously adaptive baselines;
- the original fixed safe-update gate;
- dual-timescale safe adaptation;
- persistence and drift-severity sensitivity experiments;
- conformal entropy operating-point analysis;
- a commissioning-calibrated Page-Hinkley trigger baseline;
- paired bootstrap and Wilcoxon comparisons.

## Public benchmark

The raw benchmark is **not redistributed in this repository**.

Download the public dataset from Zenodo:

**A Real Controlled-Fault Benchmark for Power-Anomaly Detection and Identification on IoT Edge Nodes (Three-Station 1 Hz Dataset)**  
DOI: https://doi.org/10.5281/zenodo.20565892

Place the downloaded archive in:

```text
data/ZENODO_UPLOAD_v1.0.zip
```

The expected archive contains:

```text
controlled-fault-benchmark-v1.0/
  station1_normal_features.csv
  station2_medium_features.csv
  station3_high_features.csv
```

## Python environment

Tested with:

- Python 3.11+
- NumPy 2.3.5
- pandas 2.2.3
- scikit-learn 1.8.0
- SciPy 1.17.0
- Matplotlib 3.10.8

Install:

```bash
python -m pip install -r requirements.txt
```

## Reproducibility and paths

The scripts are portable. By default they infer the repository root from their own location and read the benchmark from:

```text
data/ZENODO_UPLOAD_v1.0.zip
```

Two optional environment variables are available:

```text
DT_RLR_ROOT       # override the repository working root
DT_RLR_DATA_ZIP   # override the benchmark ZIP path
```

The numerical outputs used in the manuscript are already preserved under `results/`, so the exact reported summaries can be audited without rerunning the full benchmark.

## Recommended experiment order

### Step 1 — Static, naive adaptive, and fixed safe adaptive RLR

```bash
python scripts/paper2_step1_adaptive_rlr.py
```

This generates the reference event representations and entropy thresholds used by downstream experiments.

### Step 2 — Dual-timescale safe adaptation

```bash
python scripts/paper2_step2_dual_fast.py
```

### Step 3 — Electrically interpretable controlled perturbations

Run once per held-out target station:

```bash
python scripts/paper2_step3_station.py --station S1
python scripts/paper2_step3_station.py --station S2
python scripts/paper2_step3_station.py --station S3
```

### Step 4 — DT-RLR drift-triggered adaptation

```bash
python scripts/paper2_step4_triggered_rlr.py --station S1
python scripts/paper2_step4_triggered_rlr.py --station S2
python scripts/paper2_step4_triggered_rlr.py --station S3
```

### Step 5 — Persistence and drift-severity sensitivity

```bash
python scripts/paper2_step5_stress.py --station S1
python scripts/paper2_step5_stress.py --station S2
python scripts/paper2_step5_stress.py --station S3
```

### Step 6 — Conformal operating-point analysis

```bash
python scripts/paper2_step6_conformal_multiscore.py --station S1
python scripts/paper2_step6_conformal_multiscore.py --station S2
python scripts/paper2_step6_conformal_multiscore.py --station S3
```

### Major-revision audit and paired statistics

```bash
python scripts/paper2_major_revision_audit.py
```

### Page-Hinkley drift-trigger baseline

```bash
python scripts/paper2_page_hinkley_station.py --station S1
python scripts/paper2_page_hinkley_station.py --station S2
python scripts/paper2_page_hinkley_station.py --station S3
```

## Main frozen DT-RLR configuration

- Block duration: 1 hour
- Trigger persistence: `P = 3` blocks
- Safe fast-context gate: `z_gate = 4`
- Slow-reference update factor: `eta = 0.20`

Commissioning-derived thresholds:

```text
gamma_on  = max(Q0.95(D_ref), 0.50)
gamma_off = max(Q0.75(D_ref), 0.70 * gamma_on)
```

If `gamma_off >= gamma_on`, the implementation applies:

```text
gamma_off = 0.85 * gamma_on
```

All threshold statistics are derived exclusively from Day-1 commissioning data.

## Dataset split

For each target station:

- Day 1: unlabeled target commissioning/reference only.
- Days 2–7: source development or causal target evaluation.
- Target labels are used only for retrospective performance evaluation and contamination auditing.

The open-set benchmark combines:
- 3 leave-one-station-out target folds; and
- 13 leave-one-fault-type-out unknown-fault folds,

for 39 matched scenarios per configuration.

## Main input variables

The 13 features are:

```text
voltage_V
current_A
power_W
frequency_Hz
power_factor
PE_m3
PE_m4
PE_m5
WPE
MSPE
CV
delta_PF
SII
```

## Controlled temporal perturbations

The principal controlled perturbations are:

- Voltage ramp: voltage only
- Load ramp: current and active power
- Sensor gain: voltage, current, and active power
- Frequency offset: frequency only

They are intended as controlled, electrically interpretable perturbations and are **not** claimed to be full electromechanical simulations of naturally occurring ageing or long-term drift.

## Key files

`results/paper2_master_results.csv`  
Consolidated results from the main experiment stages.

`results/threshold_audit.csv`  
Commissioning-derived trigger thresholds by station.

`results/paired_statistics.csv`  
Paired bootstrap and Wilcoxon comparisons across matched open-set scenarios.

`results/page_hinkley_vs_dt_stats.csv`  
Paired DT-RLR versus Page-Hinkley trigger comparison.

`figures/paper2_major_revision_temporal_trigger.png`  
Representative temporal trace of the DT-RLR trigger.

## Interpretation

The experiments do **not** support a claim that DT-RLR universally outperforms static or continuously adaptive calibration in H-score. The supported conclusion is narrower:

> Adapt the reference when persistent drift indicates obsolescence, and update it only with observations compatible with the recent local context.

The method is intended to balance drift tracking, diagnostic preservation, and reference-contamination protection.

## Reuse

Software source code is released under the MIT License.

Derived tables and figures produced by this repository are released under CC BY 4.0, subject to citation of the associated manuscript and this archived repository.

The underlying public benchmark remains governed by the terms specified by its original Zenodo deposit.

## Citation

Please cite the archived Zenodo record for this repository once a DOI is assigned. A `CITATION.cff` file is included and can be updated with the Zenodo DOI after publication of the repository release.
