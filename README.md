# Multipath Background Learning for Target Detection and Tracking

This repository contains MATLAB code for **multipath background learning, target detection, and track-before-detect (TrBD) processing in time-varying underwater multipath environments**.

The repository contains two main components:

1. **Multipath Background Learning and SLRT Detection**
   Learning and tracking of the time-varying multipath background, estimation of model hyperparameters, and target detection using the Sequential Likelihood Ratio Test (SLRT).

2. **Background-Aware Track-Before-Detect (TrBD)**
   A particle-based target tracking framework that operates directly on the measurement data while accounting for the time-varying multipath background.


# Dataset

The measurement datasets required to reproduce the experiments are provided separately as **GitHub Release assets** because of their size.

## Download

1. Go to the **Releases** page of this repository:

   https://github.com/ASHKoul/Bcg_SLRT/releases

2. Download the dataset archive associated with the experiment.

3. Extract the downloaded archive before running the MATLAB scripts.

---

## Extracting the Dataset

### Windows

Right-click the downloaded `.zip` file and select:

```text
Extract All...
```

Place the extracted data in the working directory used by the corresponding MATLAB script.

### macOS / Linux

The archive can be extracted graphically or from the terminal:

```bash
unzip dataset.zip
```

---

# Requirements

The code is written in **MATLAB**.

Make sure that:

* MATLAB can access all files in the repository.
* The downloaded BELLHOP datasets have been extracted.
* The required `.mat`, configuration, and waveform files are available in the working directory.
* All required source folders have been added to the MATLAB path.

For example:

```matlab
addpath(genpath(pwd));
```

---

# 1. Multipath Background Learning and SLRT Detection

The first part of the repository evaluates multipath-background learning and target detection using the SLRT.

## Step 1 — Compute the SLRT Threshold

Run:

```matlab
Threshold_computation_for_SLRT_test
```

This computes the threshold required for the SLRT detection experiment.

---

## Step 2 — Run the SLRT Detection Experiment

Run:

```matlab
SLRT_for_bellhop_data
```

This script evaluates the detector using the BELLHOP-generated measurements.

The resulting analysis includes quantities such as:

* probability of detection,
* detection-delay-related performance,
* multipath-background adaptation behavior,
* and the corresponding performance plots.

---

## Step 3 — Estimate Multipath Hyperparameters

To estimate the parameters governing the multipath-background model, run:

```matlab
CIR_bellhop_basis_3param_MLE_all_SNR_bayesian_optimization
```

The script performs hyperparameter estimation using Bayesian optimization and evaluates the corresponding negative log-likelihood (NLL).

The estimated parameters characterize the evolution and uncertainty of the time-varying multipath background.

---

# 2. Background-Aware Track-Before-Detect (TrBD)

The `TrBD` folder contains the target-tracking implementation.

Unlike conventional detect-before-track processing, the TrBD framework uses the measurement information directly to recursively estimate:

* the probability that the target exists, and
* the target kinematic state.

The time-varying multipath background is recursively tracked and incorporated into the measurement model used by the target tracker.

---


## Running the TrBD Code

### Step 1 — Enter the TrBD Folder

In MATLAB, change the current working directory to:

```text
TrBD
```

or run:

```matlab
cd TrBD
```

It is recommended to add the folder and all of its subfolders to the MATLAB path:

```matlab
addpath(genpath(pwd));
```

---


## Background Tracking

Background tracking can be enabled using:

```matlab
s.bg_tracking = true;
```

The background model recursively estimates the time-varying multipath component before evaluating the target likelihood.

The relevant background-processing routines include:

```text
background_tracking.m
mult_background_processing.m
generate_basis_cov_matrix.m
power_calibrations.m
```

---

## Step 4 — Run the TrBD Experiment

After configuring the experiment, run:

```matlab
main
```

The script performs the complete simulation chain:

1. Loads the BELLHOP multipath measurements.
2. Loads the estimated multipath hyperparameters.
3. Generates the time-varying synthetic multipath background.
4. Generates the target trajectory and target return.
5. Combines the target signal with the background and measurement noise.
6. Configures the tracking filter.
7. Runs the background-aware TrBD filter for each Monte Carlo realization.
8. Stores the resulting target-state estimates and target-existence probabilities.

---

# TrBD Outputs

After `main.m` completes, the principal results are stored in:

```matlab
trbd_filter_out
```

The estimated target states are available through:

```matlab
trbd_filter_out.x
```

and the posterior target-existence probabilities through:

```matlab
trbd_filter_out.q
```

For each Monte Carlo realization,

```matlab
trbd_filter_out.q{imc}
```

contains the evolution of the posterior target-existence probability.

Similarly,

```matlab
trbd_filter_out.x{imc}
```

contains the corresponding target-state estimates.

These outputs can subsequently be used to compute tracking-performance quantities such as target-confirmation probability, confirmation time, and localization error.

---

## Plotting

Diagnostic plotting is controlled in `get_settings.m` using:

```matlab
s.do_plot = true;
```

For large Monte Carlo experiments, plotting can be disabled using:

```matlab
s.do_plot = false;
```

to reduce computational overhead.

---

# Recommended Workflow

For reproducing the complete set of experiments, the suggested order is:

```text
1. Download and extract the dataset
                |
                v
2. Estimate multipath hyperparameters
                |
                v
3. Compute the SLRT threshold
                |
                v
4. Run the SLRT detection experiment
                |
                v
5. Configure the TrBD experiment
                |
                v
6. Run TrBD/main.m
                |
                v
7. Evaluate detection and tracking performance
```

---

# Notes

* Large measurement datasets are distributed through **GitHub Releases** rather than being stored directly in the repository.
* Check the file paths in the main MATLAB scripts if the dataset is stored outside the repository directory.
* The scripts assume that the required waveform, configuration, BELLHOP data, and estimated hyperparameter files are available on the MATLAB path or in the working directory.
* Monte Carlo experiments can be computationally intensive. The number of realizations and particles can be reduced in the configuration files for preliminary testing.

---

# License

This project is distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.
