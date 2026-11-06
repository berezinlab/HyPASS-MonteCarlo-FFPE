# HyPASS-MonteCarlo-FFPE

# HyPASS: Monte Carlo Modeling of SWIR Photon Transport in FFPE Tissue

## Overview

This repository contains MATLAB Monte Carlo simulations supporting the development of HyPASS (Hyperspectral Pathology Scanning System), an experimental imaging platform for direct, stain-free hyperspectral imaging of intact formalin-fixed, paraffin-embedded (FFPE) tissue blocks.

The simulations investigate how tissue thickness and assumed optical scattering properties influence detector-collected diffuse short-wave infrared (SWIR) reflectance.

## Optical Model

The model represents an exposed surface of paraffin-infiltrated tissue supported by an underlying layer of embedding paraffin. There is no overlying paraffin layer.

Photon transport incorporates absorption, scattering, Fresnel reflection, refraction, and detector collection geometry.

Two scattering scenarios are evaluated:

1. Literature-based reference coefficients used by Zhang et al. (2016), originally adopted from published rat-skin optical properties.
2. A higher-scattering sensitivity scenario based on an assumed wavelength-dependent power-law model.

Neither scenario represents experimentally measured optical properties of FFPE tissue.

## Simulation Parameters

| Parameter | Value |
|---|---|
| Wavelengths | 900, 1300, 1500, 1600 nm |
| Exposed tissue thickness | 5–5000 µm |
| Incident photons per condition | 500,000 |
| Maximum photon interactions | 15,000 |
| Tissue refractive index | 1.44 (assumed) |
| Paraffin refractive index | 1.46 (assumed) |
| Tissue scattering anisotropy | 0.90 |
| Detector numerical aperture | 0.20 (assumed) |
| Underlying paraffin thickness | 10 mm |

## Repository Contents

- `HyPASS_FFPE_Thickness_2016_Comparison.m` — MATLAB simulation and figure generation.
- `HyPASS_Wavelength_Optical_Properties.csv` — wavelength-dependent optical input data.
- `HyPASS_Objective1_Optical_Inputs.csv` — optical coefficients used in the simulations.
- `HyPASS_Objective1_Thickness_Comparison.csv` — numerical Monte Carlo outputs.
- `HyPASS_Objective1_Thickness.pdf` — comparison of the two scattering scenarios.
- `HyPASS_Objective1_Settings.mat` — saved model settings.

## Reproducibility

Place the MATLAB script and wavelength-dependent optical-properties CSV in the same directory. Run the script using MATLAB R2020b or later.

The simulation uses a fixed random seed (`rng(61,'twister')`) and exports numerical results and figures. Depending on the script version, outputs are written to its directory or a designated temporary results folder.

## Limitations

This preliminary model assumes plane-parallel, laterally infinite tissue and paraffin layers. It does not explicitly represent tissue heterogeneity, irregular tissue–paraffin boundaries, wavelength-dependent refractive-index dispersion, or microscope spatial resolution.

The FFPE and paraffin optical properties are provisional. The results provide a sensitivity analysis and physical rationale for block-face imaging, not experimental validation of FFPE optical behavior.

## Reference

Zhang H, Salo D, Kim DM, Komarov S, Tai YC, Berezin MY. Penetration depth of photons in biological tissues from hyperspectral imaging in shortwave infrared in transmission and reflection geometries. *Journal of Biomedical Optics*. 2016;21(12):126006.

DOI: https://doi.org/10.1117/1.JBO.21.12.126006

## Project

HyPASS — Hyperspectral Pathology Scanning System

Berezin Laboratory, Washington University in St. Louis.

## Availability

Research code is maintained for scientific reproducibility. Repository access and any public distribution are subject to institutional intellectual-property review.
