# NeuralDataAnalysis_MATLAB

This repository collects MATLAB functions developed for published papers and for future research work, with the goal of maintaining a single, organized, and reusable codebase.

## Purpose

The main purpose of this repository is to provide a unique place where MATLAB functions can be stored, documented, and reused across different papers and research projects.  
It is intended to grow over time as new methods, utilities, and supporting code are developed.

## Contents

This repository will include:

- MATLAB functions used in published papers
- MATLAB functions developed for ongoing or future papers
- Reusable code modules for research activities

## Usage and License

The code in this repository is free to use without requesting permission.

However, if you use any function from this repository, please cite:

- the paper for which the code was originally developed, when specified in the header of the code, or
- this repository, if no specific paper is indicated in the code header

## Disclaimer

The code is provided as is, without any warranty or responsibility for the correctness of the results.

Please always check your results carefully and validate them independently before using them in research, publications, or other applications.

**Check your results twice.**

## Citation

If no specific related paper is indicated in the function header, please cite this repository.

## Notes

This repository is intended as a personal and evolving collection of research code.  
Documentation and organization may improve over time as new functions are added.

## Functions

- **`analyzeGazeData.m`** – Analyzes gaze position data to detect and characterize saccadic and fixation events based on velocity thresholds.

- **`computeDAMindex.m`** – Computes a data-availability/adjacency-missingness (DAM) index to quantify missing data in 2-D matrices, accounting for spatial clustering of missing elements.

- **`computePI.m`** – Computes the Preference Index (PI) as a measure of population neural response favorability (Moody et al., 1998, Hadjidimitrakis et al., 2011).

- **`computedPI.m`** – Computes the depth Preference Index (dPI) as a measure of directional selectivity in neural population responses (Moody et al., 1998, Hadjidimitrakis et al., 2011).

- **`computeMatchedFanoFactor.m`** – Computes the Fano Factor for two datasets with matched mean distributions through iterative resampling, a measure of noise characteristics in neural populations (Churchland et al., 2010).

- **`matchDistributions.m`** – Resamples two distributions to match their empirical probability distributions across percentiles.

- **`vectorCorrelation.m`** – Computes 2D vector correlation between two paired vector sets, providing correlation coefficient, scale factor, rotation/reflection angle and p-value (Hanson et al., 1992).
